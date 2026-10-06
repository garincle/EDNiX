#!/usr/bin/env bash
# Open the EDNiX interface running on a remote server (for example the GPU server) in the
# local browser, through an SSH tunnel. Runs keep going on the server when you close it.
#
#   bash GUI/open_remote.sh                 # host and folder from ednix_gui.conf
#   bash GUI/open_remote.sh gpu-07-01       # host given here, EDNiX folder from ednix_gui.conf
#
# Settings (ednix_gui.conf, same folder):
#   EDNIX_REMOTE_HOST=gpu-07-01                         ssh host (an alias of ~/.ssh/config works)
#   EDNIX_REMOTE_REPO=/home/cgarin/PycharmProjects/EDNiX EDNiX folder on that host
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF="$DIR/ednix_gui.conf"
if [ -f "$CONF" ]; then
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line#$'\xef\xbb\xbf'}"; line="${line%$'\r'}"; line="${line#"${line%%[![:space:]]*}"}"
    case "$line" in
      EDNIX_REMOTE_HOST=*|EDNIX_REMOTE_REPO=*)
        key="${line%%=*}"; val="${line#*=}"; val="${val%\"}"; val="${val#\"}"; val="${val%\'}"; val="${val#\'}"
        export "$key=$val";;
    esac
  done < "$CONF"
fi
HOST="${1:-$EDNIX_REMOTE_HOST}"
REMOTE_REPO="${EDNIX_REMOTE_REPO:-\$HOME/PycharmProjects/EDNiX}"
if [ -z "$HOST" ]; then
  echo "No server given. Set EDNIX_REMOTE_HOST in $CONF or run: bash GUI/open_remote.sh <ssh host>"
  exit 1
fi
command -v ssh >/dev/null || { echo "ssh is not installed on this computer."; exit 1; }

SOCK="$(mktemp -u "${TMPDIR:-/tmp}/ednix_ssh.XXXXXX")"
LOG="$(mktemp "${TMPDIR:-/tmp}/ednix_remote.XXXXXX")"
cleanup() { ssh -S "$SOCK" -O exit "$HOST" >/dev/null 2>&1; rm -f "$LOG"; }
trap cleanup EXIT INT TERM

echo "Connecting to $HOST ..."
# one shared SSH connection (ControlMaster): password / key asked once
ssh -M -S "$SOCK" -fnNT -o ControlPersist=no "$HOST" || { echo "Could not connect to $HOST."; exit 1; }
# start the interface on the server; it prints its address (with its access key)
ssh -S "$SOCK" -tt "$HOST" "bash $REMOTE_REPO/GUI/launch_EDNiX_GUI.sh --no-browser" > "$LOG" 2>&1 &
REMOTE_PID=$!
URL=""
for _ in $(seq 1 90); do
  URL="$(grep -o 'http://127\.0\.0\.1:[0-9]*/#token=[A-Za-z0-9_-]*' "$LOG" | head -1)"
  [ -n "$URL" ] && break
  kill -0 "$REMOTE_PID" 2>/dev/null || break
  sleep 1
done
if [ -z "$URL" ]; then
  echo "The interface did not start on $HOST:"; cat "$LOG"; exit 1
fi
RPORT="$(echo "$URL" | sed -E 's#http://127\.0\.0\.1:([0-9]+)/.*#\1#')"
TOKEN="${URL#*#token=}"
LPORT="$(python3 -c 'import socket;s=socket.socket();s.bind(("127.0.0.1",0));print(s.getsockname()[1])' 2>/dev/null || echo 8766)"
ssh -S "$SOCK" -O forward -L "$LPORT:127.0.0.1:$RPORT" "$HOST" >/dev/null || { echo "Could not open the tunnel."; exit 1; }
LOCAL_URL="http://127.0.0.1:$LPORT/#token=$TOKEN"
echo "EDNiX interface of $HOST: $LOCAL_URL"
if command -v xdg-open >/dev/null; then xdg-open "$LOCAL_URL" >/dev/null 2>&1
elif command -v open >/dev/null; then open "$LOCAL_URL"
fi
echo "Keep this terminal open while you use the interface. Ctrl+C closes it; runs keep going on $HOST."
wait "$REMOTE_PID"
