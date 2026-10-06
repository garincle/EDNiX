#!/usr/bin/env bash
# Starts the EDNiX interface (local server + browser). Double-click the EDNiX icon or run:
#   bash launch_EDNiX_GUI.sh
# Settings (optional) in ednix_gui.conf, same folder: EDNIX_PYTHON=, EDNIX_REPO=, EDNIX_ATLAS_LIBRARY=
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATE="${EDNIX_GUI_STATE:-$HOME/.ednix_gui}"
mkdir -p "$STATE"
LOG="$STATE/launch.log"
# Started by a double-click there is no terminal: everything goes to a log file.
if [ ! -t 1 ]; then exec >"$LOG" 2>&1; fi
echo "[EDNiX] $(date '+%F %T')  launcher: $DIR"

alert() {   # message visible even without a terminal
  echo "$1"
  [ -t 1 ] && return
  if command -v zenity >/dev/null; then zenity --error --title=EDNiX --width=520 --text="$1" 2>/dev/null
  elif command -v kdialog >/dev/null; then kdialog --title EDNiX --error "$1" 2>/dev/null
  elif command -v osascript >/dev/null; then osascript -e "display alert \"EDNiX\" message \"$1\"" 2>/dev/null
  elif command -v notify-send >/dev/null; then notify-send -u critical "EDNiX" "$1"
  elif command -v xmessage >/dev/null; then xmessage -center "EDNiX : $1"
  elif command -v xdg-open >/dev/null; then xdg-open "$LOG"
  fi
}

# ---- settings: only KEY=value lines are read (the file is never executed).
# Comments, a UTF-8 BOM, Windows line endings, quotes and a leading ~ are tolerated.
CONF="$DIR/ednix_gui.conf"
if [ -f "$CONF" ]; then
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line#$'\xef\xbb\xbf'}"; line="${line%$'\r'}"
    line="${line#"${line%%[![:space:]]*}"}"
    case "$line" in
      EDNIX_PYTHON=*|EDNIX_REPO=*|EDNIX_ATLAS_LIBRARY=*)
        key="${line%%=*}"; val="${line#*=}"
        val="${val%\"}"; val="${val#\"}"; val="${val%\'}"; val="${val#\'}"
        [ "${val:0:1}" = "~" ] && val="$HOME${val:1}"
        export "$key=$val"; echo "setting: $key=$val";;
    esac
  done < "$CONF"
fi

is_ednix() { [ -d "$1/modalities" ] && [ -d "$1/Tools" ] && [ -d "$1/Statistics" ]; }
if [ -n "$EDNIX_REPO" ] && ! is_ednix "$EDNIX_REPO"; then
  echo "EDNIX_REPO=$EDNIX_REPO is not the EDNiX folder (it must contain modalities/, Tools/, Statistics/): ignored."
  unset EDNIX_REPO
fi
REPO="${EDNIX_REPO:-}"
[ -z "$REPO" ] && is_ednix "$(dirname "$DIR")" && REPO="$(dirname "$DIR")"
[ -n "$REPO" ] && export EDNIX_REPO="$REPO"
echo "repository: ${REPO:-not found}"

# ---- Python: first one with the full EDNiX stack, else one able to run the interface
has() {   # modules found without importing them (fast)
  [ -n "$1" ] && [ -x "$1" ] && "$1" -c "import importlib.util as u, sys
sys.exit(any(u.find_spec(m.strip()) is None for m in '$2'.split(',')))" >/dev/null 2>&1
}
CANDS=("$EDNIX_PYTHON")
[ -n "$CONDA_PREFIX" ] && CANDS+=("$CONDA_PREFIX/bin/python")
for p in "$REPO"/venv*/bin/python "$REPO"/.venv*/bin/python \
         "$HOME"/PycharmProjects/venv*/bin/python "$HOME"/PycharmProjects/*/venv*/bin/python \
         "$HOME"/PycharmProjects/*/.venv/bin/python "$HOME"/PycharmProjects/*/bin/python \
         "$HOME"/.virtualenvs/*/bin/python \
         "$HOME"/miniconda3/envs/*/bin/python "$HOME"/anaconda3/envs/*/bin/python "$HOME"/miniforge3/envs/*/bin/python \
         "$(command -v python3)" "$(command -v python)"; do
  CANDS+=("$p")
done
PY=""
for mods in "numpy, nibabel, scipy, skimage, pandas, ants, nilearn" "numpy, nibabel, scipy, skimage, pandas"; do
  for c in "${CANDS[@]}"; do
    if has "$c" "$mods"; then PY="$c"; break 2; fi
  done
done
if [ -z "$PY" ]; then
  alert "No Python with numpy, nibabel, scipy, scikit-image and pandas was found.
Set the Python of your EDNiX environment in
$CONF
for example: EDNIX_PYTHON=$HOME/PycharmProjects/EDNiX/venv/bin/python
(Pythons tried: see $LOG)"
  printf 'tried: %s\n' "${CANDS[@]}"
  exit 1
fi
has "$PY" "ants, nilearn" || echo "WARNING: $PY has no ants/nilearn: the interface works but runs will fail."
export EDNIX_PYTHON="$PY"
echo "python: $PY"

"$PY" "$DIR/ednix_gui_server.py" "$@"
rc=$?
if [ $rc -ne 0 ]; then
  alert "The EDNiX interface stopped with an error (code $rc).
Details: $LOG"
fi
exit $rc
