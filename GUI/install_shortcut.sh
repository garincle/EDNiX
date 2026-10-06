#!/usr/bin/env bash
# Creates a double-clickable "EDNiX" icon (Linux: menu + desktop; macOS: EDNiX.command on the Desktop).
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
chmod +x "$DIR/launch_EDNiX_GUI.sh"
if [ "$(uname)" = "Darwin" ]; then
  TARGET="$HOME/Desktop/EDNiX.command"
  printf '#!/bin/bash\n"%s/launch_EDNiX_GUI.sh"\n' "$DIR" > "$TARGET"
  chmod +x "$TARGET"
  echo "Shortcut created: $TARGET"
  exit 0
fi
ENTRY="[Desktop Entry]
Type=Application
Name=EDNiX
Comment=3D atlases and EDNiX pipeline launcher
Exec=\"$DIR/launch_EDNiX_GUI.sh\"
Icon=$DIR/web/icon.svg
Terminal=false
Categories=Science;Education;"
mkdir -p "$HOME/.local/share/applications"
printf '%s\n' "$ENTRY" > "$HOME/.local/share/applications/ednix.desktop"
chmod +x "$HOME/.local/share/applications/ednix.desktop"
DESK="$(xdg-user-dir DESKTOP 2>/dev/null || echo "$HOME/Desktop")"
if [ -d "$DESK" ]; then
  printf '%s\n' "$ENTRY" > "$DESK/EDNiX.desktop"
  chmod +x "$DESK/EDNiX.desktop"
  command -v gio >/dev/null && gio set "$DESK/EDNiX.desktop" metadata::trusted true 2>/dev/null
  echo "Icon created: $DESK/EDNiX.desktop"
fi
echo "Menu entry created: ~/.local/share/applications/ednix.desktop"
