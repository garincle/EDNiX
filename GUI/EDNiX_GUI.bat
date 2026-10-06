@echo off
REM Windows: double-click to start the EDNiX interface (runs need the Linux pipeline tools; the atlas viewer works).
cd /d "%~dp0"
if defined EDNIX_PYTHON (set PY=%EDNIX_PYTHON%) else (set PY=python)
"%PY%" ednix_gui_server.py
pause
