# EDNiX interface

3D atlases (several at once, volume and cortical surface), launcher editing without code for the anatomical, functional and PET pipelines, background runs with live logs.

## Installation (once)

1. Copy the `GUI` folder into the EDNiX folder (next to `modalities/`, `Tools/`...).
2. In a terminal (a double-click on a `.sh` file often opens it in a text editor instead of running it):

       bash /home/cgarin/PycharmProjects/EDNiX/GUI/install_shortcut.sh

   An EDNiX icon appears on the desktop and in the applications menu.
3. Usually nothing else is needed: the launcher finds the Python of the EDNiX environment (venv in `~/PycharmProjects`, `<EDNiX>/venv`, conda environments). Otherwise set it in `ednix_gui.conf`: one line `EDNIX_PYTHON=/path/to/bin/python`, **without `#` in front**.

## If nothing starts

Run it once from a terminal to see the messages:

    bash /home/cgarin/PycharmProjects/EDNiX/GUI/launch_EDNiX_GUI.sh

After a double-click the same messages are in `~/.ednix_gui/launch.log`, and an error window opens if the launcher fails.

## Use

- **Atlases**: pick a species, a reference, a source (volume atlas or cortical surface), an atlas, then *Add to view*. Add several layers to compare levels, atlases or species; *Side by side* and *Same size* help comparisons. Click a region to select it, double-click to centre the view, drag to rotate, scroll to zoom. *View* buttons give standard orientations; *Fix orientation* rotates everything by 90° if a template is stored in an unusual orientation. The *Challenge* asks you to find 10 regions. *Use as seeds* sends the selection to the functional launcher.
- **Launchers** (Anatomy, Functional, PET): open an example or your own launcher; every parameter shows its description (from the pipeline docstrings) and the launcher comment; *Skip_step* is a checklist of the pipeline steps. *Check*, *Save*, *Run in the background*.
- **Runs**: status, current pipeline step and live log of every run. Runs keep going if you close the interface; reopen it to follow them.
- **Quality control**, for a dataset:
  - images of `<BIDS>/QC` per check, and the reports of every subject: functional (`func/QC`: QC summary of step 12, connectivity specificity of step 13) and anatomical (`anat/QC_anat`, step 16), with a subject filter and a full-screen viewer (arrow keys);
  - the QC metrics of all runs (and of all subjects for the anatomy) and the connectivity specificity gathered in one table each, sortable, with values more than 3.5 robust standard deviations from the median highlighted; clicking a row shows the images of that run;
  - the errors and warnings of the processing logs of every subject (most errors first), and each log with its errors highlighted, alone or with two lines of context.
- **Runs**: the log keeps the colours of the pipeline messages, highlights errors and warnings, counts them, and can show them only or jump from one error to the next.
- **Linked selections** (Atlases, Layers box): selecting a region also selects the region with the same name (left/right prefixes ignored) or the same label ID in the other layers, to compare a region across species, levels or surfaces.

## Interface on a remote server

To use the interface of a server (for example the GPU server) from your own computer, set once in `GUI/ednix_gui.conf` of your computer:

    EDNIX_REMOTE_HOST=gpu-07-01
    EDNIX_REMOTE_REPO=/home/cgarin/PycharmProjects/EDNiX

then run from the EDNiX folder of your computer:

    bash GUI/open_remote.sh

It starts the interface on the server, opens an SSH tunnel and the local browser. Keep the terminal open while you use it; Ctrl+C closes the interface, runs keep going on the server.

## Notes

- The server only listens on this computer (127.0.0.1) and requires a random access key given by the launcher: other users of a shared machine cannot use it.
- Volume atlas surfaces are computed at the first display and cached in `~/.ednix_gui/cache`. Cortical surfaces are read from `<species>/<reference>/surfaces/Native_resol`.
- Region names come from `<species>/label_code/<atlas>_label.txt` (or `_StatsLUT.txt`, `_ITK.txt`); for cortical surfaces, from the label table of the `.label.gii` files.
- Run logs are in `~/.ednix_gui/`.
- Runs need the pipeline tools (Singularity containers, AFNI, ANTs...): on Windows only the atlas viewer is usable.
- three.js (MIT licence) is bundled in `web/vendor` so the interface works offline.
