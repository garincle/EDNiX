"""Where EDNiX and the data are, without any path written in the launchers.

MAIN_PATH  the EDNiX folder (this repository, which also holds Atlases_library/ and Singularity_library/),
           found from the location of this file; ends with a path separator, as the pipeline expects.
DATA_ROOT  the folder that holds the BIDS datasets of the example launchers (<DATA_ROOT>/<Species>/<BIDS folder>):
           environment variable EDNIX_DATA, default <MAIN_PATH>/data.

In a launcher, ``opj(DATA_ROOT, 'Rat/BIDS_Grandjean/')`` gives the dataset folder; an absolute path in place of
'Rat/BIDS_Grandjean/' is used as is (os.path.join keeps the last absolute component), so a launcher can also point
anywhere without setting EDNIX_DATA.
"""
import os

MAIN_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + os.sep
DATA_ROOT = os.environ.get('EDNIX_DATA') or os.path.join(MAIN_PATH, 'data')


def resolve_data_path(path):
    """Dataset folder of a launcher written relative to DATA_ROOT (absolute paths are returned unchanged)."""
    if not path:
        return path
    return os.path.join(DATA_ROOT, os.path.expanduser(str(path)))
