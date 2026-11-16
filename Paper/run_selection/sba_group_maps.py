"""SBA group maps of the paper datasets ('All' category): group_maps.main() with the paper configuration.

    python Paper/run_selection/sba_group_maps.py [--species Rat Mouse ...] [--min-n 4] [--output <folder>]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from Paper.run_selection import config_ednix  # noqa: E402,F401  (configures group_maps)
from Statistics.Group_fMRI.run_selection import group_maps  # noqa: E402

if __name__ == '__main__':
    group_maps.main()
