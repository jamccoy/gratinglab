"""
Command-line entry point.

Mode selection still comes from config.py, which is the interface this project
has always had and which the README documents. This module exists so there is one
dispatch shared by `python main.py` and the `afm-analysis` console script,
instead of the routing living in a root-level script that only works from a
checkout.
"""
import sys

import matplotlib.pyplot as plt


def main(argv=None):
    """Run the analysis mode selected in config.py"""
    from .config import ANALYSIS_MODE
    from .workflows import (
        run_single_file_analysis,
        run_multiple_file_analysis,
        run_comparison_analysis,
        run_boundary_profile_export,
        run_icc_report,
    )

    modes = {
        'single': run_single_file_analysis,
        'multiple': run_multiple_file_analysis,
        'compare': run_comparison_analysis,
        'ggp': run_boundary_profile_export,
        'icc': run_icc_report,
    }

    plt.close("all")
    plt.ion()
    plt.rcParams.update({"text.usetex": False, "font.family": "sans-serif"})

    run = modes.get(ANALYSIS_MODE)
    if run is None:
        print(f"Unknown ANALYSIS_MODE: {ANALYSIS_MODE}")
        print(f"Set ANALYSIS_MODE in config.py to one of: "
              f"{', '.join(sorted(modes))}")
        return 1

    run()
    plt.show()
    print("\n✓ Analysis complete!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
