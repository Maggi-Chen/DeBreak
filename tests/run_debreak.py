"""Test launcher for the `debreak` script.

`debreak` has no `if __name__ == "__main__"` guard, so under the `spawn`
start method (default on macOS) every worker re-executes the whole pipeline.
This launcher forces `fork` (Linux default) and runs the script unchanged.
"""
import multiprocessing
import os
import runpy
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    sys.path.insert(0, REPO)
    sys.argv = ["debreak"] + sys.argv[1:]
    runpy.run_path(os.path.join(REPO, "debreak"), run_name="__main__")
