"""Dev entry point for the desktop app (used by ``pixi run cafe`` / conda).

The installable ``cafe`` command lives in ``cafe_app.launcher``; this thin
wrapper simply makes the package importable from a source checkout and
delegates to it, so there is a single launch code path.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from cafe_app.launcher import main

if __name__ == "__main__":
    main()
