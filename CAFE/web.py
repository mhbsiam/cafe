"""Dev entry point for the hosted (web) layout — used by ``pixi run web``.

Thin wrapper around ``cafe_app.launcher`` so there is a single launch code path.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "src"))

from cafe_app.launcher import main_web

if __name__ == "__main__":
    main_web()
