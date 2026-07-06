"""Console entry points for CAFE.

Installing this package exposes the ``cafe`` command (and ``cafe-web`` for the
hosted layout). Both boot the Streamlit app the same way the old ``cafe.py`` /
``web.py`` wrappers did, but resolve the entry script relative to this file so
the command works from any working directory.

The theme and telemetry settings that used to live in ``.streamlit/config.toml``
are passed here as config overrides, since that file is only picked up when
Streamlit is launched from the project directory.
"""
import os

from streamlit.web import cli

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))

# Mirrors the former .streamlit/config.toml so the installed command looks
# identical no matter where it is launched from.
_CONFIG_FLAGS = [
    "--theme.primaryColor=#0f5070",
    "--theme.backgroundColor=#ffffff",
    "--theme.secondaryBackgroundColor=#f8f9fa",
    "--theme.textColor=#1a262a",
    "--theme.font=sans serif",
    "--browser.gatherUsageStats=false",
]


def _launch(script_name, max_upload_mb):
    script = os.path.join(_PKG_DIR, script_name)
    cli.main_run(
        [script, "--server.maxUploadSize", str(max_upload_mb), *_CONFIG_FLAGS]
    )


def main():
    """Desktop app — the ``cafe`` command."""
    _launch("run.py", 3000)


def main_web():
    """Hosted deployment layout — the ``cafe-web`` command."""
    _launch("run_w.py", 2000)


if __name__ == "__main__":
    main()
