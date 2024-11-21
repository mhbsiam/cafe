import os
from streamlit.web import cli

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cafe_path = os.path.join(current_dir, "bin", "run.py")
    cli.main_run([cafe_path, "--server.maxUploadSize", "3000"])
