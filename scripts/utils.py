
from pathlib import Path

def find_project_root(start_path=None):
    """
    Find the root directory of the project by searching for pyproject.toml file.

    :param start_path: The directory to start the search from. Defaults to the current working directory.
    :return: The path to the project root directory as a pathlib.Path object or None if pyproject.toml is not found.
    """
    if start_path is None:
        start_path = Path.cwd()
    else:
        start_path = Path(start_path)

    current_path = start_path

    while True:
        if (current_path / 'pyproject.toml').is_file():
            return current_path

        parent_path = current_path.parent
        if parent_path == current_path:
            # Reached the root of the filesystem
            return None

        current_path = parent_path
