
import shutil
from pathlib import Path
from subprocess import check_output

def clean_dir(directory: Path | str, dry_run: bool = False) -> None:
    path = Path(directory)
    if dry_run:
        message = f"[DRY RUN] Cleaned directory {path}."
    elif not path.exists() or not path.is_dir():
        message = f"Directory {path} does not exist, skipping clean."
    else:
        for item in path.iterdir():
            if item.is_file() or item.is_symlink():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
        message = f"Cleaned directory: {path}."
    print(message)


def remove_dir(directory: Path | str, dry_run: bool = False) -> None:
    path = Path(directory)
    if dry_run:
        message = f"[DRY RUN] Removed directory {path}."
    elif not path.exists() or not path.is_dir():
        message = f"Directory {path} does not exist, skipping clean."
    else:
        shutil.rmtree(path)
        message = f"Removed directory {path}."
        print(message)


def delete_file(file_path: Path | str, dry_run: bool = False) -> None:
    path = Path(file_path)
    if dry_run:
        message = f"[DRY RUN] Deleted file {path}."
    elif not path.exists():
        raise FileNotFoundError(f"Missing file: {path}, skipping delete.")
    else:
        path.unlink()
        message = f"Deleted file {path}."
    print(message)


def count_lines(file_path: Path | str, dry_run: bool = False) -> int:
    path = Path(file_path)
    if dry_run:
        return 0
    elif not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    else:
        result = check_output(['wc', '-l', str(path)], text=True)
        return int(result.split()[0])