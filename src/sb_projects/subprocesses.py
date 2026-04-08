
"""
Samuel Joseph Bryson
Copyright 2026
"""

import subprocess
from pathlib import Path
import logging
from typing import Tuple, Optional

def run_subprocess(formatted_command: str) -> Tuple[int, str, str]:
    """Run and capture everything. Use ONLY for small outputs (e.g. --version)."""
    result = subprocess.run(
        formatted_command, 
        shell=True, 
        capture_output=True, # Shorthand for PIPE on both
        text=True,
        check=False
    )
    return result.returncode, result.stdout, result.stderr

def run_check_call(
        formatted_command: str, 
        devnull = False,
        dry_run = False,
        log_file: Optional[Path | str] = None,
        logger: Optional[logging.Logger] = None
        ) -> None:
    """Standard runner. If log_file is provided, streams output to disk.
    If devnull is set to true, no logs and all outputs to DEVNULL.
    """
    if dry_run:
        msg = f"[DRY RUN]: {formatted_command}"
        if logger:
            logger.info(msg)
        else:
            print(msg)
        return 

    elif devnull:
        subprocess.check_call(
            formatted_command, 
            shell=True, 
            stdout=subprocess.DEVNULL, 
            stderr=subprocess.STDOUT
        )
    elif log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a") as f:
            # use STDOUT to capture errors in the same file
            subprocess.check_call(
                formatted_command, 
                shell=True, 
                stdout=f, 
                stderr=subprocess.STDOUT
            )
    else:
        subprocess.check_call(
            formatted_command, 
            shell=True
        )


def run_check_call_devnull (formatted_command: str) -> None:
    """Run check_call, route stderr & stdout to subprocess.DEVNULL."""
    subprocess.check_call(
        formatted_command, 
        shell=True, 
        stdout=subprocess.DEVNULL, 
        stderr=subprocess.STDOUT
    )

def run_check_output_to_str(
        formatted_command: str,
        dry_run = False,
    ) -> str:    
    """Returns stdout as string. Good for getting a single string, path, or short json string."""
    if dry_run:
        msg = f"[DRY RUN]: {formatted_command}"
        print(msg)
        return "DRY RUN"
    else:
        # check_output returns bytes; using text=True handles decoding automatically
        return subprocess.check_output(formatted_command, shell=True, text=True).strip()


def run_and_log(formatted_command: str, logger: logging.Logger):
    """Runs a command and streams its stdout/stderr line-by-line into 
    the logger.
    """
    logger.info(f"Running command: {formatted_command}")
    
    process = subprocess.Popen(
        formatted_command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, # Merge stderr into stdout
        text=True,
        bufsize=1 # Line buffered
    )

    # Stream output to the logger as it happens
    if process.stdout:
        for line in process.stdout:
            logger.info(line.strip())

    return_code = process.wait()
    
    if return_code != 0:
        logger.error(f"Command failed with return code {return_code}")
        raise subprocess.CalledProcessError(return_code, formatted_command)