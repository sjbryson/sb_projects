#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026
"""

import logging
from pathlib import Path
from typing import Optional
from contextlib import contextmanager

def setup_logger (
        logger: Optional[logging.Logger] = None,
        level: int = logging.INFO,
        log_file: Optional[str | Path] = None,
        debug: bool = False,
        name: Optional[str] = '_default_'
    ) -> logging.Logger:

    if logger is None or not isinstance(logger, logging.Logger):
        logger = logging.getLogger(name)

    logger.setLevel(logging.DEBUG if debug else level)
    
    if not logger.handlers:
        if log_file is not None:
            log_file = Path(log_file)
        else:
            log_file = Path.cwd() / f"{name}_log_file.txt"

    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    file_handler = logging.FileHandler(log_file, delay=True)
    file_handler.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # If debug is True, set a console handler with DEBUG level
    if debug:
        stream_handler = logging.StreamHandler()
        stream_formatter = logging.Formatter('%(levelname)s - %(asctime)s - %(message)s')
        stream_handler.setFormatter(stream_formatter)
        stream_handler.setLevel(logging.DEBUG)
        logger.addHandler(stream_handler)

    return logger

def shutdown_logger(logger: logging.Logger):
    """Closes all handlers and removes them from the logger."""
    for handler in logger.handlers[:]:  # Iterate over a copy
        handler.close()
        logger.removeHandler(handler)

def switch_logger(old_logger: logging.Logger, new_name: str, new_log_file: str | Path) -> logging.Logger:
    """Convenience to kill one logger and start a fresh one."""
    shutdown_logger(old_logger)
    return setup_logger(name=new_name, log_file=new_log_file)

#def get_log_file_path(logger: logging.Logger) -> Path | None:
#    
#    for handler in logger.handlers:
#        if isinstance(handler, logging.FileHandler):
#            return Path(handler.baseFilename)
#    return None

@contextmanager
def log_session(name: str, log_file: str | Path, level: int = logging.INFO, debug: bool = False):
    """
    Context manager to safely create, use, and shutdown a logger.
    """
    logger = setup_logger(name=name, log_file=log_file, level=level, debug=debug)
    try:
        yield logger
    finally:
        shutdown_logger(logger)


'''
Example usage:

from src.config_manager import ConfigManager
from src.wrappers import SomeBioTool
from src.logging_utilities import log_session
from src.subprocess_utilities import run_cmd

def run_pipeline():
    cfg = ConfigManager("manifest.tsv")

    for row in cfg:
        sample_id = row['sample_id']
        log_path = f"logs/{sample_id}_run.log"

        # The 'with' block handles the start and stop of the logger
        with log_session(name=sample_id, log_file=log_path) as log:
            log.info(f"--- Starting Analysis for {sample_id} ---")
            
            try:
                # Build the command
                tool = SomeBioTool(**row)
                cmd = tool.build()
                log.debug(f"Generated Command: {cmd}")
                
                # subprocess utility is called here
                # run_cmd(cmd, logger=log)
                
                log.info(f"--- Finished Analysis for {sample_id} ---")

            except Exception as e:
                # This will be logged to the sample-specific log
                log.error(f"Failed to process {sample_id}: {str(e)}")
                # Depending on your needs, you can 'continue' to next sample 
                # or 'raise' to stop the whole pipeline
                continue 

# Once the loop exits the 'with' block, the file handler is released.
'''