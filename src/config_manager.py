#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026
"""

import pandas as pd
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict, Any, Iterator

@dataclass
class ConfigManager:
    config_file: Path | str
    delimiter: str = '\t'
    config_fieldnames: Optional[List[str]] = None
    index_col: Optional[str] = None
    df: pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self):
        
        # 1. Convert config_file to Path and check existence
        if isinstance(self.config_file, str):
            self.config_file = Path(self.config_file)
        
        if not self.config_file.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_file}")

        # 2. Load the config_file into a dataframe
        self.df = pd.read_csv(
            self.config_file, 
            sep=self.delimiter, 
            names=self.config_fieldnames,
            index_col=self.index_col,
            header=0 if self.config_fieldnames is None else None
        )

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        """Yields each row of the dataframe as a dictionary."""
        # to_dict('records') returns a list of dicts, one for each row
        for record in self.df.to_dict(orient='records'):
            yield record

    def add_column(self, name: str, default_value: Any = None):
        """Creates a new field/column in the config."""
        self.df[name] = default_value

    def update_row(self, index: Any, column: str, value: Any):
        """Updates a specific row's value for a field. 
        Requires index_col to have been set during init."""
        self.df.at[index, column] = value

    def save(self, path: Optional[Path | str] = None):
        """Saves the current state back to disk."""
        out_path = path if path else self.config_file
        self.df.to_csv(out_path, sep=self.delimiter, index=(self.index_col is not None))