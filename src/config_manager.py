
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
    """A lightweight wrapper for using a Pandas df for workflow data."""
    config_file: Path | str
    delimiter:   str = '\t'
    config_df:   pd.DataFrame = field(init=False, repr=False)

    def __post_init__(self):
        
        # 1. Convert config_file to Path and check existence
        if isinstance(self.config_file, str):
            self.config_file = Path(self.config_file)
        
        if not self.config_file.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_file}")

        # 2. Load the config_file into a dataframe
        self.config_df = pd.read_csv(
            self.config_file,
            sep = self.delimiter,
            header=0,
            index_col=False
            )
        
    def __iter__(self) -> Iterator[Dict[str, Any]]:
        """Yields each row of the dataframe as a dictionary."""
        # to_dict('records') returns a list of dicts, one for each row
        for index, row in self.config_df.iterrows():
            yield index, row.to_dict()

    def add_column(self, name: str, default_value: Any = None):
        """Creates a new field/column in the config."""
        self.config_df[name] = default_value

    def update_row(self, index: Any, column: str, value: Any):
        """Updates a specific row's value for a field. """
        self.config_df.at[index, column] = value

    def save(self, path: Optional[Path | str] = None):
        """Saves the current state back to disk."""
        out_path = path if path else self.config_file
        self.config_df.to_csv(out_path, sep=self.delimiter, index=False, na_rep="")