
"""
Samuel Joseph Bryson
Copyright 2026
"""

from dataclasses import dataclass, asdict, field, fields, InitVar
from pathlib import Path
from typing import Optional

@dataclass(kw_only=True)
class Wrapper():
    dry_run: InitVar[bool] = False
    cmd:     str = field(default="", init=False, repr=False)

    def __post_init__(self, dry_run: bool):
        for f in fields(self):
            metadata_type = f.metadata.get('type', '')
            value = getattr(self, f.name)

            if value is None:
                continue

            # --- Value Flags (e.g., threads, quality scores) --- #
            if metadata_type == 'value_flag':
                flag_fmt = f.metadata.get('flag_fmt') # (e.g., "-t {value}"| allows piped cmds with dif. thread flags, e.g. t, w, @)
                if value is not None and value is not False:
                    setattr(self, f.name, flag_fmt.format(value=value))
                else:
                    setattr(self, f.name, "")

            # --- Bool Flags --- #
            elif metadata_type == 'flag':
                flag_string = f.metadata.get('option')
                # If the attribute is True, use the flag; otherwise, empty string
                setattr(self, f.name, flag_string if value else "")

            # --- IO Paths --- #
            elif metadata_type in ['input_file', 'output_file'] and value:
                path_obj = Path(value)
                setattr(self, f.name, path_obj)
                #if metadata_type == 'input_file' and not path_obj.exists():
                #    raise FileNotFoundError(f"Missing: {path_obj}")
                if metadata_type == 'input_file' and not dry_run:
                    if not path_obj.exists():
                        raise FileNotFoundError(f"Missing: {path_obj}")
                if metadata_type == 'output_file':
                    path_obj.parent.mkdir(parents=True, exist_ok=True)
    
    def build(self) -> str:
        """Format self.cmd string with attribute values."""
        return " ".join(self.cmd.format(**asdict(self)).split())
