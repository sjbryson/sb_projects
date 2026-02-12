#!/usr/bin/env python3

"""
Samuel Joseph Bryson
Copyright 2026
"""

from pathlib import Path
import gzip
import subprocess

def fasta_to_dict (fasta: str) -> dict:
    """Open a fasta file (possibly compressed) and add each sequence to a dictionary.
    - removes the ">" from the sequence identifier
    - removes sequence identifier text after first white space
    - removes "*" from end of sequence if present
    Args:
        input_fasta: (str) path to fasta file
    Returns:
        fasta_dict: dict[str:str] sequence_id to sequence dictionary
    """
    try:
        fasta = Path(fasta).resolve()
        if fasta.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
        fasta_dict = dict()
        with open_func(fasta, mode) as f:
            seq_id = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    seq_id = line[1:].split()[0]
                    fasta_dict[seq_id] = ''
                else:
                    if line.endswith('*'):
                        line = line[:-1]
                    fasta_dict[seq_id] += line
    except FileNotFoundError as e:
        #message = f"Fasta file {fasta} not found - {e}."
        raise
    except gzip.BadGzipFile as e:
        #message = f"Corrupted or bad gzip file {fasta} - {e}."
        raise
    return fasta_dict


def write_fasta_from_dict (fasta_dict: dict, fasta: str, gzip_output: bool = False) -> None:
        """Write the sequences in a fasta dictionary to a file, with optional gzip compression.
        Args:
            fasta_dict: (dict[str:str]) keys=sequence identifiers, values=sequences
            fasta: (str) path for fasta file to write
            gzip_output: (bool) whether to compress the fasta file with gzip
        """
        try:
            fasta = Path(fasta).resolve(strict=False)
       
            if gzip_output and fasta.suffix != '.gz':
                output_fasta = fasta.with_suffix(fasta.suffix + '.gz')
           
            if gzip_output:
                open_func = gzip.open
                mode = 'wt'
            else:
                open_func = open
                mode = 'w'

            with open_func(fasta, mode) as f:
                for k, v in fasta_dict.items():
                    f.write(f">{k}\n{v}\n")
        except Exception as e:
            #message = f"write_fasta_from_dict() failed - {e}."
            raise

def fastq_count_reads(file_path: Path | str, dry_run: bool = False) -> int:
    """Count reads in a FASTQ file by dividing line count by 4."""
    path = Path(file_path)
    
    if dry_run:
        return 0
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    is_gz = '.gz' in path.suffixes
    cmd_base = "gzip -dc" if is_gz else "cat"
    
    try:
        cmd = f"{cmd_base} {path} | wc -l"
        output = subprocess.check_output(cmd, shell=True, text=True)
        line_count = int(output.strip())
        return line_count // 4
    
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to count reads for {path}. Error: {e}")