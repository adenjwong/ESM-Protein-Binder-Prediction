from pathlib import Path
from Bio import SeqIO


def load_fasta_sequence(fasta_path: str | Path) -> tuple[str, str]:
    """
    Load a single protein sequence from a FASTA file.

    Returns:
        header: FASTA record ID/name
        sequence: amino acid sequence as a string
    """
    fasta_path = Path(fasta_path)

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    records = list(SeqIO.parse(str(fasta_path), "fasta"))

    if len(records) == 0:
        raise ValueError(f"No FASTA records found in: {fasta_path}")

    if len(records) > 1:
        raise ValueError(
            f"Expected one FASTA record in {fasta_path}, but found {len(records)}"
        )

    record = records[0]
    header = record.id
    sequence = str(record.seq).upper()

    return header, sequence