from argparse import ArgumentParser
from pathlib import Path

from src.fasta import load_fasta_sequence


def test_fasta_loading(fasta_paths: list[str | Path], show_sequence: bool = False) -> None:
    for fasta_path in fasta_paths:
        fasta_path = Path(fasta_path)

        print("=" * 80)
        print("FASTA:", fasta_path)

        header, sequence = load_fasta_sequence(fasta_path)

        print("Header:", header)
        print("Length:", len(sequence))

        if show_sequence:
            print("Sequence:", sequence)

        print()


def parse_args():
    parser = ArgumentParser(
        description="Test loading one or more protein FASTA files."
    )

    parser.add_argument(
        "fastas",
        nargs="+",
        help="One or more FASTA files to load.",
    )

    parser.add_argument(
        "--show-sequence",
        action="store_true",
        help="Print the full protein sequence.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    test_fasta_loading(
        fasta_paths=args.fastas,
        show_sequence=args.show_sequence,
    )


if __name__ == "__main__":
    main()