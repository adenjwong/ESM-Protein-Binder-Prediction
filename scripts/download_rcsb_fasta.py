from argparse import ArgumentParser
from pathlib import Path
from urllib.request import urlopen


def download_rcsb_fasta(pdb_id: str, output_path: str | Path) -> None:
    pdb_id = pdb_id.upper()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

    print("Downloading RCSB FASTA")
    print("PDB ID:", pdb_id)
    print("URL:", url)
    print("Output:", output_path)
    print()

    with urlopen(url) as response:
        fasta_text = response.read().decode("utf-8")

    if not fasta_text.strip().startswith(">"):
        raise ValueError("Downloaded file does not look like FASTA.")

    output_path.write_text(fasta_text)

    print("Saved FASTA to:")
    print(output_path)
    print()
    print(fasta_text)


def parse_args():
    parser = ArgumentParser(description="Download FASTA records for an RCSB PDB entry.")

    parser.add_argument(
        "--pdb-id",
        required=True,
        help="RCSB PDB ID, for example 2BBN.",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output FASTA path.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    download_rcsb_fasta(
        pdb_id=args.pdb_id,
        output_path=args.out,
    )


if __name__ == "__main__":
    main()