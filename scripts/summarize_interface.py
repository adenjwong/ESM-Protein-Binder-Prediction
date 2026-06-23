from argparse import ArgumentParser
from pathlib import Path

import pandas as pd


def summarize_interface_contacts(
    contacts_csv: str | Path,
    output_csv: str | Path,
    top_n: int = 20,
) -> None:
    contacts_csv = Path(contacts_csv)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if not contacts_csv.exists():
        raise FileNotFoundError(f"Contact CSV not found: {contacts_csv}")

    contacts = pd.read_csv(contacts_csv)

    if contacts.empty:
        print("Contact file is empty. No residue-level contacts to summarize.")
        contacts.to_csv(output_csv, index=False)
        return

    required_columns = [
        "chain_a",
        "residue_a_id",
        "residue_a_name",
        "chain_b",
        "residue_b_id",
        "residue_b_name",
        "distance_angstrom",
    ]

    missing_columns = [
        column for column in required_columns if column not in contacts.columns
    ]

    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    residue_summary = (
        contacts.groupby(
            [
                "chain_a",
                "residue_a_id",
                "residue_a_name",
                "chain_b",
                "residue_b_id",
                "residue_b_name",
            ],
            as_index=False,
        )
        .agg(
            closest_distance_angstrom=("distance_angstrom", "min"),
            mean_distance_angstrom=("distance_angstrom", "mean"),
            atom_contact_count=("distance_angstrom", "count"),
        )
        .sort_values(
            ["closest_distance_angstrom", "atom_contact_count"],
            ascending=[True, False],
        )
    )

    residue_summary["residue_pair"] = (
        residue_summary["chain_a"].astype(str)
        + ":"
        + residue_summary["residue_a_name"].astype(str)
        + residue_summary["residue_a_id"].astype(str)
        + " - "
        + residue_summary["chain_b"].astype(str)
        + ":"
        + residue_summary["residue_b_name"].astype(str)
        + residue_summary["residue_b_id"].astype(str)
    )

    columns = [
        "residue_pair",
        "chain_a",
        "residue_a_id",
        "residue_a_name",
        "chain_b",
        "residue_b_id",
        "residue_b_name",
        "closest_distance_angstrom",
        "mean_distance_angstrom",
        "atom_contact_count",
    ]

    residue_summary = residue_summary[columns]
    residue_summary.to_csv(output_csv, index=False)

    chain_a_interface = (
        residue_summary[["chain_a", "residue_a_id", "residue_a_name"]]
        .drop_duplicates()
        .sort_values(["chain_a", "residue_a_id"])
    )

    chain_b_interface = (
        residue_summary[["chain_b", "residue_b_id", "residue_b_name"]]
        .drop_duplicates()
        .sort_values(["chain_b", "residue_b_id"])
    )

    print("Residue-level interface summary completed.")
    print("Input atom-contact CSV:", contacts_csv)
    print("Output residue summary CSV:", output_csv)
    print()
    print("Residue-residue contact pairs:", len(residue_summary))
    print("Chain A interface residues:", len(chain_a_interface))
    print("Chain B interface residues:", len(chain_b_interface))
    print()

    print(f"Top {min(top_n, len(residue_summary))} closest residue-residue contacts:")
    print(
        residue_summary[
            [
                "residue_pair",
                "closest_distance_angstrom",
                "mean_distance_angstrom",
                "atom_contact_count",
            ]
        ]
        .head(top_n)
        .to_string(index=False)
    )


def parse_args():
    parser = ArgumentParser(
        description="Summarize atom-level interface contacts into residue-level contact pairs."
    )

    parser.add_argument(
        "--contacts",
        required=True,
        help="Path to atom-level interface contact CSV.",
    )

    parser.add_argument(
        "--out",
        default="outputs/interface_contacts/residue_interface_summary.csv",
        help="Output CSV path for residue-level interface summary.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of closest residue-residue contacts to print.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    summarize_interface_contacts(
        contacts_csv=args.contacts,
        output_csv=args.out,
        top_n=args.top_n,
    )


if __name__ == "__main__":
    main()