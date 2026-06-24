from argparse import ArgumentParser
from pathlib import Path
import json

import pandas as pd


def make_interface_quality_report(
    contacts_csv: str | Path,
    residue_summary_csv: str | Path,
    output_json: str | Path,
    clash_cutoff: float = 2.0,
    close_contact_cutoff: float = 3.5,
    top_n: int = 10,
) -> None:
    contacts_csv = Path(contacts_csv)
    residue_summary_csv = Path(residue_summary_csv)
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    if not contacts_csv.exists():
        raise FileNotFoundError(f"Contact CSV not found: {contacts_csv}")

    if not residue_summary_csv.exists():
        raise FileNotFoundError(f"Residue summary CSV not found: {residue_summary_csv}")

    contacts = pd.read_csv(contacts_csv)
    summary = pd.read_csv(residue_summary_csv)

    if contacts.empty or summary.empty:
        report = {
            "interface_detected": False,
            "notes": ["No contacts were found between the selected chains."],
        }

        with open(output_json, "w") as f:
            json.dump(report, f, indent=2)

        print("No interface detected.")
        print("Saved report to:", output_json)
        return

    min_distance = float(contacts["distance_angstrom"].min())

    clash_contacts = contacts[contacts["distance_angstrom"] < clash_cutoff]
    close_contacts = contacts[contacts["distance_angstrom"] < close_contact_cutoff]

    chain_a_interface_residues = (
        summary[["chain_a", "residue_a_id", "residue_a_name"]]
        .drop_duplicates()
        .sort_values(["chain_a", "residue_a_id"])
    )

    chain_b_interface_residues = (
        summary[["chain_b", "residue_b_id", "residue_b_name"]]
        .drop_duplicates()
        .sort_values(["chain_b", "residue_b_id"])
    )

    top_by_contact_count = (
        summary.sort_values(
            ["atom_contact_count", "closest_distance_angstrom"],
            ascending=[False, True],
        )
        .head(top_n)
    )

    top_by_closest_distance = (
        summary.sort_values(
            ["closest_distance_angstrom", "atom_contact_count"],
            ascending=[True, False],
        )
        .head(top_n)
    )

    notes = []

    if len(clash_contacts) > 0:
        notes.append(
            f"Possible clashes detected: {len(clash_contacts)} atom contacts are below {clash_cutoff:.1f} Å."
        )

    if min_distance < clash_cutoff:
        notes.append(
            f"Minimum atom distance is very short: {min_distance:.3f} Å."
        )

    if len(summary) < 5:
        notes.append("Interface appears sparse based on residue-residue contact count.")

    if len(summary) >= 5:
        notes.append("Interface contacts were detected between the two chains.")

    report = {
        "interface_detected": True,
        "minimum_atom_distance_angstrom": min_distance,
        "atom_contacts_total": int(len(contacts)),
        "residue_residue_contact_pairs": int(len(summary)),
        "chain_a_interface_residues": int(len(chain_a_interface_residues)),
        "chain_b_interface_residues": int(len(chain_b_interface_residues)),
        "possible_clash_contacts_below_cutoff": int(len(clash_contacts)),
        "close_contacts_below_cutoff": int(len(close_contacts)),
        "clash_cutoff_angstrom": clash_cutoff,
        "close_contact_cutoff_angstrom": close_contact_cutoff,
        "top_pairs_by_contact_count": top_by_contact_count[
            [
                "residue_pair",
                "closest_distance_angstrom",
                "mean_distance_angstrom",
                "atom_contact_count",
            ]
        ].to_dict(orient="records"),
        "top_pairs_by_closest_distance": top_by_closest_distance[
            [
                "residue_pair",
                "closest_distance_angstrom",
                "mean_distance_angstrom",
                "atom_contact_count",
            ]
        ].to_dict(orient="records"),
        "notes": notes,
    }

    with open(output_json, "w") as f:
        json.dump(report, f, indent=2)

    print("Interface quality report")
    print("------------------------")
    print(f"Minimum atom distance: {min_distance:.3f} Å")
    print("Total atom contacts:", len(contacts))
    print("Residue-residue contact pairs:", len(summary))
    print("Chain A interface residues:", len(chain_a_interface_residues))
    print("Chain B interface residues:", len(chain_b_interface_residues))
    print(f"Possible clash contacts < {clash_cutoff:.1f} Å:", len(clash_contacts))
    print(f"Close contacts < {close_contact_cutoff:.1f} Å:", len(close_contacts))
    print()

    print("Notes:")
    for note in notes:
        print("-", note)

    print()
    print("Top residue pairs by contact count:")
    print(
        top_by_contact_count[
            [
                "residue_pair",
                "closest_distance_angstrom",
                "atom_contact_count",
            ]
        ].to_string(index=False)
    )

    print()
    print("Saved report to:")
    print(output_json)


def parse_args():
    parser = ArgumentParser(
        description="Create a quick quality report for a predicted protein-protein interface."
    )

    parser.add_argument(
        "--contacts",
        required=True,
        help="Path to atom-level contact CSV.",
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Path to residue-level interface summary CSV.",
    )

    parser.add_argument(
        "--out",
        default="outputs/interface_contacts/interface_quality_report.json",
        help="Output JSON report path.",
    )

    parser.add_argument(
        "--clash-cutoff",
        type=float,
        default=2.0,
        help="Atom distance below this value is flagged as a possible clash.",
    )

    parser.add_argument(
        "--close-contact-cutoff",
        type=float,
        default=3.5,
        help="Atom distance below this value is counted as a close contact.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of top residue pairs to include in the report.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    make_interface_quality_report(
        contacts_csv=args.contacts,
        residue_summary_csv=args.summary,
        output_json=args.out,
        clash_cutoff=args.clash_cutoff,
        close_contact_cutoff=args.close_contact_cutoff,
        top_n=args.top_n,
    )


if __name__ == "__main__":
    main()