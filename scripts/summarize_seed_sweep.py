from argparse import ArgumentParser
from pathlib import Path
import json
from collections import Counter

import pandas as pd


def residue_label(chain: str, residue_name: str, residue_id: int) -> str:
    return f"{chain}:{residue_name}{residue_id}"


def summarize_seed_sweep(interface_dir: str | Path, output_csv: str | Path) -> None:
    interface_dir = Path(interface_dir)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    quality_reports = sorted(interface_dir.glob("*_quality_report.json"))

    if not quality_reports:
        raise FileNotFoundError(f"No quality report JSON files found in {interface_dir}")

    rows = []
    chain_a_residue_counter = Counter()
    chain_b_residue_counter = Counter()
    residue_pair_counter = Counter()

    for report_path in quality_reports:
        name = report_path.name.replace("_quality_report.json", "")
        seed = name.replace("complex_seed_", "")

        summary_path = interface_dir / f"{name}_residue_summary.csv"

        with open(report_path, "r") as f:
            report = json.load(f)

        row = {
            "name": name,
            "seed": seed,
            "interface_detected": report.get("interface_detected", False),
            "minimum_atom_distance_angstrom": report.get("minimum_atom_distance_angstrom"),
            "atom_contacts_total": report.get("atom_contacts_total"),
            "residue_residue_contact_pairs": report.get("residue_residue_contact_pairs"),
            "chain_a_interface_residues": report.get("chain_a_interface_residues"),
            "chain_b_interface_residues": report.get("chain_b_interface_residues"),
            "possible_clash_contacts_below_cutoff": report.get("possible_clash_contacts_below_cutoff"),
            "close_contacts_below_cutoff": report.get("close_contacts_below_cutoff"),
        }

        rows.append(row)

        if summary_path.exists():
            summary = pd.read_csv(summary_path)

            for residue in (
                summary[["chain_a", "residue_a_name", "residue_a_id"]]
                .drop_duplicates()
                .itertuples(index=False)
            ):
                chain_a_residue_counter[
                    residue_label(residue.chain_a, residue.residue_a_name, residue.residue_a_id)
                ] += 1

            for residue in (
                summary[["chain_b", "residue_b_name", "residue_b_id"]]
                .drop_duplicates()
                .itertuples(index=False)
            ):
                chain_b_residue_counter[
                    residue_label(residue.chain_b, residue.residue_b_name, residue.residue_b_id)
                ] += 1

            for pair in summary["residue_pair"].drop_duplicates():
                residue_pair_counter[pair] += 1

    results = pd.DataFrame(rows).sort_values("seed")
    results.to_csv(output_csv, index=False)

    print("Seed sweep summary saved to:")
    print(output_csv)
    print()

    print("Per-seed summary:")
    print(results.to_string(index=False))
    print()

    print("Chain A residues recurring across seeds:")
    for residue, count in chain_a_residue_counter.most_common():
        if count >= 2:
            print(f"  {residue}: {count} seeds")

    print()
    print("Chain B residues recurring across seeds:")
    for residue, count in chain_b_residue_counter.most_common():
        if count >= 2:
            print(f"  {residue}: {count} seeds")

    print()
    print("Residue-residue pairs recurring across seeds:")
    recurring_pairs = [
        (pair, count)
        for pair, count in residue_pair_counter.most_common()
        if count >= 2
    ]

    if not recurring_pairs:
        print("  No exact residue-residue pairs appeared in 2+ seeds.")
    else:
        for pair, count in recurring_pairs:
            print(f"  {pair}: {count} seeds")


def parse_args():
    parser = ArgumentParser(
        description="Summarize ESMFold2 complex seed-sweep interface results."
    )

    parser.add_argument(
        "--interface-dir",
        default="outputs/seed_sweep/interface_contacts",
        help="Directory containing seed-sweep interface outputs.",
    )

    parser.add_argument(
        "--out",
        default="outputs/seed_sweep/seed_sweep_summary.csv",
        help="Output CSV path.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    summarize_seed_sweep(
        interface_dir=args.interface_dir,
        output_csv=args.out,
    )


if __name__ == "__main__":
    main()