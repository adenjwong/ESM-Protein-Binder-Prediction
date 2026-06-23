from argparse import ArgumentParser
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_contact_map(
    residue_summary_csv: str | Path,
    output_png: str | Path,
    value_column: str = "atom_contact_count",
) -> None:
    residue_summary_csv = Path(residue_summary_csv)
    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    if not residue_summary_csv.exists():
        raise FileNotFoundError(f"Residue summary CSV not found: {residue_summary_csv}")

    summary = pd.read_csv(residue_summary_csv)

    if summary.empty:
        raise ValueError("Residue summary CSV is empty. No contact map to plot.")

    required_columns = [
        "residue_a_id",
        "residue_a_name",
        "residue_b_id",
        "residue_b_name",
        value_column,
    ]

    missing_columns = [
        column for column in required_columns if column not in summary.columns
    ]

    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    contact_matrix = summary.pivot_table(
        index="residue_a_id",
        columns="residue_b_id",
        values=value_column,
        aggfunc="sum",
        fill_value=0,
    )

    chain_a_labels = {
        row.residue_a_id: f"{row.residue_a_name}{row.residue_a_id}"
        for row in summary.itertuples(index=False)
    }

    chain_b_labels = {
        row.residue_b_id: f"{row.residue_b_name}{row.residue_b_id}"
        for row in summary.itertuples(index=False)
    }

    y_labels = [
        chain_a_labels[residue_id]
        for residue_id in contact_matrix.index
    ]

    x_labels = [
        chain_b_labels[residue_id]
        for residue_id in contact_matrix.columns
    ]

    plt.figure(figsize=(10, 8))
    plt.imshow(contact_matrix.values, aspect="auto", origin="upper")

    plt.xticks(
        ticks=range(len(x_labels)),
        labels=x_labels,
        rotation=90,
    )

    plt.yticks(
        ticks=range(len(y_labels)),
        labels=y_labels,
    )

    plt.xlabel("Chain B residues")
    plt.ylabel("Chain A residues")
    plt.title(f"Predicted interface contact map ({value_column})")

    colorbar = plt.colorbar()
    colorbar.set_label(value_column)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()

    print("Contact map saved to:")
    print(output_png)
    print()
    print("Matrix shape:", contact_matrix.shape)
    print("Chain A residues shown:", len(y_labels))
    print("Chain B residues shown:", len(x_labels))


def parse_args():
    parser = ArgumentParser(
        description="Plot a residue-level contact map from an interface summary CSV."
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Path to residue-level interface summary CSV.",
    )

    parser.add_argument(
        "--out",
        default="outputs/interface_contacts/contact_map.png",
        help="Output PNG path for the contact map.",
    )

    parser.add_argument(
        "--value-column",
        default="atom_contact_count",
        choices=[
            "atom_contact_count",
            "closest_distance_angstrom",
            "mean_distance_angstrom",
        ],
        help="Column to visualize in the contact map.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    plot_contact_map(
        residue_summary_csv=args.summary,
        output_png=args.out,
        value_column=args.value_column,
    )


if __name__ == "__main__":
    main()