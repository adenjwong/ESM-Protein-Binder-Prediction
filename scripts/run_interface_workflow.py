from argparse import ArgumentParser
from pathlib import Path

from scripts.analyze_interface import analyze_interface
from scripts.summarize_interface import summarize_interface_contacts
from scripts.plot_contact_map import plot_contact_map
from scripts.interface_quality_report import make_interface_quality_report


def run_interface_workflow(
    cif_path: str | Path,
    name: str | None = None,
    output_dir: str | Path = "outputs/interface_contacts",
    chain_a: str = "A",
    chain_b: str = "B",
    cutoff: float = 5.0,
) -> None:
    cif_path = Path(cif_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if name is None:
        name = cif_path.stem

    contacts_csv = output_dir / f"{name}_contacts.csv"
    residue_summary_csv = output_dir / f"{name}_residue_summary.csv"
    contact_map_png = output_dir / f"{name}_contact_map.png"
    quality_report_json = output_dir / f"{name}_quality_report.json"

    print("Running full interface workflow...")
    print("CIF:", cif_path)
    print("Name:", name)
    print("Output directory:", output_dir)
    print()

    print("=" * 80)
    print("Step 1: Atom-level interface contact analysis")
    print("=" * 80)

    analyze_interface(
        cif_path=cif_path,
        output_csv=contacts_csv,
        chain_a=chain_a,
        chain_b=chain_b,
        cutoff=cutoff,
        heavy_atoms_only=True,
    )

    print()
    print("=" * 80)
    print("Step 2: Residue-level interface summary")
    print("=" * 80)

    summarize_interface_contacts(
        contacts_csv=contacts_csv,
        output_csv=residue_summary_csv,
        top_n=20,
    )

    print()
    print("=" * 80)
    print("Step 3: Contact map plot")
    print("=" * 80)

    plot_contact_map(
        residue_summary_csv=residue_summary_csv,
        output_png=contact_map_png,
        value_column="atom_contact_count",
    )

    print()
    print("=" * 80)
    print("Step 4: Interface quality report")
    print("=" * 80)

    make_interface_quality_report(
        contacts_csv=contacts_csv,
        residue_summary_csv=residue_summary_csv,
        output_json=quality_report_json,
        clash_cutoff=2.0,
        close_contact_cutoff=3.5,
        top_n=10,
    )

    print()
    print("=" * 80)
    print("Workflow completed")
    print("=" * 80)
    print()
    print("Generated files:")
    print("Atom contacts:", contacts_csv)
    print("Residue summary:", residue_summary_csv)
    print("Contact map:", contact_map_png)
    print("Quality report:", quality_report_json)


def parse_args():
    parser = ArgumentParser(
        description="Run the full interface analysis workflow for a predicted complex CIF."
    )

    parser.add_argument(
        "--cif",
        required=True,
        help="Path to predicted complex .cif file.",
    )

    parser.add_argument(
        "--name",
        default=None,
        help="Base name for output files. Defaults to the CIF filename stem.",
    )

    parser.add_argument(
        "--out-dir",
        default="outputs/interface_contacts",
        help="Directory where interface outputs should be saved.",
    )

    parser.add_argument(
        "--chain-a",
        default="A",
        help="First chain ID.",
    )

    parser.add_argument(
        "--chain-b",
        default="B",
        help="Second chain ID.",
    )

    parser.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="Atom distance cutoff in Angstroms for defining contacts.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    run_interface_workflow(
        cif_path=args.cif,
        name=args.name,
        output_dir=args.out_dir,
        chain_a=args.chain_a,
        chain_b=args.chain_b,
        cutoff=args.cutoff,
    )


if __name__ == "__main__":
    main()