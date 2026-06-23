from argparse import ArgumentParser
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def as_list(value):
    """
    MMCIF2Dict sometimes returns a string if there is one value,
    or a list if there are many. Normalize everything to a list.
    """
    if isinstance(value, list):
        return value
    return [value]


def load_atoms_from_cif(cif_path: str | Path) -> pd.DataFrame:
    """
    Load atom coordinates and residue information from an mmCIF/CIF file.

    This uses MMCIF2Dict instead of Bio.PDB.MMCIFParser because ESMFold2
    output may omit fields like _atom_site.occupancy.
    """
    cif_path = Path(cif_path)

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    mmcif = MMCIF2Dict(str(cif_path))

    required_keys = [
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_comp_id",
        "_atom_site.auth_atom_id",
        "_atom_site.type_symbol",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    ]

    missing_keys = [key for key in required_keys if key not in mmcif]
    if missing_keys:
        raise KeyError(f"Missing required CIF fields: {missing_keys}")

    chains = as_list(mmcif["_atom_site.auth_asym_id"])
    residue_ids = as_list(mmcif["_atom_site.auth_seq_id"])
    residue_names = as_list(mmcif["_atom_site.auth_comp_id"])
    atom_names = as_list(mmcif["_atom_site.auth_atom_id"])
    elements = as_list(mmcif["_atom_site.type_symbol"])
    x_coords = as_list(mmcif["_atom_site.Cartn_x"])
    y_coords = as_list(mmcif["_atom_site.Cartn_y"])
    z_coords = as_list(mmcif["_atom_site.Cartn_z"])

    rows = []

    for chain, residue_id, residue_name, atom_name, element, x, y, z in zip(
        chains,
        residue_ids,
        residue_names,
        atom_names,
        elements,
        x_coords,
        y_coords,
        z_coords,
    ):
        if residue_id in {".", "?"}:
            continue

        rows.append(
            {
                "chain": chain,
                "residue_id": int(residue_id),
                "residue_name": residue_name,
                "atom_name": atom_name,
                "element": element,
                "x": float(x),
                "y": float(y),
                "z": float(z),
            }
        )

    return pd.DataFrame(rows)


def analyze_interface(
    cif_path: str | Path,
    output_csv: str | Path,
    chain_a: str = "A",
    chain_b: str = "B",
    cutoff: float = 5.0,
    heavy_atoms_only: bool = True,
) -> None:
    atoms = load_atoms_from_cif(cif_path)

    if heavy_atoms_only:
        atoms = atoms[atoms["element"].str.upper() != "H"].copy()

    atoms_a = atoms[atoms["chain"] == chain_a].copy()
    atoms_b = atoms[atoms["chain"] == chain_b].copy()

    if atoms_a.empty:
        raise ValueError(f"No atoms found for chain {chain_a}")

    if atoms_b.empty:
        raise ValueError(f"No atoms found for chain {chain_b}")

    coords_a = atoms_a[["x", "y", "z"]].to_numpy()
    coords_b = atoms_b[["x", "y", "z"]].to_numpy()

    contacts = []
    min_distance = float("inf")

    print("Analyzing interface...")
    print("CIF:", cif_path)
    print("Chain A:", chain_a, "atoms:", len(atoms_a))
    print("Chain B:", chain_b, "atoms:", len(atoms_b))
    print("Cutoff:", cutoff, "Å")
    print()

    for i, coord_a in enumerate(coords_a):
        distances = np.linalg.norm(coords_b - coord_a, axis=1)

        closest_j = int(np.argmin(distances))
        closest_distance = float(distances[closest_j])
        min_distance = min(min_distance, closest_distance)

        contact_indices = np.where(distances <= cutoff)[0]

        for j in contact_indices:
            atom_a = atoms_a.iloc[i]
            atom_b = atoms_b.iloc[int(j)]
            distance = float(distances[int(j)])

            contacts.append(
                {
                    "chain_a": atom_a["chain"],
                    "residue_a_id": atom_a["residue_id"],
                    "residue_a_name": atom_a["residue_name"],
                    "atom_a": atom_a["atom_name"],
                    "chain_b": atom_b["chain"],
                    "residue_b_id": atom_b["residue_id"],
                    "residue_b_name": atom_b["residue_name"],
                    "atom_b": atom_b["atom_name"],
                    "distance_angstrom": distance,
                }
            )

    contacts_df = pd.DataFrame(contacts)

    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if contacts_df.empty:
        contacts_df.to_csv(output_csv, index=False)

        print("No atom contacts found within cutoff.")
        print(f"Minimum chain-chain atom distance: {min_distance:.3f} Å")
        print()
        print("Saved empty contact table to:")
        print(output_csv)
        return

    contacts_df = contacts_df.sort_values("distance_angstrom")
    contacts_df.to_csv(output_csv, index=False)

    interface_residues_a = (
        contacts_df[["residue_a_id", "residue_a_name"]]
        .drop_duplicates()
        .sort_values("residue_a_id")
    )

    interface_residues_b = (
        contacts_df[["residue_b_id", "residue_b_name"]]
        .drop_duplicates()
        .sort_values("residue_b_id")
    )

    residue_pair_count = (
        contacts_df[["residue_a_id", "residue_b_id"]]
        .drop_duplicates()
        .shape[0]
    )

    print("Interface analysis completed.")
    print(f"Minimum chain-chain atom distance: {min_distance:.3f} Å")
    print(f"Atom contacts within {cutoff:.1f} Å:", len(contacts_df))
    print("Unique residue-residue contact pairs:", residue_pair_count)
    print("Chain A interface residues:", len(interface_residues_a))
    print("Chain B interface residues:", len(interface_residues_b))
    print()

    print("Chain A interface residues:")
    print(
        ", ".join(
            f"{row.residue_a_name}{row.residue_a_id}"
            for row in interface_residues_a.itertuples(index=False)
        )
    )
    print()

    print("Chain B interface residues:")
    print(
        ", ".join(
            f"{row.residue_b_name}{row.residue_b_id}"
            for row in interface_residues_b.itertuples(index=False)
        )
    )
    print()

    print("Saved contact table to:")
    print(output_csv)


def parse_args():
    parser = ArgumentParser(
        description="Analyze atom contacts between two chains in an ESMFold2 complex CIF."
    )

    parser.add_argument(
        "--cif",
        required=True,
        help="Path to predicted complex .cif file.",
    )

    parser.add_argument(
        "--out",
        default="outputs/interface_contacts.csv",
        help="Output CSV path for atom-level contacts.",
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
        help="Distance cutoff in Angstroms for defining atom contacts.",
    )

    parser.add_argument(
        "--include-hydrogens",
        action="store_true",
        help="Include hydrogen atoms. By default, hydrogens are excluded.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    analyze_interface(
        cif_path=args.cif,
        output_csv=args.out,
        chain_a=args.chain_a,
        chain_b=args.chain_b,
        cutoff=args.cutoff,
        heavy_atoms_only=not args.include_hydrogens,
    )


if __name__ == "__main__":
    main()