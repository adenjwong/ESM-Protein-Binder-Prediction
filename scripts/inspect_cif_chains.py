from argparse import ArgumentParser
from pathlib import Path
from collections import defaultdict

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def as_list(value):
    """
    MMCIF2Dict sometimes returns a string if there is one value,
    or a list if there are many. Normalize everything to a list.
    """
    if isinstance(value, list):
        return value
    return [value]


def inspect_cif(cif_path: str | Path) -> None:
    cif_path = Path(cif_path)

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    mmcif = MMCIF2Dict(str(cif_path))

    print("Inspecting:", cif_path)
    print()

    print("Available atom-site fields:")
    atom_site_fields = sorted(key for key in mmcif.keys() if key.startswith("_atom_site."))
    for field in atom_site_fields:
        print(" ", field)
    print()

    chain_key = None
    for candidate_key in [
        "_atom_site.auth_asym_id",
        "_atom_site.label_asym_id",
        "_atom_site.auth_chain_id",
    ]:
        if candidate_key in mmcif:
            chain_key = candidate_key
            break

    if chain_key is None:
        raise KeyError("Could not find a chain ID field in the CIF file.")

    residue_key = None
    for candidate_key in [
        "_atom_site.auth_seq_id",
        "_atom_site.label_seq_id",
    ]:
        if candidate_key in mmcif:
            residue_key = candidate_key
            break

    if residue_key is None:
        raise KeyError("Could not find a residue index field in the CIF file.")

    residue_name_key = None
    for candidate_key in [
        "_atom_site.label_comp_id",
        "_atom_site.auth_comp_id",
    ]:
        if candidate_key in mmcif:
            residue_name_key = candidate_key
            break

    chains = as_list(mmcif[chain_key])
    residue_ids = as_list(mmcif[residue_key])
    residue_names = as_list(mmcif[residue_name_key]) if residue_name_key else ["UNK"] * len(chains)

    chain_to_residues = defaultdict(dict)

    for chain_id, residue_id, residue_name in zip(chains, residue_ids, residue_names):
        if residue_id in {".", "?"}:
            continue

        chain_to_residues[chain_id][residue_id] = residue_name

    print("Chain field used:", chain_key)
    print("Residue field used:", residue_key)
    print()
    print("Number of chains:", len(chain_to_residues))
    print()

    for chain_id in sorted(chain_to_residues.keys()):
        residues = chain_to_residues[chain_id]

        residue_indices = list(residues.keys())

        print(f"Chain {chain_id}")
        print("  Residues:", len(residue_indices))

        if residue_indices:
            first_id = residue_indices[0]
            last_id = residue_indices[-1]
            print("  First residue:", residues[first_id], first_id)
            print("  Last residue:", residues[last_id], last_id)

        print()


def parse_args():
    parser = ArgumentParser(
        description="Inspect chains and residue counts in an ESMFold2 mmCIF/CIF file."
    )

    parser.add_argument(
        "--cif",
        required=True,
        help="Path to .cif or .mmcif file.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    inspect_cif(args.cif)


if __name__ == "__main__":
    main()