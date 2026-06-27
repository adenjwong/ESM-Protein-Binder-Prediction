import os

# Needed on Apple Silicon because some ESMFold2 operations are not implemented on MPS.
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

from argparse import ArgumentParser
from pathlib import Path

from scripts.predict_complex_structure import predict_complex_structure
from scripts.run_interface_workflow import run_interface_workflow


def run_seed_sweep(
    chain_a: str | Path,
    chain_b: str | Path,
    seeds: list[int],
    output_dir: str | Path = "outputs/seed_sweep",
    num_loops: int = 4,
    num_sampling_steps: int = 20,
) -> None:
    output_dir = Path(output_dir)
    structures_dir = output_dir / "structures"
    interface_dir = output_dir / "interface_contacts"

    structures_dir.mkdir(parents=True, exist_ok=True)
    interface_dir.mkdir(parents=True, exist_ok=True)

    print("Running complex prediction seed sweep")
    print("Chain A:", chain_a)
    print("Chain B:", chain_b)
    print("Seeds:", seeds)
    print("Output directory:", output_dir)
    print()

    for seed in seeds:
        name = f"complex_seed_{seed}"
        complex_cif = structures_dir / f"{name}.cif"

        print("=" * 80)
        print(f"Predicting complex for seed {seed}")
        print("=" * 80)

        predict_complex_structure(
            chain_a_fasta=chain_a,
            chain_b_fasta=chain_b,
            output_path=complex_cif,
            num_loops=num_loops,
            num_sampling_steps=num_sampling_steps,
            seed=seed,
        )

        print()
        print("=" * 80)
        print(f"Analyzing interface for seed {seed}")
        print("=" * 80)

        run_interface_workflow(
            cif_path=complex_cif,
            name=name,
            output_dir=interface_dir,
            chain_a="A",
            chain_b="B",
            cutoff=5.0,
        )

    print()
    print("Seed sweep completed.")
    print("Structures saved to:", structures_dir)
    print("Interface outputs saved to:", interface_dir)


def parse_args():
    parser = ArgumentParser(
        description="Run multiple ESMFold2 complex predictions with different seeds."
    )

    parser.add_argument(
        "--chain-a",
        required=True,
        help="Path to Chain A FASTA file.",
    )

    parser.add_argument(
        "--chain-b",
        required=True,
        help="Path to Chain B FASTA file.",
    )

    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=[0, 1, 2],
        help="Random seeds to run.",
    )

    parser.add_argument(
        "--out-dir",
        default="outputs/seed_sweep",
        help="Output directory for seed sweep results.",
    )

    parser.add_argument(
        "--num-loops",
        type=int,
        default=4,
        help="Number of ESMFold2 refinement loops.",
    )

    parser.add_argument(
        "--num-sampling-steps",
        type=int,
        default=20,
        help="Number of diffusion sampling steps.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    run_seed_sweep(
        chain_a=args.chain_a,
        chain_b=args.chain_b,
        seeds=args.seeds,
        output_dir=args.out_dir,
        num_loops=args.num_loops,
        num_sampling_steps=args.num_sampling_steps,
    )


if __name__ == "__main__":
    main()