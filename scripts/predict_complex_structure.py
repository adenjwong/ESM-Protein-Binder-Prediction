from argparse import ArgumentParser
from pathlib import Path

import torch

from esm.models.esmfold2 import (
    ESMFold2InputBuilder,
    ProteinInput,
    StructurePredictionInput,
)
from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

from src.fasta import load_fasta_sequence


def get_device() -> str:
    """
    Choose the best available device.

    On Apple Silicon, this should return 'mps'.
    On NVIDIA systems, this should return 'cuda'.
    Otherwise, it returns 'cpu'.
    """
    if torch.cuda.is_available():
        return "cuda"

    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return "mps"

    return "cpu"


def predict_complex_structure(
    chain_a_fasta: str | Path,
    chain_b_fasta: str | Path,
    output_path: str | Path,
    num_loops: int = 4,
    num_sampling_steps: int = 20,
    seed: int = 0,
) -> None:
    chain_a_fasta = Path(chain_a_fasta)
    chain_b_fasta = Path(chain_b_fasta)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    chain_a_header, chain_a_sequence = load_fasta_sequence(chain_a_fasta)
    chain_b_header, chain_b_sequence = load_fasta_sequence(chain_b_fasta)

    print("Running ESMFold2 two-chain complex prediction...")
    print("Chain A FASTA:", chain_a_fasta)
    print("Chain A header:", chain_a_header)
    print("Chain A length:", len(chain_a_sequence))
    print()
    print("Chain B FASTA:", chain_b_fasta)
    print("Chain B header:", chain_b_header)
    print("Chain B length:", len(chain_b_sequence))
    print()
    print("Output:", output_path)
    print()

    device = get_device()
    print("Using device:", device)
    print()

    print("Loading ESMFold2 model...")
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").eval()

    if device == "cuda":
        model = model.cuda()
    elif device == "mps":
        # ESMFold2 mostly runs on Apple Silicon MPS, but some PyTorch ops are not
        # implemented on MPS yet. Run with PYTORCH_ENABLE_MPS_FALLBACK=1 so those
        # unsupported ops can fall back to CPU instead of crashing.
        model = model.to("mps")
    else:
        model = model.cpu()

    structure_input = StructurePredictionInput(
        sequences=[
            ProteinInput(id="A", sequence=chain_a_sequence),
            ProteinInput(id="B", sequence=chain_b_sequence),
        ]
    )

    print("Running complex folding...")
    print("num_loops:", num_loops)
    print("num_sampling_steps:", num_sampling_steps)
    print()

    result = ESMFold2InputBuilder().fold(
        model,
        structure_input,
        num_loops=num_loops,
        num_sampling_steps=num_sampling_steps,
        num_diffusion_samples=1,
        seed=seed,
    )

    print("Complex folding completed.")
    print(f"pLDDT mean: {float(result.plddt.mean()):.3f}")
    print(f"pTM: {float(result.ptm):.3f}")
    print(f"ipTM: {float(result.iptm):.3f}")

    with open(output_path, "w") as f:
        f.write(result.complex.to_mmcif())

    print()
    print("Saved predicted complex structure to:")
    print(output_path)


def parse_args():
    parser = ArgumentParser(
        description="Predict a two-chain protein complex using ESMFold2."
    )

    parser.add_argument(
        "--chain-a",
        required=True,
        help="Path to FASTA file for Chain A.",
    )

    parser.add_argument(
        "--chain-b",
        required=True,
        help="Path to FASTA file for Chain B.",
    )

    parser.add_argument(
        "--out",
        default="outputs/structures/predicted_complex.cif",
        help="Output .cif file path.",
    )

    parser.add_argument(
        "--num-loops",
        type=int,
        default=4,
        help="Number of ESMFold2 refinement loops. Lower is faster.",
    )

    parser.add_argument(
        "--num-sampling-steps",
        type=int,
        default=20,
        help="Number of diffusion sampling steps. Lower is faster.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for complex prediction.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    predict_complex_structure(
        chain_a_fasta=args.chain_a,
        chain_b_fasta=args.chain_b,
        output_path=args.out,
        num_loops=args.num_loops,
        num_sampling_steps=args.num_sampling_steps,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()