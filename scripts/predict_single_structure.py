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


def predict_structure(
    fasta_path: str | Path,
    output_path: str | Path,
    num_loops: int = 4,
    num_sampling_steps: int = 20,
    seed: int = 0,
) -> None:
    fasta_path = Path(fasta_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    header, sequence = load_fasta_sequence(fasta_path)

    print("Running ESMFold2 single-chain structure prediction...")
    print("FASTA:", fasta_path)
    print("Header:", header)
    print("Length:", len(sequence))
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
            ProteinInput(id="A", sequence=sequence),
        ]
    )

    print("Running folding...")
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

    print("Folding completed.")
    print(f"pLDDT mean: {float(result.plddt.mean()):.3f}")
    print(f"pTM: {float(result.ptm):.3f}")
    print(f"ipTM: {float(result.iptm):.3f}")

    with open(output_path, "w") as f:
        f.write(result.complex.to_mmcif())

    print()
    print("Saved predicted structure to:")
    print(output_path)


def parse_args():
    parser = ArgumentParser(
        description="Predict a single-chain protein structure using ESMFold2."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to input FASTA file containing one protein sequence.",
    )

    parser.add_argument(
        "--out",
        default="outputs/structures/single_chain_prediction.cif",
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
        help="Random seed for structure prediction.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    predict_structure(
        fasta_path=args.fasta,
        output_path=args.out,
        num_loops=args.num_loops,
        num_sampling_steps=args.num_sampling_steps,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()