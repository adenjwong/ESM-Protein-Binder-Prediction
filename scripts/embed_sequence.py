from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

from src.fasta import load_fasta_sequence


def embed_sequence(
    fasta_path: str | Path,
    output_dir: str | Path,
    model_name: str = "biohub/ESMC-6B",
) -> None:
    fasta_path = Path(fasta_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Running ESMC embedding generation...")
    print("FASTA:", fasta_path)
    print("Model:", model_name)
    print("Output directory:", output_dir)
    print()

    header, sequence = load_fasta_sequence(fasta_path)

    print("Loaded sequence:")
    print("Header:", header)
    print("Length:", len(sequence))
    print()

    print("Loading tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    print("Tokenizing sequence...")
    inputs = tokenizer([sequence], return_tensors="pt", padding=True)

    print("Loading model...")
    model = AutoModelForMaskedLM.from_pretrained(
        model_name,
        device_map="auto",
    ).eval()

    device = next(model.parameters()).device
    print("Model device:", device)

    inputs = {key: value.to(device) for key, value in inputs.items()}

    print("Running inference...")
    with torch.inference_mode():
        output = model(**inputs)

    print("Inference completed successfully.")
    print()

    logits = output.logits.detach().float().cpu().numpy()
    hidden_states = output.last_hidden_state.detach().float().cpu().numpy()

    print("Logits shape:", logits.shape)
    print("Last hidden state shape:", hidden_states.shape)

    token_embeddings = hidden_states[0]

    residue_embeddings = token_embeddings[1 : len(sequence) + 1]

    protein_embedding = residue_embeddings.mean(axis=0)

    safe_header = header.replace("/", "_").replace(" ", "_")

    residue_embedding_path = output_dir / f"{safe_header}_residue_embeddings.npy"
    protein_embedding_path = output_dir / f"{safe_header}_protein_embedding.npy"

    np.save(residue_embedding_path, residue_embeddings)
    np.save(protein_embedding_path, protein_embedding)

    print()
    print("Saved residue-level embeddings to:")
    print(residue_embedding_path)
    print("Shape:", residue_embeddings.shape)

    print()
    print("Saved protein-level embedding to:")
    print(protein_embedding_path)
    print("Shape:", protein_embedding.shape)


def parse_args():
    parser = ArgumentParser(
        description="Generate ESMC embeddings for a protein sequence in FASTA format."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to input FASTA file containing one protein sequence.",
    )

    parser.add_argument(
        "--out",
        default="outputs/esmc_embeddings",
        help="Directory where embedding files should be saved.",
    )

    parser.add_argument(
        "--model",
        default="biohub/ESMC-6B",
        help="Hugging Face model name to use.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    embed_sequence(
        fasta_path=args.fasta,
        output_dir=args.out,
        model_name=args.model,
    )


if __name__ == "__main__":
    main()