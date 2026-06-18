from pathlib import Path

import numpy as np
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

from src.fasta import load_fasta_sequence


def main():
    print("Testing ESMC model loading and inference...")
    print()

    model_name = "biohub/ESMC-6B"
    fasta_path = "data/sequences/chain_a.fasta"

    output_dir = Path("outputs/esmc_embeddings")
    output_dir.mkdir(parents=True, exist_ok=True)

    header, sequence = load_fasta_sequence(fasta_path)

    print("Loaded sequence:")
    print("Header:", header)
    print("Length:", len(sequence))
    print()

    print("Loading tokenizer...")
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    print("Tokenizing sequence...")
    inputs = tokenizer([sequence], return_tensors="pt", padding=True)

    print("Loading model:")
    print(model_name)
    print()

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

    print()
    print("Inference completed successfully.")

    logits = output.logits.detach().float().cpu().numpy()
    hidden_states = output.last_hidden_state.detach().float().cpu().numpy()

    print("Logits shape:", logits.shape)
    print("Last hidden state shape:", hidden_states.shape)

    token_embeddings = hidden_states[0]

    residue_embeddings = token_embeddings[1 : len(sequence) + 1]

    protein_embedding = residue_embeddings.mean(axis=0)

    token_embedding_path = output_dir / f"{header}_token_embeddings.npy"
    protein_embedding_path = output_dir / f"{header}_protein_embedding.npy"

    np.save(token_embedding_path, residue_embeddings)
    np.save(protein_embedding_path, protein_embedding)

    print()
    print("Saved residue-level embeddings to:")
    print(token_embedding_path)
    print("Shape:", residue_embeddings.shape)

    print()
    print("Saved protein-level embedding to:")
    print(protein_embedding_path)
    print("Shape:", protein_embedding.shape)


if __name__ == "__main__":
    main()