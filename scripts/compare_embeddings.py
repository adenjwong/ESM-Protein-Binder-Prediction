from argparse import ArgumentParser
from pathlib import Path

import numpy as np


def cosine_similarity(vector_a: np.ndarray, vector_b: np.ndarray) -> float:
    """
    Compute cosine similarity between two vectors.

    Returns a value between -1 and 1:
        1  = very similar direction
        0  = unrelated/orthogonal direction
       -1  = opposite direction
    """
    numerator = np.dot(vector_a, vector_b)
    denominator = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)

    if denominator == 0:
        raise ValueError("Cannot compute cosine similarity with a zero vector.")

    return float(numerator / denominator)


def load_embedding(path: str | Path) -> np.ndarray:
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Embedding file not found: {path}")

    embedding = np.load(path)

    if embedding.ndim != 1:
        raise ValueError(
            f"Expected a protein-level embedding with shape (features,), "
            f"but got shape {embedding.shape} from {path}"
        )

    return embedding


def parse_args():
    parser = ArgumentParser(
        description="Compare two protein-level ESMC embeddings using cosine similarity."
    )

    parser.add_argument(
        "--embedding-a",
        required=True,
        help="Path to first protein-level embedding .npy file.",
    )

    parser.add_argument(
        "--embedding-b",
        required=True,
        help="Path to second protein-level embedding .npy file.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    embedding_a = load_embedding(args.embedding_a)
    embedding_b = load_embedding(args.embedding_b)

    similarity = cosine_similarity(embedding_a, embedding_b)

    print("Embedding A:", args.embedding_a)
    print("Shape A:", embedding_a.shape)
    print()
    print("Embedding B:", args.embedding_b)
    print("Shape B:", embedding_b.shape)
    print()
    print("Cosine similarity:", similarity)


if __name__ == "__main__":
    main()