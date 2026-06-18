# ESM-Protein-Binder-Prediction

A small project for experimenting with Biohub ESM models for protein-protein binder prediction, protein embedding generation, and eventually protein-protein interface analysis.

## Current status

This project currently supports:

1. Creating a working Python environment for Biohub ESM models.
2. Loading protein sequences from FASTA files.
3. Running ESMC-6B inference locally.
4. Saving residue-level and protein-level embeddings.
5. Comparing two protein-level embeddings using cosine similarity.

## Project structure

```text
ESM-Protein-Binder-Prediction/
├── README.md
├── .gitignore
├── environment.yml
├── data/
│   └── sequences/
│       ├── chain_a.fasta
│       └── chain_b.fasta
├── outputs/
│   └── esmc_embeddings/
├── scripts/
│   ├── compare_embeddings.py
│   ├── embed_sequence.py
│   ├── test_esmc_model.py
│   ├── test_fasta_loading.py
│   └── test_install.py
├── src/
│   ├── __init__.py
│   └── fasta.py
└── notebooks/
```

## Setup

Create the conda environment:

```bash
conda env create -f environment.yml
```

Activate the environment:

```bash
conda activate esm-binder
```

Add the environment to Jupyter:

```bash
python -m ipykernel install --user --name esm-binder --display-name "Python (esm-binder)"
```

## Test the installation

Run:

```bash
python scripts/test_install.py
```

A successful test should show that `torch` and `esm` both import correctly.

On Apple Silicon Macs, it is good if you see:

```text
MPS available: True
```

This means PyTorch can use Apple’s Metal acceleration.

## Test FASTA loading

Run:

```bash
python scripts/test_fasta_loading.py
```

This loads the example sequences from:

```text
data/sequences/chain_a.fasta
data/sequences/chain_b.fasta
```

## Run a first ESMC model test

To prevent your Mac from sleeping while the model runs, use `caffeinate`:

```bash
caffeinate -dis python scripts/test_esmc_model.py
```

This script loads `biohub/ESMC-6B`, runs inference on the example Chain A sequence, and saves embeddings.

Expected output files:

```text
outputs/esmc_embeddings/chain_a_example_token_embeddings.npy
outputs/esmc_embeddings/chain_a_example_protein_embedding.npy
```

## Generate embeddings for any FASTA file

Use:

```bash
caffeinate -dis python scripts/embed_sequence.py \
  --fasta data/sequences/chain_a.fasta
```

For Chain B:

```bash
caffeinate -dis python scripts/embed_sequence.py \
  --fasta data/sequences/chain_b.fasta
```

This saves two types of embeddings:

```text
*_residue_embeddings.npy
*_protein_embedding.npy
```

The residue-level embedding has one vector per amino acid residue.

The protein-level embedding is the average of all residue embeddings and gives one vector for the whole protein.

## Compare two protein embeddings

After generating embeddings for Chain A and Chain B, run:

```bash
python scripts/compare_embeddings.py \
  --embedding-a outputs/esmc_embeddings/chain_a_example_protein_embedding.npy \
  --embedding-b outputs/esmc_embeddings/chain_b_example_protein_embedding.npy
```

Example output:

```text
Cosine similarity: 0.6095994114875793
```

Cosine similarity ranges from `-1` to `1`.

```text
1.0   very similar embedding direction
0.0   unrelated or orthogonal direction
-1.0  opposite embedding direction
```

## Notes

Some warnings may appear when loading ESMC-6B locally.

The following warnings are okay for now:

- `torchvision ... libjpeg`
- `transformer_engine is not installed`
- `xformers nor flash-attn is installed`
- `No checkpoint found... auto_docstring`

These are not blocking the current workflow.

## Current mini pipeline

```text
FASTA sequence
→ ESMC-6B inference
→ residue-level embedding
→ protein-level embedding
→ cosine similarity comparison
```

## Planned next steps

1. Clean up the embedding scripts into reusable utilities.
2. Add support for embedding multiple FASTA files at once.
3. Try a first ESMFold2 structure prediction.
4. Try a first protein-protein complex prediction.
5. Analyze predicted interface residues between two chains.
6. Eventually rank candidate binders against a target protein.