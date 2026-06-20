# ESM Protein Binder Prediction

A small project for experimenting with Biohub ESM models for protein embedding generation, protein structure prediction, and eventually protein-protein binder/interface analysis.

The long-term goal is to explore whether two protein sequences form a plausible protein-protein complex and to identify possible interface residues between them.

This project is currently in an early prototype stage.

---

## Current Status

This project currently supports:

* Creating a working Python environment for Biohub ESM models
* Loading protein sequences from FASTA files
* Running ESMC-6B inference locally
* Saving residue-level and protein-level embeddings
* Comparing two protein-level embeddings using cosine similarity
* Running ESMFold2 single-chain structure prediction
* Saving predicted structures as `.cif` files

The project does not yet support:

* Full protein-protein complex prediction
* Interface residue detection
* Binder ranking
* Docking score calculation
* Validated biological binding prediction

---

## Current Mini Pipeline

```text
FASTA sequence
    ↓
Load protein sequence
    ↓
Run ESMC-6B
    ↓
Save residue-level embedding
    ↓
Save protein-level embedding
    ↓
Compare protein embeddings
```

The structure prediction pipeline is:

```text
FASTA sequence
    ↓
Load protein sequence
    ↓
Run ESMFold2
    ↓
Generate predicted structure
    ↓
Save .cif structure file
```

---

## Folder Structure

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
│   ├── esmc_embeddings/
│   └── structures/
├── scripts/
│   ├── compare_embeddings.py
│   ├── embed_sequence.py
│   ├── predict_single_structure.py
│   ├── test_esmc_model.py
│   ├── test_fasta_loading.py
│   └── test_install.py
├── src/
│   ├── __init__.py
│   └── fasta.py
└── notebooks/
```

---

## Installation

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

---

## Environment Notes

This project currently uses:

```text
Python 3.12
PyTorch
Biohub esm package
Hugging Face model loading
Biopython
NumPy
Pandas
Matplotlib
```

On Apple Silicon Macs, PyTorch may show:

```text
MPS available: True
```

This means PyTorch can use Apple Metal acceleration.

---

## Hugging Face Access

The Biohub ESM models are loaded through Hugging Face.

If needed, log in with:

```bash
huggingface-cli login
```

or:

```bash
hf auth login
```

You may need a Hugging Face account and model access before downloading model weights.

---

## Test Installation

Run:

```bash
python scripts/test_install.py
```

A successful test should show that both `torch` and `esm` import correctly.

Example successful output:

```text
Testing ESM Protein Binder Prediction setup...

torch imported successfully
torch version: ...
CUDA available: False
MPS available: True

esm imported successfully
esm location: ...

Setup test completed successfully.
```

On Mac, it is normal for CUDA to be unavailable:

```text
CUDA available: False
```

---

## Test FASTA Loading

Run:

```bash
python scripts/test_fasta_loading.py
```

This loads the example sequences from:

```text
data/sequences/chain_a.fasta
data/sequences/chain_b.fasta
```

Expected output includes the FASTA headers, sequence lengths, and sequences.

---

## Example FASTA Files

Current example files:

```text
data/sequences/chain_a.fasta
data/sequences/chain_b.fasta
```

Each file should contain one protein sequence.

Example format:

```fasta
>chain_a_example
MKKAVINGEQIRSISDLHQTLKKELALPEYYGENLDALWDCLTGWVEYPLVLEWRQFEQSKQLTENGAESVLQVFREAKAEGADITIILS
```

---

## Run a First ESMC Model Test

To run the basic ESMC test script:

```bash
caffeinate -dis python scripts/test_esmc_model.py
```

The `caffeinate` command prevents the Mac from sleeping while the model is running.

This script currently:

* Loads the Chain A example sequence
* Loads `biohub/ESMC-6B`
* Runs inference
* Saves residue-level embeddings
* Saves a protein-level embedding

Expected output files:

```text
outputs/esmc_embeddings/chain_a_example_token_embeddings.npy
outputs/esmc_embeddings/chain_a_example_protein_embedding.npy
```

---

## Generate ESMC Embeddings for Any FASTA File

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

By default, outputs are saved to:

```text
outputs/esmc_embeddings/
```

The script saves two embedding files:

```text
*_residue_embeddings.npy
*_protein_embedding.npy
```

The residue-level embedding has one vector per amino acid residue.

The protein-level embedding is the average of all residue embeddings and gives one vector for the whole protein.

Example output shapes:

```text
Residue-level embedding: (90, 2560)
Protein-level embedding: (2560,)
```

---

## Compare Two Protein Embeddings

After generating embeddings for Chain A and Chain B, run:

```bash
python scripts/compare_embeddings.py \
  --embedding-a outputs/esmc_embeddings/chain_a_example_protein_embedding.npy \
  --embedding-b outputs/esmc_embeddings/chain_b_example_protein_embedding.npy
```

Example output:

```text
Embedding A: outputs/esmc_embeddings/chain_a_example_protein_embedding.npy
Shape A: (2560,)

Embedding B: outputs/esmc_embeddings/chain_b_example_protein_embedding.npy
Shape B: (2560,)

Cosine similarity: 0.6095994114875793
```

Cosine similarity ranges from `-1` to `1`.

```text
1.0   very similar embedding direction
0.0   unrelated or orthogonal direction
-1.0  opposite embedding direction
```

For the current toy example, the exact biological meaning is not important yet. The goal is simply to confirm that the embedding pipeline works.

---

## Predict a Single Protein Structure with ESMFold2

This project can also run a first single-chain ESMFold2 structure prediction.

On Apple Silicon, ESMFold2 may require PyTorch MPS fallback because some operations are not implemented on MPS yet.

Run:

```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 caffeinate -dis python scripts/predict_single_structure.py \
  --fasta data/sequences/chain_a.fasta \
  --out outputs/structures/chain_a_example.cif
```

Expected output:

```text
Running ESMFold2 single-chain structure prediction...
FASTA: data/sequences/chain_a.fasta
Header: chain_a_example
Length: 90
Output: outputs/structures/chain_a_example.cif

Using device: mps

Running folding...
num_loops: 4
num_sampling_steps: 20

Folding completed.
pLDDT mean: ...
pTM: ...
ipTM: ...

Saved predicted structure to:
outputs/structures/chain_a_example.cif
```

For single-chain prediction, `ipTM` is expected to be `0.000` because there is no inter-chain interface.

---

## Predict Chain B Structure

Run:

```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 caffeinate -dis python scripts/predict_single_structure.py \
  --fasta data/sequences/chain_b.fasta \
  --out outputs/structures/chain_b_example.cif
```

Expected output file:

```text
outputs/structures/chain_b_example.cif
```

---

## Viewing Predicted Structures

Predicted `.cif` files can be opened in molecular visualization tools such as:

* ChimeraX
* PyMOL
* Mol*
* UCSF Chimera

On macOS, you can try:

```bash
open outputs/structures/chain_a_example.cif
```

or open the file manually from ChimeraX/PyMOL.

---

## Current Scripts

### `scripts/test_install.py`

Checks that the environment is working and confirms that `torch` and `esm` import correctly.

### `scripts/test_fasta_loading.py`

Loads the example FASTA files and prints their headers, lengths, and sequences.

### `scripts/test_esmc_model.py`

Runs a basic ESMC-6B inference test on Chain A and saves embeddings.

### `scripts/embed_sequence.py`

Reusable script for generating ESMC embeddings from any single-sequence FASTA file.

Example:

```bash
caffeinate -dis python scripts/embed_sequence.py \
  --fasta data/sequences/chain_a.fasta
```

### `scripts/compare_embeddings.py`

Compares two protein-level embeddings using cosine similarity.

Example:

```bash
python scripts/compare_embeddings.py \
  --embedding-a outputs/esmc_embeddings/chain_a_example_protein_embedding.npy \
  --embedding-b outputs/esmc_embeddings/chain_b_example_protein_embedding.npy
```

### `scripts/predict_single_structure.py`

Runs ESMFold2 single-chain structure prediction and saves a `.cif` file.

Example:

```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 caffeinate -dis python scripts/predict_single_structure.py \
  --fasta data/sequences/chain_a.fasta \
  --out outputs/structures/chain_a_example.cif
```

---

## Common Warnings

Some warnings may appear when loading ESMC or ESMFold2 locally.

The following warnings are okay for now:

```text
torchvision ... libjpeg
```

This is related to image loading support in `torchvision`. This project is not currently using image functions.

```text
transformer_engine is not installed
```

This means optimized fused operations are unavailable. The model falls back to regular PyTorch.

```text
neither xformers nor flash-attn is installed
```

This means optimized attention kernels are unavailable. The model still runs with PyTorch attention.

```text
No checkpoint found for ESMCForSequenceClassification.forward
```

This appears to be a model/package documentation warning and does not block the current workflow.

---

## Apple Silicon / MPS Notes

For ESMC embedding generation, MPS works well on this setup.

For ESMFold2 structure prediction, MPS may need fallback:

```bash
PYTORCH_ENABLE_MPS_FALLBACK=1
```

This lets unsupported MPS operations fall back to CPU instead of crashing.

Without this, ESMFold2 may fail with an error like:

```text
NotImplementedError: The operator 'aten::scatter_reduce.two_out' is not currently implemented for the MPS device.
```

---

## Git / Output Notes

Generated model outputs are ignored by Git.

These folders and file types should not be committed:

```text
outputs/
*.npy
*.npz
*.pdb
*.cif
*.mmcif
*.pt
*.pth
```

This is because model outputs and structure files can become large or frequently regenerated.

---

## Current Working Commands

Full current workflow:

```bash
conda activate esm-binder
```

Test installation:

```bash
python scripts/test_install.py
```

Test FASTA loading:

```bash
python scripts/test_fasta_loading.py
```

Generate Chain A embedding:

```bash
caffeinate -dis python scripts/embed_sequence.py \
  --fasta data/sequences/chain_a.fasta
```

Generate Chain B embedding:

```bash
caffeinate -dis python scripts/embed_sequence.py \
  --fasta data/sequences/chain_b.fasta
```

Compare embeddings:

```bash
python scripts/compare_embeddings.py \
  --embedding-a outputs/esmc_embeddings/chain_a_example_protein_embedding.npy \
  --embedding-b outputs/esmc_embeddings/chain_b_example_protein_embedding.npy
```

Predict Chain A structure:

```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 caffeinate -dis python scripts/predict_single_structure.py \
  --fasta data/sequences/chain_a.fasta \
  --out outputs/structures/chain_a_example.cif
```

Predict Chain B structure:

```bash
PYTORCH_ENABLE_MPS_FALLBACK=1 caffeinate -dis python scripts/predict_single_structure.py \
  --fasta data/sequences/chain_b.fasta \
  --out outputs/structures/chain_b_example.cif
```

---

## Planned Next Steps

The next planned steps are:

1. Predict structures for both example chains.
2. Add a two-chain ESMFold2 complex prediction script.
3. Save predicted complex structures as `.cif` files.
4. Identify residues close to the opposite chain.
5. Generate a simple interface residue table.
6. Generate a contact map between Chain A and Chain B.
7. Eventually rank possible binder candidates against a target protein.

Longer-term possible directions:

* Alpha-synuclein binder exploration
* KLK6 binder exploration
* Mutation effect analysis using ESMC embeddings
* Protein-switch design sandbox
* Comparison of ESMFold2 complex predictions with docking tools
