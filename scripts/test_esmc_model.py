from src.fasta import load_fasta_sequence


def main():
    print("Testing ESMC model loading and inference...")
    print()

    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer

    model_name = "biohub/ESMC-6B"

    header, sequence = load_fasta_sequence("data/sequences/chain_a.fasta")

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

    if hasattr(output, "logits"):
        print("Logits shape:", tuple(output.logits.shape))

    if hasattr(output, "hidden_states") and output.hidden_states is not None:
        print("Hidden states:", len(output.hidden_states))

    if hasattr(output, "last_hidden_state"):
        print("Last hidden state shape:", tuple(output.last_hidden_state.shape))


if __name__ == "__main__":
    main()