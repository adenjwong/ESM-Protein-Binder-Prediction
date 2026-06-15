def main():
    print("Testing ESM Protein Binder Prediction setup...")
    print()

    try:
        import torch

        print("torch imported successfully")
        print("torch version:", torch.__version__)
        print("CUDA available:", torch.cuda.is_available())

        if hasattr(torch.backends, "mps"):
            print("MPS available:", torch.backends.mps.is_available())
        else:
            print("MPS available: not supported by this PyTorch build")

    except Exception as e:
        print("Failed to import torch")
        raise e

    print()

    try:
        import esm

        print("esm imported successfully")
        print("esm location:", esm.__file__)

    except Exception as e:
        print("Failed to import esm")
        raise e

    print()
    print("Setup test completed successfully.")


if __name__ == "__main__":
    main()