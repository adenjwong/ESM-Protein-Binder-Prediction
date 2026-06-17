from src.fasta import load_fasta_sequence


def main():
    chain_a_header, chain_a_sequence = load_fasta_sequence("data/sequences/chain_a.fasta")
    chain_b_header, chain_b_sequence = load_fasta_sequence("data/sequences/chain_b.fasta")

    print("Loaded Chain A")
    print("Header:", chain_a_header)
    print("Length:", len(chain_a_sequence))
    print("Sequence:", chain_a_sequence)
    print()

    print("Loaded Chain B")
    print("Header:", chain_b_header)
    print("Length:", len(chain_b_sequence))
    print("Sequence:", chain_b_sequence)


if __name__ == "__main__":
    main()