#!/usr/bin/env python3
# Juliana A

import sys
from Bio import SeqIO
#test
def calculate_gene_lengths(fasta_file):
    gene_lengths = []
    gene_length_dict = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_length = len(record.seq)
        gene_lengths.append(gene_length)
        gene_length_dict[record.id] = gene_length

    if gene_lengths:
        mean_length = sum(gene_lengths) / len(gene_lengths)
        min_length = min(gene_lengths)
        max_length = max(gene_lengths)
    else:
        mean_length = 0
        min_length = 0
        max_length = 0

    return gene_length_dict, mean_length, min_length, max_length

def main(fasta_file):
    gene_length_dict, mean_length, min_length, max_length = calculate_gene_lengths(fasta_file)

    # Output gene lengths
    print("Gene lengths:")
    for gene_id, length in gene_length_dict.items():
        print(f"{gene_id}: {length}")

    # Output statistics
    print("\nStatistics:")
    print(f"Mean gene length: {mean_length}")
    print(f"Minimum gene length: {min_length}")
    print(f"Maximum gene length: {max_length}")

    # Calculate percentages relative to mean length
    if mean_length > 0:
        min_percent = ((mean_length - min_length) / mean_length) * 100
        max_percent = ((max_length - mean_length) / mean_length) * 100
    else:
        min_percent = 0
        max_percent = 0

    # Output percentages
    print("\nPercentage relative to mean length:")
    print(f"Percentage of minimum length: {min_percent:.2f}%")
    print(f"Percentage of maximum length: {max_percent:.2f}%")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_gene_length.py PATH/TO/FASTA/FILE")
        sys.exit(1)

    fasta_file = sys.argv[1]
    main(fasta_file)

