import os
import matplotlib.pyplot as plt
from collections import Counter
from Bio import SeqIO

# Function to read FASTA file


def read_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# Function to count nucleotide frequencies


def count_nucleotides(sequence):
    return Counter(sequence)

# Function to compute GC content


def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

# Function to get reverse complement


def reverse_complement(sequence):
    complement_map = str.maketrans("ATCG", "TAGC")
    return sequence[::-1].translate(complement_map)

# Function to transcribe DNA to RNA


def transcribe(sequence):
    return sequence.replace('T', 'U')

# Function to visualize nucleotide distribution


def plot_nucleotide_distribution(sequence):
    counts = count_nucleotides(sequence)
    plt.bar(counts.keys(), counts.values(), color=[
            'blue', 'green', 'red', 'orange'])
    plt.xlabel("Nucleotide")
    plt.ylabel("Frequency")
    plt.title("Nucleotide Distribution in DNA Sequence")
    plt.show()

# Main function to execute analysis


def analyze_fasta(file_path):
    sequences = read_fasta(file_path)

    for i, seq in enumerate(sequences):
        print(f"\nSequence {i+1}:")
        print(f"Length: {len(seq)}")
        print(f"Nucleotide counts: {count_nucleotides(seq)}")
        print(f"GC Content: {gc_content(seq):.2f}%")
        print(f"Reverse Complement: {reverse_complement(seq)}")
        print(f"Transcribed RNA: {transcribe(seq)}")

        plot_nucleotide_distribution(seq)


# Example usage
if __name__ == "__main__":
    fasta_file = "rcsb_pdb_6BY7.fasta"
    if os.path.exists(fasta_file):
        analyze_fasta(fasta_file)
    else:
        print("FASTA file not found. Please provide a valid file path.")
