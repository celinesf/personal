import time
import sys
from difflib import SequenceMatcher
from Bio import SeqIO

def read_fasta(filename):
    """Reads a FASTA file and returns a dictionary of sequences.

    Args:
        filename: The path to the FASTA file.

    Returns:
        A dictionary where keys are sequence IDs and values are sequences.
    """
    sequences = {}
    with open(filename, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = str(record.seq)
    return sequences




### My solution
def find_longest_common_motifs(seq1, seq2):
    """
    Finds all longest common motifs between two DNA sequences.

    Args:
        seq1: The first DNA sequence.
        seq2: The second DNA sequence.

    Returns:
        A list of longest common motifs.
    """

    matrix = [([0] * (len(seq2) + 1)) for _ in range(len(seq1) + 1)]
    longest_length = 0

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if seq1[i - 1] == seq2[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1] + 1
                longest_length = max(longest_length, matrix[i][j])

    motifs = []
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if matrix[i][j] == longest_length:
                motifs.append(seq1[i - longest_length:i])

    # Remove duplicate motifs
    return list(set(motifs))


# Example usage:
sequence1 = "AGGTABBCDEF"
sequence2 = "GGTABBCDEFGH"

common_motifs = find_longest_common_motifs(sequence1, sequence2)
print(f"Longest common motifs: {common_motifs}")  # Output: ['GGTABBCDEF']

sequence1 = "AGGTABBCDEFX"
sequence2 = "GGTAABBCDEFY"

common_motifs = find_longest_common_motifs(sequence1, sequence2)
print(f"Longest common motifs: {common_motifs}")  # Output: ['ABBCDEF']

sequence1 = "AGGTABBCDEFX"
sequence2 = "GGTAABBCDEFABBCDEFY"

common_motifs = find_longest_common_motifs(sequence1, sequence2)
print(f"Longest common motifs: {common_motifs}")  # Output: ['ABBCDEF']

sequence1 = "AGGTTCCAA"
sequence2 = "GGAACCC"

common_motifs = find_longest_common_motifs(sequence1, sequence2)
print(f"Longest common motifs: {common_motifs}")  # Output: ['ABBCDEF']

# def LCSM(file_name):
#
#     data = read_fasta(file_name, )  # input file
#     # print(data)
#     fout_name = file_name.removesuffix('.txt') + "_output.txt"
#     print(fout_name)
#     fout = open(fout_name, 'w')
#
#     motif = []
#     print(list(enumerate(data.items())))
#     for index ,(seq_name, seq) in list(enumerate(data.items())):
#         print(str(index)+" "+seq_name + ' '+ seq + ' '+ str(len(seq)))
#         if index == 0:
#             motif.append(seq)
#         else:
#             motif = find_motif(motif, seq)
#     print(motif)
#     if(len(motif)>0):
#         fout.write(motif[0]+ "\n")
#
#
# start_time = time.time() # report time of function
#
# input_fname= 'Sample_Dataset4.txt'
# # input_fname= 'rosalind_lcsm.txt'
#
# #sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file
#
# LCSM(input_fname) #
#
# end_time = time.time()
# execution_time = end_time - start_time
# print(f"Execution time: {execution_time} seconds")
#
#


