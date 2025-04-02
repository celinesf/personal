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
def find_motif(motif=[],seq=''):
    # print(seq)
    new_motif=[]
    max_len=0
    # print(motif)
    for mi,m in enumerate(motif):
        match = SequenceMatcher(None, m, seq).find_longest_match()
        print(match)
    #     print('here1 '+str(mi)+" "+ m)
    #     print(match)
    #     for block in match:
    #         if block.size > 0 and max_len <= block.size:
    #             max_len = block.size
    #             print('here2')
    #             print(block)
    #             print(seq[block.b:block.b + block.size])
    #             print(m[block.a:block.a + block.size])
    #             new_motif.append(seq[block.b:block.b + block.size])
    #             print(max_len)
    #             print( new_motif)
    # new_motif = [s for s in new_motif if len(s) == max_len]
    # print(max_len)
    # print(new_motif)
    return new_motif

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
    return (list(set(motifs)), longest_length)

def LCSM(file_name):

    data = read_fasta(file_name, )  # input file
    # print(data)
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')

    motif = []

    for index ,(seq_name, seq) in list(enumerate(data.items())):
        print(str(index)+" "+seq_name )#+' '+ seq + ' '+ str(len(seq)))
        if index == 0:
            motif.append(seq)
        elif len(motif) > 0:
            new_motif=[] # list of new found motifs
            max_len = len(motif[0]) # previous motifs lenghts
            new_max = 0 # max lenght of new motifs
            for mi,m in enumerate(motif):
                # print('here')
                # print(m)
                tp_motif, new_len =find_longest_common_motifs(m, seq)
                # print('here2')
                # print(tp_motif)
                # print(new_len)
                if new_len <= max_len:# if new motif lenght <= old
                    if new_len >= new_max:# add motif if  larger than previously found
                        new_max = new_len
                        for n in tp_motif:
                            new_motif.append(n)

                # print('here3')
                # print(new_motif)
                # print(new_len)
                # print(new_max)
            motif = [s for s in new_motif if len(s) == new_max]
            # print('here4')
            # print(motif)
            # print(max_len)
    print(motif)
    if(len(motif)>0):
        fout.write(motif[0]+ "\n")


start_time = time.time() # report time of function

input_fname= 'Sample_Dataset4.txt'
input_fname= 'rosalind_lcsm (1).txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file

LCSM(input_fname) #

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




