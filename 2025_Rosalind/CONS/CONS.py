import time
import sys
import numpy as np

from Bio import SeqIO

### Rosalind Solution


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
def CONS(file_name):
    data = read_fasta(file_name,)  # input file
    # print(data)
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    profile = ['A','C', 'G', 'T']
    n=0
    m=4 # A C G T
    i = 0
    for seq_id, seq in data.items():
        # print(f"Sequence ID: {seq_id}")
        # print(f"Sequence: {seq}\n")
        if n==0:
            n= len(seq)
            pmat=np.zeros((m,n), dtype=int)
        for index, nt in enumerate(list(seq)):
            # print(f"Index: {index}, Value: {nt}")
            pmat[ profile.index(nt),index] += 1
    cons = ''
    for i in np.argmax(pmat, axis=0):
        cons += (profile[i])
    print(cons)
    fout.write(cons+'\n')
    for index, nt in enumerate(list(profile)):
        print(nt+': '+ str(" ".join(map(str, pmat[index])) ) )
        fout.write(nt+': '+ str(" ".join(map(str, pmat[index])) )+ "\n")

start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt'
# input_fname= 'rosalind_cons.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file

# CONS(input_fname) #0.018140316009521484 seconds
CONS2(input_fname) #0.0004248619079589844

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




