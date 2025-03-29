import time
import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
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
def GRPH(file_name, k=3):
    data = read_fasta(file_name,)  # input file
    print(data)
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    G= nx.DiGraph()
    for s_name, s in data.items():
        # print("s "+s_name + " "+s)
        for t_name, t in data.items():
            if t_name != s_name and s[len(s)-k:len(s)] == t[0:k ]:
                # print("t ", t_name + " " + t)
                G.add_edge(s_name, t_name)
            #list(G.all_neighbors( s_name))
    print(G.edges())
    for edge in G.edges():
        print(edge[0] +" "+edge[1])
        fout.write(edge[0] +" "+edge[1]+ "\n")
    # nx.draw_networkx(G, arrows=True,with_labels = True)
    # plt.show()



start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt'
input_fname= 'rosalind_grph.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file


GRPH(input_fname) #0.0004248619079589844

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




