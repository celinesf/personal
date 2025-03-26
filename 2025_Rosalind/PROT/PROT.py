import time
import sys

def get_Codon_2_AminoAcid_table():
    fin = open('Codon2AminoAcid.txt', 'r')
    data=fin.readlines()
    Codon2AA={}
    for code in data:
        codon,aa=code.strip('\n').split( ' ')
        Codon2AA[codon] = aa
        # print(codon)
        # print(aa)
    # print(Codon2AA)
    return(Codon2AA)

### My solution
def PROT(file_name,Codon2AA):
    fin = open(file_name, 'r')  # input file
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    rna = fin.readline().strip('\n')
    prot=''
    for index in range(0, len(rna), 3):
        codon = (rna[index: index + 3])
        if(Codon2AA[codon] != 'Stop' ):
            prot += Codon2AA[codon]
    # print(prot)
    fout.write(prot + "\n")

start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt'
input_fname= 'rosalind_prot.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file
AA_table= get_Codon_2_AminoAcid_table()
PROT(input_fname,AA_table) # 0.0021431446075439453 seconds

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




