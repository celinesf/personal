import time
import sys

### My solution
def HAMM(file_name):
    fin = open(file_name, 'r')  # input file
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    seq1 = fin.readline().strip('\n')
    seq2 = fin.readline().strip('\n')
    #print(seq1 +" "+ seq2)
    count = sum(1 for a, b in zip(seq1, seq2) if a != b)
    print(count)
    fout.write(str(count) + "\n")

start_time = time.time() # report time of function

#input_fname= 'Sample_Dataset.txt'
input_fname= 'rosalind_hamm.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file
HAMM(input_fname) # 0.0021431446075439453 seconds

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




