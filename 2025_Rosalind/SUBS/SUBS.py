import time
import sys
import re

### My solution
def SUBS(file_name):
    fin = open(file_name, 'r')  # input file
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    dna = fin.readline().strip('\n')
    t = fin.readline().strip('\n')
    # print(dna)
    # print(t)

    matches = [(match.start()+1) for match in re.finditer(f'(?=({t}))', dna)]

    print( ' '.join(map(str, matches)))
    fout.write(' '.join(map(str, matches))+ "\n")


start_time = time.time() # report time of function

input_fname= 'Sample_Dataset2.txt'
# input_fname= 'rosalind_subs.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file

SUBS(input_fname,) #

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




