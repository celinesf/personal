import time
import sys
import itertools
import math

### My solution
def IEV(file_name):
    fin = open(file_name, 'r')  # input file
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    pairs = fin.readline().strip('\n').split(' ')
    c = list(map(int, pairs))
    print((c))

    Pair_A={'AAAA' : 2*int(c[0]),'AAAa' :  2*int(c[1]), 'AAaa' : 2 * int(c[2]),
              'AaAa': .75 *2*int(c[3]), 'Aaaa': .5 *2* int(c[4]), 'aaaa': 0* 2* int(c[5])}
    Pair_a = {'AAAA': 0 *2* int(c[0]), 'AAAa': 0*2 * int(c[1]), 'AAaa': 0 *2* int(c[2]),
              'AaAa': .25 *2* int(c[3]), 'Aaaa': .5 *2*int(c[4]), 'aaaa': 1*2 * int(c[5])}

    p= (Pair_A['AAAA']+Pair_A['AAAa']+Pair_A['AAaa']+Pair_A['AaAa']+Pair_A['Aaaa'])
    total= (Pair_A['AAAA']+Pair_A['AAAa']+Pair_A['AAaa']+Pair_A['AaAa']+Pair_A['Aaaa']+
            Pair_a['AaAa'] +Pair_a['Aaaa'] +Pair_a['aaaa'] )/ (sum(c)*2)
    print(int(p))
    print(total)
    # print(round( 1-(Pair_a['AaAa'] +Pair_a['Aaaa'] +Pair_a['aaaa'] )/math.comb(k+m+n, 2) ,6))
    fout.write(str(int(p)) + "\n")

start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt'
input_fname= 'rosalind_iev.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file

IEV(input_fname) # 0.0021431446075439453 seconds

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




