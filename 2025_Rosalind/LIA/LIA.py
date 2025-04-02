import time
import sys
import itertools
import math

### My solution
# k individuals are homozygous dominant for a factor,
# m  are heterozygous, and
# n  are homozygous recessive.
def IPRB(file_name):
    fin = open(file_name, 'r')  # input file
    fout_name = file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')
    k, m, n = fin.readline().strip('\n').split(' ')
    k = int(k)
    m = int(m)
    n = int(n)
    print(str(k) +" "+ str(m)+" "+ str(n))

    Pair_A={'AAAA' : 1*math.comb(k, 2),'AAAa' :  1*k*m, 'AAaa' : 1 * k*n,
              'AaAa': .75 *math.comb(m, 2), 'Aaaa': .5 * m*n, 'aaaa': 0* math.comb(n, 2)}
    Pair_a = {'AAAA': 0 * math.comb(k, 2), 'AAAa': 0 * k * m, 'AAaa': 0 * k * n,
              'AaAa': .25 * math.comb(m, 2), 'Aaaa': .5 * m * n, 'aaaa': 1 * math.comb(n, 2)}

    p= (Pair_A['AAAA']+Pair_A['AAAa']+Pair_A['AAaa']+Pair_A['AaAa']+Pair_A['Aaaa'])/math.comb(k+m+n, 2)
    total= (Pair_A['AAAA']+Pair_A['AAAa']+Pair_A['AAaa']+Pair_A['AaAa']+Pair_A['Aaaa']+
            Pair_a['AaAa'] +Pair_a['Aaaa'] +Pair_a['aaaa'] )/math.comb(k+m+n, 2)#should equal1
    print(round(p,6))
    print(total)
    print(round( 1-(Pair_a['AaAa'] +Pair_a['Aaaa'] +Pair_a['aaaa'] )/math.comb(k+m+n, 2) ,6))
    # fout.write(str(count) + "\n")

start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt''
# input_fname= 'rosalind_lia.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file

LIA(input_fname) # 0.0021431446075439453 seconds

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




