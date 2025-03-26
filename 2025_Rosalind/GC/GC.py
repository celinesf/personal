
import time
import sys

def my_function():
    time.sleep(1)  # Simulate some work
    return

### Rosalind Soulution
def GC_Rosalind(file_name):
    dnafile = open(file_name, 'r')
    raw = dnafile.read()
    d = {}
    for seqblock in raw.split(">")[1:]:
        parts = seqblock.split("\n")
        fasta = parts[0]
        seq = ''.join(parts[1:])
        gc = 100 * (seq.count("G") + seq.count("C")) / float(len(seq))
        d[gc] = fasta
    print( d[max(d)], max(d))

### My Solution
def GC_percentage(seq):
    return round((seq.count('G')+seq.count('C'))*100/len(seq),6)

def GC(file_name):
    fin = open(file_name, 'r')  # input file

    fout_name=file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')  # output file

    seq_name=''
    seq=''
    max_seq_name =''
    max_GC=0
    for line in fin:
        #print(line)
        if '>' in line:
            if seq != '' :
                if(GC_percentage(seq) > max_GC):
                    max_seq_name =seq_name
                    max_GC = GC_percentage(seq)

            seq_name = line.strip(">\n")
            seq = ''

        else:
            seq += line.strip("\n")
            #print(seq)
        #GC_count = round((seq.count('G')+seq.count('C'))*100/len(seq),6)
    if (GC_percentage(seq) > max_GC):
        max_seq_name = seq_name
        max_GC = GC_percentage(seq)

    print(max_seq_name + "\n" + str(max_GC) + "\n")
    fout.write(max_seq_name + "\n" + str(max_GC) + "\n")


start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt'
# input_fname= 'rosalind_gc.txt'
input_fname= 'rosalind_gc (4).txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file
#GC(input_fname) #0.0005459785461425781 seconds
GC_Rosalind(input_fname) #0.0003268718719482422 seconds

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




