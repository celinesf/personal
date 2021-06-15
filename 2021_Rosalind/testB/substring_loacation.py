
def find_index(string, motif):
    enumerated = [index+1  for index, value in enumerate(string) if string[index:index+(len(motif))] == motif]
    return enumerated

fin= open("input.txt", 'r') # input file
fout= open("output.txt", 'w') # output file

ntests = int(fin.readline().strip()) # header/ number of sequences/motif pair

for count in range(0,ntests):
    seq= fin.readline().strip()
    motif= fin.readline().strip()
    indexes= find_index(seq, motif)
    fout.write(" ".join(str(x) for x in indexes) + "\n")
