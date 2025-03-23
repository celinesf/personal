
fin= open("input.txt", 'r') # input file
fout= open("output.txt", 'w') # output file

data = fin.readlines()
ntests= 0 # header: number of tests/lines in input file

count = 0

for line in data:
    count += 1
    if(count==1):
        ntests = int(line.strip())
    else:
        #elements= list((line.strip().split(' ')))
        elements= list(map(int, line.strip().split(' ')))
        print(elements)
        print(sum(elements))
        fout.write( str(sum(elements)) + "\n")