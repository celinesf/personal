
import time
import sys

### Cheated on solution here and googled it. took me a minute to rap my head around the concept
def fib(n, k):#months, age
    print('n=' + str(n), " k=", str(k))
    ages = [1] + [0] * (k - 1) # array of size k # rabit of each month [0,1,2] at generation 1
    #print(ages)
    for i in range(n - 1): #generation 2, i=0 to generation n, i=4
        ages = [sum(ages[1:])] + ages[:-1] # age[0] = sum of rabit age 1 month ro death, age[1:k-1] = previous [0:k-2]
        #print(str(i)+ " "+ str(ages[1:]) + " "+ str(ages[:-1]) )
        print(ages)
    return sum(ages)


def FIBD(file_name):
    fin = open(file_name, 'r')  # input file

    fout_name=file_name.removesuffix('.txt') + "_output.txt"
    print(fout_name)
    fout = open(fout_name, 'w')  # output file
    data = fin.readlines()
    #print(str(data.split()))
    for line in data:
        n,k= line.split()
        Fn=fib(int(n),int( k))
        print(Fn)
        fout.write(str(Fn)+"\n")


start_time = time.time() # report time of function

input_fname= 'Sample_Dataset.txt'
input_fname= 'rosalind_fibd.txt'

#sys.stdout = open(input_fname.strip('.txt')+".out", mode='w') # redirect stdout to file
FIBD(input_fname) # 0.0021431446075439453 seconds

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution time: {execution_time} seconds")




