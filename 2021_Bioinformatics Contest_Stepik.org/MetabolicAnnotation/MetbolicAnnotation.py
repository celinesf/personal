import pandas as pd
import numpy as np
from itertools import product
from decimal import *
import time
from functools import wraps

t0 = time.time()
getcontext().prec = 6

# file = "0"
# modsi= 5
# modjk = 10
#Total time running 0.05088210105895996
file="1"
modsi= 5
modjk = 10
#Total time running 0.04101991653442383
#file="2"
# file="3"
# modsi= 999
# modjk = 999
# file="4"
# modsi= 99
# modjk = 999999
# file="5"
# modsi= 99
# modjk = 99999

fin = open(file + ".txt", 'r')  # input file
fout = open(file + "_output.txt", 'w')  # output file

ntests = int(fin.readline().strip())  # header/ number of tests t
print(ntests)
for count in range(0, ntests):
    print("test # " + str(count))
    #  M neutral metabolites with mass m_i.
    # K potential adducts with masses a_i
    # N measured signals * s_i > 0 *
    M, K, N = list(map(int, fin.readline().strip().split(' ')))
    print(str(M) + " " + str(K) + " " + str(N))
    mj = np.round([float(i) for i in (fin.readline().strip().split(' '))], 7)
    ak = np.round([float(i) for i in (fin.readline().strip().split(' '))], 7)
    si = np.round([float(i) for i in (fin.readline().strip().split(' '))], 7)
    # print(mj)
    # print(ak)
    # print(si)

    # combinatorics of mi and a1
    df1 = pd.DataFrame.from_records(list(i for i in product(mj, ak)), columns=['mj', 'ak'])
    df = pd.DataFrame.from_records(list(i for i in product(list(range(1, M + 1)), list(range(1, K + 1)))),
                                   columns=['j', 'k'])
    df['mj'] = df1['mj']
    df['ak'] = df1['ak']
    df1 = None
    # Sum of mj+ak
    df["mi+ak"] = np.round(df["mj"] + df["ak"], 7)
    # print(df)
    df = df[df["mi+ak"] > 0]
    # print(df)

    # loop on si
    for i in range(0, N):
        if i % modsi == 0:
            print("signal # " + str(i) + " out of "+ str (N))
        df['delta' + str(i)] = None
        # print('i '+ str(i) + " si " + str( si[i]))
        min = 10000000
        minjk = None
        for jk in df.index:
            if jk % modjk == 0 and jk > 0:
                print("jk # " + str(jk) + " out of " + str(len(df.index)))
            # print('jk '+ str(jk))
            # print(df.at[jk, 'mi+ak'])
            # print(si[i])

            df.at[jk, 'delta' + str(i)] = np.round(abs(df.at[jk, 'mi+ak'] - si[i]), 7)
            # print(df.at[jk, 'delta'+ str(i)])
            if df.at[jk, 'delta' + str(i)] < min:
                min = np.round(df.at[jk, 'delta' + str(i)], 7)
                minjk = jk

                # print("here minjk " + str(jk) + " min "+ str(min) )
            if min == 0:
                break

        # print(" ".join(str(x) for x in df.loc[ minjk, list(df.columns[0:2])]) + "\n")
        fout.write(" ".join(str(x) for x in df.loc[minjk, list(df.columns[0:2])]) + "\n")

t1 = time.time()
print ("Total time running " + str(t1-t0))