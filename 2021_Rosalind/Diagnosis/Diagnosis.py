import pandas as pd
import numpy as np
from itertools import product
from decimal import *

getcontext().prec = 6

#file = "0"
#file="1"
#file="2"
file="3"
# file="4"
# file="5"

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
    # print(str(M) + " " + str(K) + " " + str(N))
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
        df['delta' + str(i)] = None
        # print('i '+ str(i) + " si " + str( si[i]))
        min = 1000000
        minjk = None
        for jk in df.index:
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
