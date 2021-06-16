import pandas as pd
import numpy as np

#file="0"
#file="1"
file="2"

fin= open(file +".txt", 'r') # input file
fout= open(file+"_output.txt", 'w') # output file

ntests= int(fin.readline().strip()) # header/ number of tests t

for count in range(0, ntests):
    nseq,length = list(map(int, fin.readline().strip().split(' '))) # read nseq and lenght
    df= pd.DataFrame({})
    # populate a dataframe with sequences as columns
    for s in range(0,nseq):
        seq= pd.Series(list(fin.readline().strip()),dtype="string")
        df[s] = seq
    # create a col of each state per position
    df['state'] = df.agg(''.join, axis=1)
    df= df.drop(df.columns[0:nseq],1)
    df['pos'] = df.index +1
    # group by states
    df['states'] = df['state'].map( df.groupby(['state'])['pos'].apply(list).to_dict())
    # get 1st position for each state
    df.loc[:, 's'] = df.states.map(lambda x: x[0])

    states = df['s'].unique()
    print(states)
    fout.write(str(len(states))+"\n")

    # fix state number based on position
    for nums in (range(len(states))):

        if states[nums] != nums+1:
            df.loc[df['s'] == states[nums], 's'] = nums+1

    fout.write(" ".join(str(x) for x in df['s']) + "\n")
