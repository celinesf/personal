import numpy as np
import json
import networkx as nx
import matplotlib.pyplot as plt
import time
t0 = time.time()
# file = "0"
# modnode= 9
# file="1" #94853
# modnode= 999
# file="2" # 94877
# modnode= 999
# file="3"
# modnode= 9999
file="4"
modnode= 9999
# file="5"
# modnode= 9999

fin = open("test" + file, 'r')  # input file
fout = open(file + "_output.txt", 'w')  # output file

npheno = int(fin.readline().strip())  # header/ number of tests t
print("Number phenotypes: " + str(npheno)  )
### phenotype vertices parents
P = fin.readline().strip().split(' ')  # .astype(np.int)
# print(P)
print("Done reading Phenotype " + str( time.time() - t0))
t1= time.time()
IC = fin.readline().strip().split(' ')  # .astype(np.int)
# print(IC)
print("Done reading IC - staring graph "+ str( time.time() - t1) )

t1= time.time()
### phenotype graph
GP = nx.DiGraph()
DIC = {1: {'IC': int(IC[0])}}
for v in range(0, npheno - 1):
    GP.add_edge(int(P[v]), v + 2)
    DIC[v + 2] = {'IC': int(IC[v + 1])}
#print(json.dumps(DIC, indent=4))
nx.set_node_attributes(GP, DIC)
# print( GP.edges())
# print(nx.is_directed_acyclic_graph(GP))
# print(json.dumps(nx.to_dict_of_dicts(GP), indent=4))

#print(nx.get_node_attributes(GP,'IC'))

DCI = None  # release memory
print("Done graph -  start LCA/ diseases reading "+ str( time.time() - t1) )

t1= time.time()
LCA ={}
### Disease set - calculate LCA to all node pairs considered - recoring IC for LCA
mdiseases = int(fin.readline().strip())  # num disease
print("Num diseases " + str (mdiseases)  )
DM=np.zeros((mdiseases, npheno), int) # list of diseasese phenotype vetices

### only look at those pairs of nodes
pairsnodes = []
pairsdic = {}
for i in range(0, mdiseases):
    if int(v) % modnode == 0 & int(v) != 0:
        print(" disease v " + str(v))
    MI = fin.readline().strip().split(' ')
    for v in MI[1: int(MI[0])+1]:
        if int(v) % modnode ==0 & int(v) !=0:
            print(" disease v " + str(v))
        DM[int(i)][int(v)-1] = v
        # print(i)
        # print(v)
        for key in GP:
            pair = ''
            if int(key) < int(v) :
                pair = str(key) + "-" +str(v)
            else:
                pair =str(v) +"-" +str(key)
            ##print(pair)
            if pair not in pairsdic:
                pairsdic[pair] = 1
                pairsnodes.append((int(key), int(v)))
                #LCA[pair]['lca'] = nx.lowest_common_ancestor(GP, int(key), int(v), default=None)
                #LCA[pair]['IC'] = GP.nodes[ LCA[pair]['lca'] ]['IC']
# print(pairsdic)
# print(pairsnodes)
print(" Before lca")
lca= nx.tree_all_pairs_lowest_common_ancestor(GP, root = None, pairs = pairsnodes )
# print(lca)

print(" Before LCA" )
LCA={}
for (u, v), l in lca:
    # print('u:' + str(u), 'v:' + str(v) ,'l:' + str(l) )
    pair = ''
    if int(u) < int(v):
        pair = str(u) + "-" + str(v)
    else:
        pair = str(v) + "-" + str(u)
    # print(pair)
    if pair not in LCA:
        LCA[pair] = {}
        LCA[pair]['lca'] = l
        LCA[pair]['IC'] = GP.nodes[ l ]['IC']

print(" after LCA")
GP.clear()
# print(DM)

#print(json.dumps(LCA, indent=4))

# print(nx.lowest_common_ancestor(GP,5,4))
#
# np.random.seed(8)
# plt.subplot(121)
# nx.draw(GP, with_labels = True)
# plt.show()

### Disease set
nqpatients = int(fin.readline().strip())  # num patients
print("Number phenotypes: " + str(npheno) + "\n" + "Number diseases: " + str(
    mdiseases) + "\n" + "Number Patients: " + str(nqpatients) + " - "+ str( time.time() - t1) )

# print(DM)

# loop on patients
for p in range(0, nqpatients):
    t1 = time.time()
    CQ = fin.readline().strip().split(' ')
    # print(CQ)
    ### loop on patient p phenotype set Qp
    pdiagnosis = 0 ## disease of choice
    pmaxsumic = 0 ## max sum of IC for patient p
    for mi in range(0, mdiseases):
        # print('mi ' + str (mi))
        # print(DM[mi])
        pmsumic = 0 ## sum of max IC for patient p and disease mi
        for q in CQ[1:len(CQ)]:
            # print(q)
            maxdic = 0 ## max IC for patient p and disease pgenotype d and patient phenotype q
            for d in DM[mi][np.where(DM[mi] > 1)]:
                # print(d)
                # get LCA for (q,d) pair
                pair = ''
                if int(d) < int(q):
                    pair = str(d) + "-" + str(q)
                else:
                    pair = str(q) + "-" + str(d)

                # print(pair)
                # print(LCA[pair])
                if LCA[pair]['IC'] > maxdic:
                    maxdic = LCA[pair]['IC']
            # print("maxdic " + str(maxdic))
            pmsumic = pmsumic + maxdic
        # print('pmsumic' + str(pmsumic))
        if pmsumic > pmaxsumic:
            pmaxsumic = pmsumic
            pdiagnosis = int(mi) + 1
    print("patient " +str(p) + " diagnosis " + str(pdiagnosis)  + " pmaxsumic " + str(pmaxsumic) + " - "+ str( time.time() - t1) )
    fout.write(str(pdiagnosis) + "\n")