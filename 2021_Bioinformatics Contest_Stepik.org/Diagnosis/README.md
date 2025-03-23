# Diagnosis

In this problem, you are given seven different tests for the proposed problem in seven different tabs. You should submit
only answers in plain text format.

The points are distributed as follows: *100* points for the first test, *250* points for tests 2-5, and *500* points for
tests 6 and 7.

All tests are artificially generated.

Author: German Demidov Tests: Aleksandra Drozdova

## Diagnosis

Determining the correct diagnosis is a crucial step in patient treatment. However, it can be far from trivial, as the
disease can have multiple phenotypic manifestations with different degrees of specificity, partly overlapping with other
conditions. To formalize the diagnostics process, Human Phenotype Ontology (HPO) was invented, using which both the
disease manifestations and patient phenotypic traits can be described. In this problem, your task will be to identify
the patients' diseases given their clinical phenotypes.

In this problem, you are given a human phenotype tree. Each vertex of this tree corresponds to a phenotypic abnormality.
The abnormalities on the lower levels are more specific, with the specificity of the vertex *v* defined as information
content *IC(v)* value. Additionally, you are given descriptions of several diseases, each defined as a set of
abnormalities, that is vertices of the phenotype tree. Finally, you are given descriptions of patients, similarly
defined as sets of phenotype tree vertices, describing their clinical phenotypes.

Your task is to find for each of the patients their most likely disease. More precisely, for every patient *p* with the
phenotype set *Q_p* find the disease *m* with the phenotype set *D_m* which maximizes the value of  
*∑(    max IC(LCA(q,d)))*,  
q∈Q_p d∈D_m  

where *LCA(q,d)* is the lowest common ancestor of phenotype vertex *q* and
phenotype vertex *d* and *IC(v)* is the information content of vertex *v*. If there are several diseases with the same
maximal value, then any of them can be returned.

###Input Format 
The first line of the input file contains one integer *n* — the number of vertices in the phenotype tree.
The vertices are identified by numbers *1, 2, ..., n* with vertex *1* being the root of the tree.

The second line contains *n - 1* integers — the parent identifiers for the vertices *2, 3, ..., n*.

The third line contains *n* integers — information content values of the corresponding vertices *1, 2, ..., n*.
It is guaranteed that for every vertex *v* its information content *IC(v)* is greater than the information content of 
its parent vertex.

The fourth line contains one integer *m* — the number of diseases.

Next *m* lines contain descriptions of diseases. The *i*-th line contains an integer *cm_i* — the number of 
vertices in the phenotypes tree describing the *i*-th disease, followed by *cm_i* different integers — the 
identifiers of vertices describing the *i*-th disease.

The next line contains one integer *nq* — the number of patients. 

Next *nq* lines contain descriptions of patients. The *i*-th line contains an integer *cq_i* — the number of
vertices in the phenotypes tree describing the *i*-th patient, followed by *cq_i* different integers — the
identifiers of vertices describing the *i*-th patient.

### Output Format

The output file should contain *q* lines. The *i*-th line should contain one integer — the identifier of a disease of
the *i*-th patient. The indexing starts from 1.

### Scoring

For all tests, you will receive points proportional to the number of diseases you recovered correctly.

### Sample Input

10  
1 1 3 3 4 4 5 5 5  
5 7 8 13 18 14 15 21 20 29  
2  
2 4 2  
1 10  
4  
3 5 9 8  
1 6  
2 7 10  
1 10

### Sample Output

2  
1  
2  
2  
