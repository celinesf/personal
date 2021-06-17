#Metabolite Annotation
In this problem, you need to submit the answers to the tests in a plain text file.

The problem is exact. Your score is calculated using the percentage of the correct answers. For example, if there are *T* test subcases and your answer is correct in KK of them, then your score for the test will be *P.K/T* , where *P* is the maximal score.

The maximal total score for the problem is 10001000 points.

You can download all the tests here: https://stepik.org/media/attachments/lesson/541850/all.zip.

Author: Alexey Sergushichev
Tests: Grigorii Shovkoplias

Mass spectrometry is a technique that can be used to detect the presence of metabolites (biochemical compounds) in a sample. In this technique, a neutral metabolite is ionized by gaining or losing a charged fragment (adduct), and then the mass-to-charge ratio is measured for this ionized metabolite. Your task is to annotate mass-spectrometry results: find for a measured mass-to-charge ratio from which metabolite it could come from.

Formally, there is a database of *M* neutral metabolites with masses *m_i > 0* and a database of *K* potential adducts with masses *a_i* ( *a_i* can be both positive and negative). Then there are *N* measured signals *s_i > 0*. Each signal *s_i* corresponds to some metabolite *m_j* (*1 ≤ j ≤ M*) with an adduct *a_k* (*1 ≤ k ≤ K*) and some noise Δ (can be both positive and negative), that is  *s_i = m_j + a_k + Δ*, with *m_j + a_k > 0*.

Your task is to find for each of *N* signals *s_i* the pair of metabolite *m_j* and adduct *a_k* with the closest sum *m_j + a_km*.

###Input format
The first line of the input contains one integer *T (*1 ≤ T ≤ 3*) − the number of test cases.

Each test case is specified by four lines. The first line of each test case contains three integer numbers *M*, *K*, *N*. The second line contains *M* numbers *m_i* − masses of metabolites (*0 < m_i ≤1000*). The third line contains *K* numbers *a_i* − masses of adducts (*−1000 ≤ a_i ≤1000). The fourth line contains *N* numbers *s_i* − masses of signals (*0 <s_i ≤1000*). **All the masses are indicated with exactly six decimal places**.

###Output format
For each signal *s_i* of each test case, print numbers *j* and *k* such that  *s_i = m_j + a_k + Δ, m_j + a_k > 0* and an absolute value of *Δ* is smallest possible. If there are multiple numbers *j* and *k* with same absolute value of *Δ* for some signal, you can print any of them.

###Sample input
 
3  
2 2 5  
1.000002 0.000002  
0.500000 -0.500000  
0.500001 0.500002 0.500003 1.000000 0.000001  
2 2 5  
1.000002 0.000001  
0.500000 -0.500000  
0.500001 0.500002 0.500003 1.000000 0.000001  
5 4 7  
0.000001 0.000002 0.000003 0.000004 0.000005  
0.000002 0.000010 0.000001 -0.000001  
0.000001 0.000002 0.000100 0.000005 0.000020 0.000010 0.000003  

###Sample output
 
1 2  
1 2  
1 2  
1 2  
1 2  
2 1  
1 2  
1 2  
1 2  
2 1  
2 4   
1 3    
5 2  
3 1  
5 2  
1 2  
1 1  