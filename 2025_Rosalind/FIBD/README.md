# Consensus and Profile

## Problem
A **matrix** is a rectangular table of values divided into rows and columns. An $m×n$
 matrix has $m$ rows and $n$ columns. 
Given a matrix $A$
, we write $A_{i,j}$
 to indicate the value found at the intersection of row $i$ and column $j$.

Say that we have a collection of DNA strings, all having the same length $n$
. Their profile matrix is a $4×n$
 matrix $P$
 in which $P_{1,j}$
 represents the number of times that 'A' occurs in the 
$j^{th}$ position of one of the strings, 
$P_{2,j}$  represents the number of times that C occurs in the $j^{th}$ position, and so on (see below).

A consensus string $c$
 is a string of length $n$
 formed from our collection by taking the 
most common symbol at each position; 
the $j^{th}$ symbol of $c$
 therefore corresponds to the symbol having the maximum value in the $j^{th}$ column of the profile matrix. Of course, 
there may be more than one most common symbol, leading to multiple possible consensus strings.

DNA Strings<br>
A T C C A G C T<br>
G G G C A A C T<br>
A T G G A T C T<br>
A A G C A A C C<br>
T T G G A A C T<br>
A T G C C A T T<br>
A T G G C A C T<br>

Profile<br>
A   5 1 0 0 5 5 0 0<br>
C   0 0 1 4 2 0 6 1<br>
G   1 1 6 3 0 1 0 0<br>
T   1 5 0 0 0 1 1 6<br>

Consensus<br>
A T G C A A C T

**Given**: A collection of at most 10 DNA strings of equal length 
(at most 1 kbp) in FASTA format.

**Return**: A consensus string and profile matrix for the collection. 
(If several possible consensus strings exist, then you may return any one of them.)

## Sample Dataset
\>Rosalind_1<br>
ATCCAGCT<br>
\>Rosalind_2<br>
GGGCAACT<br>
\>Rosalind_3<br>
ATGGATCT<br>
\>Rosalind_4<br>
AAGCAACC<br>
\>Rosalind_5<br>
TTGGAACT<br>
\>Rosalind_6<br>
ATGCCATT<br>
\>Rosalind_7<br>
ATGGCACT<br>

## Sample Output
ATGCAACT<br>
A: 5 1 0 0 5 5 0 0<br>
C: 0 0 1 4 2 0 6 1<br>
G: 1 1 6 3 0 1 0 0<br>
T: 1 5 0 0 0 1 1 6<br>