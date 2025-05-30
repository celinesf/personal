# Finding a Shared Motif

## Problem
A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. 
For example, "CG" is a common substring of "A**CG**TA**CG**T" and "AAC**CG**TATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "A**CGTA**CGT" and "AADC**CGTA**TA".

Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".

Given: A collection of $k$
 ($k≤100$
) DNA strings of length at most 1 kbp each in FASTA format.

Return: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution.)

## Sample Dataset
\>Rosalind_1<br>
GATTACA<br>
\>Rosalind_2<br>
TAGACCA<br>
\>Rosalind_3<br>
ATACA<br>

## Sample Output
AC