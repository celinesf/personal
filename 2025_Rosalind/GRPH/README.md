# Overlap Graphs

## Problem
A graph whose nodes have all been labeled can be represented by an adjacency list, in which each row of the list contains the two node labels corresponding to a unique edge.

A **directed graph* (or digraph) is a graph containing directed edges, each of which has an orientation. That is, a directed edge is represented by an arrow instead of a line segment; the starting and ending nodes of an edge form its tail and head, respectively.
The directed edge with tail $v$
 and head $w$
 is represented by $(v,w)$
 (but not by (w,v)
). A directed loop is a directed edge of the form $(v,v)$
.

For a collection of strings and a positive integer $k$
, the overlap graph for the strings is a directed graph $O_k$
 in which each string is represented by a node, and string $s$
 is connected to string $t$
 with a directed edge when there is a length $k$
 suffix of $s$
 that matches a length $k$
 prefix of $t$
, as long as $s≠t$
; we demand $s≠t$
 to prevent directed loops in the overlap graph (although directed cycles may be present).

**Given**: A collection of DNA strings in FASTA format having total length at most 10 kbp.

**Return**: The adjacency list corresponding to $O_3$
. You may return edges in any order.

## Sample Dataset
\>Rosalind_0498<br>
AAATAAA<br>
\>Rosalind_2391<br>
AAATTTT<br>
\>Rosalind_2323<br>
TTTTCCC<br>
\>Rosalind_0442<br>
AAATCCC<br>
\>Rosalind_5013<br>
GGGTGGG<br>
## Sample Output
Rosalind_0498 Rosalind_2391<br>
Rosalind_0498 Rosalind_0442<br>
Rosalind_2391 Rosalind_2323<br>
