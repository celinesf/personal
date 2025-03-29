# Mortal Fibonacci Rabbits

## Problem
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, which followed the recurrence relation $Fn=Fn−1+Fn−2$
 and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a dynamic programming solution in the case that all rabbits die out after a fixed number of months. 
See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).

![Figure 4](mortal_rabbit_tree.png)
A figure illustrating the propagation of Fibonacci's rabbits if they die after three months.
Figure 4 of 4

**Given**: Positive integers $n≤100$ and $m≤20$.

**Return**: The total number of pairs of rabbits that will remain after the $n^{th}$ month if all rabbits live for $m$
 months.

## Sample Dataset
6 3
## Sample Output
4

age = [0,0,0] # age at 0,2,3 months
generation 1 [1,0,0]
gen2 [0,1,0] # only 2 months old rabbit reproduce
gen3 [1,0,1] # age[1]+age[2]=age[0] new rabbit pair
gen4 [1,1,0]
gen5 [1,1,1]
gen6 [2,1,1] = 4 pairs