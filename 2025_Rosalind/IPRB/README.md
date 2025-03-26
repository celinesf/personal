# Mendel's First Law

## Problem
Probability is the mathematical study of randomly occurring phenomena. We will model such a phenomenon with a random variable, which is simply a variable that can take a number of different distinct outcomes depending on the result of an underlying random process.
![img.Figure 2](balls_tree.png)
The probability of any outcome (leaf) in a probability tree diagram is given by the product of probabilities from the start of the tree to the outcome. For example, the probability that X is blue and Y is blue is equal to (2/5)(1/4), or 1/10.
Figure 2 of 2

For example, say that we have a bag containing 3 red balls and 2 blue balls. If we let $X$
 represent the random variable corresponding to the color of a drawn ball, then the probability of each of the two outcomes is given by $Pr(X=red)=35$
 and $Pr(X=blue)=25$
.

Random variables can be combined to yield new random variables. Returning to the ball example, let $Y$
 model the color of a second ball drawn from the bag (without replacing the first ball). The probability of $Y$
 being red depends on whether the first ball was red or blue. To represent all outcomes of $X$
 and $Y$
, we therefore use a probability tree diagram. This branching diagram represents all possible individual probabilities for $X$
 and $Y$
, with outcomes at the endpoints ("leaves") of the tree. The probability of any outcome is given by the product of probabilities along the path from the beginning of the tree; see [Figure 2](balls_tree.png) for an illustrative example.

An event is simply a collection of outcomes. Because outcomes are distinct, the probability of an event can be 
written as the sum of the probabilities of its constituent outcomes. For our colored ball example, let $A$
 be the event "$Y$
 is blue." $Pr(A)$
 is equal to the sum of the probabilities of two different outcomes: $Pr(X=blue$ and $Y=blue)+Pr(X=red$ and $Y=blue)$
, or 310+110=25
 (see [Figure 2](balls_tree.png) above).

**Given**: Three positive integers $k$, $m$, and $n$
, representing a population containing $k+m+n$ organisms: $k$
 individuals are homozygous dominant for a factor, $m$
 are heterozygous, and $n$
 are homozygous recessive.

**Return***: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype). Assume that any two organisms can mate.

## Sample Dataset
2 2 2
## Sample Output
0.78333

### Understanding the problem: 
We need to find the probability that two randomly selected organisms, regardless of their genotype ($k$, $m$, or $n$), 
will produce offspring with at least one dominant allele, thus displaying the dominant phenotype.

Total population: We have 
- k=2 homozygous dominant,
- m=2 heterozygous, and 
- n=2 homozygous recessive organisms, 
- totaling 6 organisms.

Total possible pairings: 
The total number of possible pairings is calculated using combinations: (6 choose 2) = 15.
Valid pairings (those that produce a dominant offspring):
- k&k (homozygous dominant x homozygous dominant): All offspring will be dominant. There is 1 pairing (2 choose 2) = 1.
- k&m (homozygous dominant x heterozygous): All offspring will be dominant. There are 2 pairings (2 * 2) = 4.
- k&n (homozygous dominant x homozygous recessive): All offspring will be dominant. There are 2 pairings (2 * 2) = 4.
- m&m (heterozygous x heterozygous): 3/4 of the offspring will be dominant. There is 1 pairing (2 choose 2) = 1, and 3/4 * 1 = 0.75.
- m&n (heterozygous x homozygous recessive): 1/2 of the offspring will be dominant. There are 2 pairings (2 * 2) = 4, and 1/2 * 4 = 2.
- n&n (homozygous recessive x homozygous recessive): All offspring will be recessive. There is 1 pairing (2 choose 2) = 1.
Total valid pairings: 1 + 4 + 4 + 0.75 + 2 = 11.75.
Probability of a dominant offspring: (Total valid pairings) / (Total possible pairings) = 11.75 / 15 = 0.7833.


### Hint
Consider simulating inheritance on a number of small test cases in order to check your solution.

## How to Calculate the Probability:
Define the Population:
Assume a population with the following genotypes:
- $k$ individuals are homozygous dominant (AA)
- $m$ individuals are heterozygous (Aa)
- $n$ individuals are homozygous recessive (aa) 

### Consider all possible mating combinations:
- AA x AA (homozygous dominant x homozygous dominant)
- AA x Aa (homozygous dominant x heterozygous)
- AA x aa (homozygous dominant x homozygous recessive)
- Aa x Aa (heterozygous x heterozygous)
- Aa x aa (heterozygous x homozygous recessive)
- aa x aa (homozygous recessive x homozygous recessive)

### Determine the probability of each mating combination:
The probability of two individuals mating is the product of their individual probabilities of being selected. 
For example, the probability of an AA individual mating with another AA individual is $(k/(k+m+n)) * (k/(k+m+n))$. 
Similarly, the probability of an AA individual mating with an Aa individual is $(k/(k+m+n)) * (m/(k+m+n))$. 
Calculate the probability of offspring with a dominant phenotype for each mating combination:
- AA x AA: All offspring will have the dominant phenotype (100% probability). 
- AA x Aa: All offspring will have the dominant phenotype (100% probability). 
- AA x aa: All offspring will have the dominant phenotype (100% probability). 
- Aa x Aa: 75% of offspring will have the dominant phenotype (AA or Aa). 
- Aa x aa: 50% of offspring will have the dominant phenotype (Aa). 
- aa x aa: All offspring will have the recessive phenotype (0% probability of dominant phenotype). 
