# Mendel's First Law

![img.Figure 2](dihybrid_cross.png)

Figure 1. Mendel's second law dictates that every one of the 16 possible assignments of parental alleles is equally likely. The Punnett square for two factors therefore places each of these assignments in a cell of a 4 X 4 table. The probability of an offspring's genome is equal to the number of times it appears in the table, divided by 16.

## Problem
![img.Figure 2](two_dice.png)

Figure 2. The probability of each outcome for the sum of the values on two rolled dice (black and white), broken down depending on the number of pips showing on each die. You can verify that 18 of the 36 equally probable possibilities result in an odd sum.

Two events $A$
 and B
 are independent if $Pr(A$ and $B)$
 is equal to $Pr(A)×Pr(B)$
. In other words, the events do not influence each other, so that we may simply calculate each of the individual probabilities separately and then multiply.

More generally, random variables $X$
 and $Y$
 are independent if whenever $A$
 and $B$
 are respective events for $X$
 and $Y$
, $A$
 and $B$
 are independent (i.e., $Pr(A$ and $B)=Pr(A)×Pr(B)$
).

As an example of how helpful independence can be for calculating probabilities, let $X$
 and $Y$
 represent the numbers showing on two six-sided dice. Intuitively, the number of pips showing on one die should not 
affect the number showing on the other die. If we want to find the probability that $X+Y$
 is odd, then we don't need to draw a tree diagram and consider all possibilities. We simply first note that for $X+Y$
 to be odd, either $X$
 is even and $Y$
 is odd or $X$
 is odd and $Y$
 is even. In terms of probability, $Pr(X+Y$ is odd$)=Pr(X$ is even and $Y$ is odd$)+Pr(X$ is odd and $Y$ is even$)$
. Using independence, this becomes $[Pr(X$ is even)×$Pr(Y$ is odd)$]+[Pr(X$ is odd$)×Pr(Y$ is even$)]$
, or $(1\2)^2+(1/2)^2=12$
. You can verify this result in Figure 2, which shows all 36 outcomes for rolling two dice.

**Given**: Two positive integers $k$
 ($k≤7$
) and $N$
 ($N≤2^k$
). In this problem, we begin with Tom, who in the 0th generation has genotype Aa Bb. Tom has two children in the 1st generation, each of whom has two children, and so on. Each organism always mates with an organism having genotype Aa Bb.

**Return**: The probability that at least $N$
 Aa Bb organisms will belong to the $k^{
th} generation of Tom's family tree (don't count the Aa Bb mates at each level). Assume that Mendel's second law 
holds for the factors.

## Sample Dataset
2 1
## Sample Output
0.684