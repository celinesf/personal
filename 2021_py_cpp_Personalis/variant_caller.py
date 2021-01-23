#!/bin/env python
"""
Implement a naive variant caller.

Given a reference sequence and a read pileup, make one of the following variant calls at each position in the reference:
1. "nocall": not enough evidence to make a call
2. "ref": only the reference allele is present
3. "[ACGT]hom: an Alt allele (A,C,G, or T) is present in >75% of reads
4. "[ACGT]het: an Alt allele (A,C,G, or T) is present in 25%--75% of reads
5. "[ACGT]low: an Alt allele (A,C,G, or T) is present in <25% of reads

If more than 1 Alt allele is present, the variant call should contain all of them, comma delimited

For this exercise we are only considering SNVs

The call_variants function takes a min_depth argument; if the total read depth is less than this,
the position is "nocall"

Pileup notation:
"." or "," : Aligned base matches the ref allele at this position (+ or - strand)
[ACGT] or [acgt] : Aligned base is a mismatch to the reference (+ or - strand)
For this exercise, we are not considering strand information, so treat uppercase and lowercase equally.

IMPORTANT NOTE: the test_pileup data structure is a *pileup*: each string in the list is the set of
aligned bases at one position on the reference.  For example, the first string is only '.' and ',',
so at the first reference position (which is base 'T'), the pileup of aligned bases contains
only T base calls.
"""

"""
@author: Celine Becquet
@creation_date: 01/17/2021
"""

__author__ = "Celine Becquet"
__email__ = "celine.becquet@gmail.com"
__status__ = "dev"
__version__ = 1.0

import doctest
import numpy as np

test_ref = 'TTTAGAGCGC'
test_pileup = ['....,...,,,,.,...,',
               ',..,,,.c..C.,,,',
               'AAaAAaaaA.,AaaaAa',
               '.C.,.g,gG.G...,G.',
               '..,C,..C,,CCCC.C,c',
               '.,..tT,,.t',
               '',
               '.Taat..,.,...,',
               '.cc',
               'aA.,,A.,...,a.,Aa...,']


def call_variants(pileup, ref, min_depth):
    """ Make a variant call at each position on ref, using the pileup
    Args:
        pileup (list): list of pileup string for each ref position
        ref (list): list of reference allele [A,C,T,G] at each position
        min_depth (int): Minimum number of read (lenght of pileup) to make a call
    Example: Test call_variants
        >>> test_ref = 'TTTAGAGCGC'
        >>> test_pileup = ['....,...,,,,.,...,',',..,,,.c..C.,,,','AAaAAaaaA.,AaaaAa','.C.,.g,gG.G...,G.','..,C,..C,,CCCC.C,c','.,..tT,,.t','', '.Taat..,.,...,', '.cc','aA.,,A.,...,a.,Aa...,']
        >>> call_set = call_variants(test_pileup, test_ref, 5)
        >>> print(*call_set)
        ref [C]low [A]hom [C]low[G]het [C]het [T]het nocall [A]low[T]low nocall [A]het
    """
    # >
    call_set = []

    for i in range(0, len(ref)):
        # homogenize pileup, i.e., ignore strand information
        p = pileup[i].upper().replace(",", ref[i]).replace(".", ref[i])

        # 1. "nocall": not enough evidence to make a call
        if (len(p) < min_depth):
            call_set.append("nocall")
        else:
            counts = list(zip(*np.unique(np.array(list(p)), return_counts=True)))
            alt = ''
            for allele in counts:
                # 2. "ref": only the reference allele is present
                if allele[0] == ref[i] and allele[1] == len(p):
                    call_set.append("ref")
                # 3. "[ACGT]hom: an Alt allele (A,C,G, or T) is present in >75% of reads
                elif allele[0] != ref[i] and allele[1] / len(p) > 0.75:
                    alt = alt + "[" + allele[0] + "]hom"
                # 4. "[ACGT]het: an Alt allele (A,C,G, or T) is present in 25%--75% of reads
                elif allele[0] != ref[i] and allele[1] / len(p) <= 0.75 and allele[1] / len(p) >= 0.25:
                    alt = alt + "[" + allele[0] + "]het"
                # 5. "[ACGT]low: an Alt allele (A,C,G, or T) is present in <25% of reads
                elif allele[0] != ref[i] and allele[1] / len(p) < 0.25:
                    alt = alt + "[" + allele[0] + "]low"
            if alt != '':
                call_set.append(alt)

    return call_set


def main():
    # Keep this line:
    call_set = call_variants(test_pileup, test_ref, 5)

    # TODO: print the ref positions (1-10), and the variant call at each position
    # (Bonus: can you use list comprehension to print the results with one line of code?)
    print("\n".join(("pos" + str(i + 1) + ": " + v for i, v in enumerate(call_set))))


if __name__ == "__main__":
    main()
    #doctest.testmod(verbose=True)
