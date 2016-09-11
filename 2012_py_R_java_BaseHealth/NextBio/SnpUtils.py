#!/usr/bin/env python

"""
     functions to deal with alleles/snps
     06/19/13- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging
from NextBioUtils import NextBioUtils

class SnpUtils():
    def __init__(self):
        self.util = NextBioUtils()

    ''' compareAlleles
    '''
    def compareAlleles(self, alleles1, alleles2):
        logging.debug(' Function: compareAlleles - alleles1: %s , alleles2: %s' % (alleles1, alleles2))  
        found = 0
        for a1 in alleles1:
            if a1 in alleles2:
                found+=1
        if found == len(alleles1):
            return True
        else:
            return False

    ''' reverse
    '''
    def reverseAlleles(self, alleles):
        logging.debug(' Function: reverseAlleles - %s' % alleles)  
        reversed_allele = ""
        if alleles is not None:
            reversed_allele = []
            for index in range(0,len(alleles)):
                allele= alleles[index ]
                if allele == 'G':
                    reversed_allele.append('C')
                elif allele == 'C':
                    reversed_allele.append('G')
                elif allele == 'A':
                    reversed_allele.append('T')
                elif allele == 'T':
                    reversed_allele.append('A')
                elif allele == '-':
                    reversed_allele.append('-')
                else :
                    self.util.warnMe('warning', ' CANT REVERSE ALLELE %s (reverseAlleles) - alleles:%s' % (allele,alleles))  
            if len(reversed_allele) == 1:
                reversed_allele = reversed_allele[0]
        return reversed_allele


    ''' isSameAllele '''
    def isSameAllele(self, allele1, allele2):
        logging.debug(' Function:  isSameAllele - allele1:%s, allele2:%s' % (allele1, allele2))
        if allele1 == allele2:
            return 1
        else:
            return 0

               
