#!/usr/bin/python

# run beagle for all chromosome

import os
import sys
import MySQLdb #@UnresolvedImport

def get_dbsnp_snpor():
    db = MySQLdb.connect(host="",user="celin",passwd="celin",db="")
    cursor = db.cursor()
    dbsnp_snpor={}
    sql= "select distinct dbsnp from SNPOR where dbsnp like 'rs%' "
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        for s in tab: 
            dbsnp_snpor[s[0]]={}
    except: 
        print "Error: unable to fecth data %s" % sql  
    return dbsnp_snpor

if __name__ == '__main__':
    
    dbsnp_snpor= get_dbsnp_snpor()
    
    # get alleles from snpinfo for each of the dbsnp (indel and major/minor)
    # focus only on rs names
    # flag AT and CG snps without indel or aleles
    # get pmid/snp rows for snp
    
    # check allele_genotype and convert
    # does convert fit AGplus (output if not)
    #  special case for haplotype -> 
            # recover snp/alleles from each haplotype
            # and convert each allele/ genotypes
    
    print dbsnp_snpor
    
    
    pass