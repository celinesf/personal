#!/usr/bin/env python

"""
    split GWAS_standard9.txt
   05/16/13 - 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) split GWAS_standard9.txt to allow reading   
"""

import NextBioConfig as config


""" MAIN function """    
if __name__ == "__main__":
    partition = "="*30         # this separates meta data from snp information
    snp_part  = "-"*30         # this separates column names from data in snp section
    fname = "GWAS_standard9"
    out_number = 0
    f = open("%s/%s%s" % (config.DATAPATH_ORI, fname,".txt"), 'r')
    out = open("%s/%s_%s%s" % (config.OUTPUTPATH,fname,out_number,".txt"), 'w')
    
    line = f.readline() 
    bioset_number = 0
    out.write(line)
    nl = 1
    while len(line)>0:
        ### Starting new bioset                    
        while partition in line and len(line)>0:     
            bioset_number += 1
            if nl >1000000 or (bioset_number % 50) == 0 :
                print out_number
                nl = 0
                out.close()
                out_number +=1
                out = open("%s%s_%s%s" % (config.OUTPUTPATH,fname,out_number,".txt"), 'w')
                out.write(line)
            line = f.readline() ### bioset title/ID
            nl +=1
            out.write(line)
            while partition not in line and len(line)>0 :
                line = f.readline()
                nl +=1
                out.write(line)
            ### starting SNP table
            if partition in line:
                line = f.readline() ### header ot table
                nl +=1
                out.write(line)
            line = f.readline() ### delimiter before SNP data
            nl +=1
            out.write(line)
            if snp_part in line:
                line = f.readline() ### 1st SNP info
                nl +=1
                out.write(line)
        line = f.readline()  
        nl +=1
        out.write(line)

    f.close()
    out.close()
    print'DONE with DMmain'
    
""" END OF main """ 