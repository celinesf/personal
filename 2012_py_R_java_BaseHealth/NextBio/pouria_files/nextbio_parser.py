#!/usr/bin/env python

"""
   NextBio batch paraser
   08/11/12 - 0.0.2 : second version after the first was unceremoniously deleted
   03/25/13 - 0.0.3 : picking up nextbio project with django frame work
"""
__author__ = "Pouria Mojabi"
__copyright__ = "Copyright 2012, Genophen.com"
__maintainer__ = "Pouria Mojabi"
__email__ = "mojabi@genophen.com"
__status__ = "dev" #

"""
   Algorithm & Flow
   
   To go through the flat text file from nextbio and tag things
"""

import json
import pymongo
import nltk

class nextbio_parser(object):
    
    def __init__(self, batchname):
        self.batchname = batchname      # batch file provided by nextbio
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        
        self.meta_sp   = [":", "="]     # this separates fields whithin meta data
        self.vs        = ["_vs_", "vs"] # this separates comparison and tag fields
                
        self.chuncks    = []             # aggregating all the sections of batch
        
     
    def mongo_init(self):
        connection = pymongo.Connection('50.112.142.17',27017)
        dbh = connection["admin"]   
        dbh.authenticate("adm1n","adm1npw0rd")
        dbh = connection["genetics"] # admin"]
        self.colxn     = dbh['nextbio_raw']
 
    def main(self):
        f = open(self.batchname)
        line = f.readline()
        while 1:
            line = f.readline()
            if not line:
                break
            chunk = {"snp":[]}
            
            # get meta data dynamically 
            while self.partition not in line:
                k=None; v=None
                for sep in self.meta_sp:
                    if sep in line:
                        try:
                            k, v = line.split(sep,1)
                        except Exception, e:
                            print line
                            print "Exception : ", str(e)
                if k : 
                    chunk[k.strip()] = v.strip()
                elif "\n" not in line :
                    print "Error Parsing: ", line
                
                # next line
                line = f.readline()
            
            # now line is all ==== so we ignore it!
            line = f.readline() # this line contains column names for snps
            ks   = line.split("\t")
            keys = [k.strip() for k in ks]
            snps = dict.fromkeys(keys)
            line = f.readline() # this line is just dash line
            line = f.readline() # this line is snp information
            
            # get snps dynamically:
            while self.partition not in line:
                s = line.split("\t")
                if len(s) != len(keys):
                    print "odd situation here, column names dont match the info"
                for i in range(len(s)):
                    snps[keys[i]] = s[i].strip()
                
                chunk["snp"].append(snps)
                    
                # next line
                line = f.readline()
                
            print json.dumps(chunk, indent=4)
            self.chuncks.append(chunk)
            
            line = f.readline()
   
    
if __name__ == "__main__":
    for batch in ["GWAS_standard1.txt"]:#, "GWAS_standard2.txt", "GWAS_standard3.txt", "GWAS_standard3_added.txt", "GWAS_standard4.txt", "GWAS_standard5.txt", "GWAS_standard6.txt", "GWAS_standard7.txt", "GWAS_standard8.txt", "GWAS_standard9.txt", "GWAS_standard10.txt", "GWAS_standard11.txt", "GWAS_standard12.txt" ]:
        n = nextbio_parser("/Users/pouria3/NextBio/" + batch)
        print batch
        n.main()
    