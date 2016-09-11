#!/usr/bin/env python


"""
   Celine - 02/15/2013
   Need a script to insert large document in genetics.rep
"""


import pymongo
import json

# tiger_dev
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genophen30"]


class NextbioPoc():
    def __init__(self,filename ):
        self.directory = "/Users/celine/Box Documents/AlgoPhen/Genophen3.0/POC_Nextbio"
        self.filename = filename
     
    def insertDoc(self):
        print(' Function: insertDoc')   
        infile = open("%s/%s"%(self.directory,self.filename), 'r')
        for line in infile:
            newdoc = json.loads(line)

            ### add/update this data
            GENOPHEN30.genetics.rep.insert(newdoc)  

### main
###     
if __name__ == "__main__":
    
    filename = "t2d_copy"
    nextbioPOC = NextbioPoc("%s.json" % filename)
    nextbioPOC.insertDoc()
    print('I am done with %s ' % filename )