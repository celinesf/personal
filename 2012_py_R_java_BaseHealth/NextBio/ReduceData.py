#!/usr/bin/env python

"""
     for Parkinson_disease
     bioset 5675 and 5614 are to large to manage.
     the full GWAS is there
     1) need to filter out the non-signficant snp from what is in nextbio_cook
     2) need to rewrite data with only significant snp not in nextbio_cook
     09/15/13- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging,  pymongo, copy, csv, string
import datetime
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from SnpUtils import SnpUtils


conn = pymongo.Connection(config.MONGO, 27017)
nextbio_db = conn["nextbio"]


class ReduceData():
    def __init__(self, bioset):
        ''' '''
        self.util = NextBioUtils()
        self.snp_util = SnpUtils()
        self.bioset = bioset

    def removeNonSignificantFromDB(self):
        data = nextbio_db["nextbio.cooked"].find({'bioset_id' :self.bioset, 'pvalue' :{'$ne':None}})
        for doc in data:
            if doc['pvalue'] > 0.05:
                self.util.warnMe('warning', 'removing %s %s' % (doc['dbsnp'], doc['pvalue']))
                nextbio_db["nextbio.cooked"].remove({'bioset_id' :self.bioset, 'dbsnp' :doc['dbsnp']})
      
    def getSignificantData(self):
        infile = open('%s/parkinsons_disease_data_v%s_QC.csv' % (config.DATAPATH, self.bioset),'rU')      
        data_file = csv.reader(infile)
        outf = open('%s/parkinsons_disease_data_v%s_REDUCE.csv' % (config.OUTPUTPATH, self.bioset),'wb')   
        out = csv.writer(outf)
        for row in data_file:
            if row[8] != "":
                if row[0] == self.bioset and float(row[8]) < 0.05:
                    snp = row[2]
                    data = nextbio_db["nextbio.cooked"].find({'bioset_id' :self.bioset, 'dbsnp' :snp})
                    nd =0
                    for doc in data:
                        nd+=1
                    if nd <2:
                        print 'signficant and missing snp', row[2], row[8]
                        out.writerow(row)                   
        infile.close()
        outf.close()
        


""" MAIN function """    
if __name__ == "__main__":
    
    bioset = '5675'
    logging.basicConfig(filename='%s/%s%s_%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'getSignificantData',bioset), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )
    
    clean = ReduceData(bioset)
#     clean.removeNonSignificantFromDB()
    clean.getSignificantData()
    
    print'DONE with ReduceData'
    
""" END OF main """ 