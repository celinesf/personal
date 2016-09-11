#!/usr/bin/env python

"""
     Convert xml dbsnp data into json
     06/14/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging, json, pymongo, copy
import DbSnpConfig as config


conn = pymongo.Connection(config.MONGO, 27017)
nextbio_db = conn["nextbio"]


class ReplaceMergedSnp():
    def __init__(self):
        ''' '''
        
    def getDbSnps(self):
        dbsnp_snp = nextbio_db['dbsnp'].find({}) 
        self.data_dbsnp = {}
        self.merged ={}
        d = 0
        for doc in dbsnp_snp:
            self.data_dbsnp[doc['dbsnp']] = doc
            if doc['merged'] is not None:
                for snp in doc['merged']:
                    if snp not in self.merged:
                        self.merged[snp] = doc
                    else:
                        print 'found another', snp
            d +=1
            if d % 10000 == 0:
                print d
    
    def replaceMerged(self):

        data_snp = nextbio_db['nextbio.cooked'].find({},{'dbsnp':1})
        for doc in data_snp:
            if doc['dbsnp'] in self.merged:
                nextbio_db['nextbio.cooked'].update({'dbsnp':doc['dbsnp']},{'$set':{'dbsnp':self.merged[doc['dbsnp']]['dbsnp']}})
                logging.info('replaceMerged  %s by %s - doc = %s ' %( doc['dbsnp'],self.merged[doc['dbsnp']]['dbsnp'],doc))
                print('replaceMerged  %s by %s - doc = %s ' %( doc['dbsnp'],self.merged[doc['dbsnp']]['dbsnp'],doc))

   

""" MAIN function """    
if __name__ == "__main__":


    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'_replace_mergered' ), filemode='w',
                        level=logging.INFO,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    rms = ReplaceMergedSnp()
    rms.getDbSnps()
    rms.replaceMerged()
    
    print'DONE with ReplaceMergedSnp'
    
""" END OF main """ 