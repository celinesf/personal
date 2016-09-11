#!/usr/bin/env python

"""
     for each diseases in genetics.rep
     - for each snp, check the alleles are in + strand as defined by dbsnp
     09/13/13- 2.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging,  pymongo, copy
import datetime
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from SnpUtils import SnpUtils

DATE = datetime.datetime.now().strftime("%m-%d-%Y")

conn = pymongo.Connection(config.MONGO, 27017)
nextbio_db = conn["nextbio"]
genophen30 = conn["genophen30"]
genetics = conn["genetics"]

class FixGeneticsRepStrands():
    def __init__(self, disease):
        ''' '''
        self.util = NextBioUtils()
        self.snp_util = SnpUtils()
        self.disease = disease

    def fixData(self):
        logging.debug(' Function:  fixData' )
        genetic_data = genophen30['genetics.rep'].find({'disease':self.disease},{'_id':0})
        self.dbsnp_data = {}
        for doc in genetic_data:
            self.version = doc['version']
            self.date = doc['date']
            self.util.warnMe('info', ' fixing disease %s, version %s, date %s' % (self.disease,self.version,self.date))
            new_doc = copy.deepcopy(doc)
            new_doc['version'] =  str(int(float(new_doc['version']))+0.1)
            new_doc['date'] = DATE
            changed = False
            changed_snps = []
            for data in new_doc['snps']:
                self.snp = data['snp']
                if self.snp not in self.dbsnp_data:
                    dbsnp_data = nextbio_db['dbsnp'].find({'dbsnp':self.snp})
                    for snp_doc in dbsnp_data:
                        self.dbsnp_data[self.snp] = copy.deepcopy(snp_doc)
                if data['gt'] != "NA":
                    if len(data['gt']) == 2:
                        if data['gt'][0] not in self.dbsnp_data[self.snp]['alleles'] or data['gt'][1] not in self.dbsnp_data[self.snp]['alleles']:  
                            if self.snp not in changed_snps:
                                changed_snps.append(self.snp)
                            changed = True
                            gt = self.snp_util.reverseAlleles(data['gt'])
                            logging.debug(' changing snp %s, genotype %s to %s%s(fixData) - dbsnp alleles %s' %(self.snp,data['gt'], gt[0],gt[1], self.dbsnp_data[self.snp]['alleles']))
                            if gt[0] in self.dbsnp_data[self.snp]['alleles'] and gt[1] in self.dbsnp_data[self.snp]['alleles']:
                                data['gt'] =  "%s%s" % (gt[0],gt[1])
                            else:
                                self.util.warnMe('critical', ' ISSUE changing snp %s, genotype %s to %s%s(fixData) - dbsnp alleles %s' %(self.snp,data['gt'], gt[0],gt[1], self.dbsnp_data[self.snp]['alleles']))
                    elif len(data['gt']) == 1:
                        if data['gt'] not in self.dbsnp_data[self.snp]['alleles'] :  
                            if self.snp not in changed_snps:
                                changed_snps.append(self.snp)
                            changed = True
                            gt = self.snp_util.reverseAlleles(data['gt'])
                            logging.debug(' changing snp %s, genotype %s to %s(fixData) - dbsnp alleles %s' %(self.snp,data['gt'], gt[0], self.dbsnp_data[self.snp]['alleles']))
                            if gt in self.dbsnp_data[self.snp]['alleles'] :
                                data['gt'] =  gt[0]
                            else:
                                self.util.warnMe('critical', ' ISSUE changing snp %s, genotype %s to %s (fixData) - dbsnp alleles %s' %(self.snp,data['gt'], gt[0], self.dbsnp_data[self.snp]['alleles']))

                    else:
                        self.util.warnMe('critical', ' snp: %s genotype %s -- fixing disease %s, version %s, date %s ' % ( self.snp , data['gt'], self.disease,self.version,self.date))
            self.util.warnMe('info', ' I fixed %s SNPs -- fixing disease %s, version %s, date %s-- list: %s' % (len(changed_snps),self.disease,self.version,self.date,changed_snps))      
            if changed == True:
                self.util.writeOutput("genetics.rep_%s" %( self.disease), new_doc)
        
    def fixDataDrugs(self):
        logging.debug(' Function:  fixDataDrugs' )
        genetic_data = genophen30['genomics.rep'].find({'gid':self.disease},{'_id':0})
        self.dbsnp_data = {}
        for doc in genetic_data:  
            self.version = doc['version']
            self.date = doc['date']
            self.util.warnMe('info', ' fixing disease %s, version %s, date %s' % (self.disease,self.version,self.date))
            new_doc = copy.deepcopy(doc)
            new_doc['version'] =  str(int(float(new_doc['version']))+0.1)
            new_doc['date'] = DATE
            changed = False
            changed_snps = []
            for comp_data in new_doc[self.disease]:
                for d in comp_data['snp']:
                    self.snp = d['name']
                    for data in d['genotype']:
                        if self.snp not in self.dbsnp_data:
                            dbsnp_data = nextbio_db['dbsnp'].find({'dbsnp':self.snp})
                            for snp_doc in dbsnp_data:
                                self.dbsnp_data[self.snp] = copy.deepcopy(snp_doc)
                        if data['genotype'] != "NA":
                            if len(data['genotype']) == 2:
                                if data['genotype'][0] not in self.dbsnp_data[self.snp]['alleles'] or data['genotype'][1] not in self.dbsnp_data[self.snp]['alleles']:  
                                    if self.snp not in changed_snps:
                                        changed_snps.append(self.snp)
                                    changed = True
                                    gt = self.snp_util.reverseAlleles(data['genotype'])
                                    logging.debug(' changing snp %s, genotype %s to %s%s(fixData) - dbsnp alleles %s' %(self.snp,data['genotype'], gt[0],gt[1], self.dbsnp_data[self.snp]['alleles']))
                                    if gt[0] in self.dbsnp_data[self.snp]['alleles'] and gt[1] in self.dbsnp_data[self.snp]['alleles']:
                                        data['genotype'] =  "%s%s" % (gt[0],gt[1])
                                    else:
                                        self.util.warnMe('critical', ' ISSUE changing snp %s, genotype %s to %s%s(fixData) - dbsnp alleles %s' %(self.snp,data['genotype'], gt[0],gt[1], self.dbsnp_data[self.snp]['alleles']))
                            elif len(data['genotype']) == 1:
                                if data['genotype'] not in self.dbsnp_data[self.snp]['alleles'] :  
                                    if self.snp not in changed_snps:
                                        changed_snps.append(self.snp)
                                    changed = True
                                    gt = self.snp_util.reverseAlleles(data['genotype'])
                                    logging.debug(' changing snp %s, genotype %s to %s(fixData) - dbsnp alleles %s' %(self.snp,data['genotype'], gt[0], self.dbsnp_data[self.snp]['alleles']))
                                    if gt in self.dbsnp_data[self.snp]['alleles'] :
                                        data['genotype'] =  gt[0]
                                    else:
                                        self.util.warnMe('critical', ' ISSUE changing snp %s, genotype %s to %s (fixData) - dbsnp alleles %s' %(self.snp,data['genotype'], gt[0], self.dbsnp_data[self.snp]['alleles']))
        
                            else:
                                self.util.warnMe('critical', ' snp: %s genotype %s -- fixing disease %s, version %s, date %s ' % ( self.snp , data['genotype'], self.disease,self.version,self.date))
            self.util.warnMe('info', ' I fixed %s SNPs -- fixing disease %s, version %s, date %s-- list: %s' % (len(changed_snps),self.disease,self.version,self.date,changed_snps))      
            if changed == True:
                self.util.writeOutput("genetics.rep_%s" %( self.disease), new_doc)

        ''' getDbSnpDataFromMongo '''
    def getDbSnpDataFromMongo(self):
        logging.debug(' Function:  getDbSnpDataFromMongo-')
        db_data = nextbio_db['dbsnp'].find()
        self.dbsnp_database = {}
#         self.dbsnp_merged = {}
        for doc in db_data:
            if doc['dbsnp'] not in self.dbsnp_database:
                self.dbsnp_database[doc['dbsnp']] = doc
                # record merged documents 
#                 if doc['merged'] is not None:
#                     for snp in doc['merged']:
#                         if snp not in self.dbsnp_merged:
#                             self.dbsnp_merged[snp] = doc
#                         else: print "MERGED TWICE?"

    def ListMissing(self):
        self.getDbSnpDataFromMongo()
        
        snpor_snps = []
        snpor = genetics['snpinfo'].find({},{'dbSNP':1,'_id':0})
        for doc in snpor:
            if 'rs' in doc['dbSNP'] and doc['dbSNP'] not in self.dbsnp_database:
#                 dbsnp_data = nextbio_db['dbsnp'].find({'dbsnp':doc['dbSNP']})
#                 n = 0
#                 for d in dbsnp_data:
#                     n += 1
#                 if n == 0:
                print doc['dbSNP']
                snpor_snps.append(doc['dbSNP'])
        print snpor_snps


""" MAIN function """    
if __name__ == "__main__":
    
    disease = 'drugs'
    logging.basicConfig(filename='%s/%s%s_%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'fixGeneticRep',disease), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )
  
    clean = FixGeneticsRepStrands(disease)
#     clean.fixData()
#     clean.fixDataDrugs()
#     clean.ListMissing()
    
    print'DONE with FixGeneticsRepStrands'
    
""" END OF main """ 