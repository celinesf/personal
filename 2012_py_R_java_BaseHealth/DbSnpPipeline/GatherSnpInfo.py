#!/usr/bin/env python

"""
     Perform quality control on GatherSnpInfo data
     06/15/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) check alleles strand
    2) check minor allele for + strand
    3) check reference allele on + strand
    4) check ancestral allele on +
"""

import logging,  copy
import DbSnpConfig as config
from DbSnpUtils import DbSnpUtils


class GatherSnpInfo():
    def __init__(self):
        self.util = DbSnpUtils()
        self.conflict = {}
        self.dbsnp_info = {}


    ''' checkStrand 
    '''
    def checkStrand(self, snp_data):
        logging.debug(' Function: checkStrand') 

        if snp_data['conflict'] == None:
            if snp_data['strand'] == "+":
                logging.info(' + strand %s (checkStrand)' % snp_data['dbsnp'])
            elif snp_data['strand'] == "-":
                self.util.warnMe('critical', ' - STRAND NO CONFLICT? %s (checkStrand)' % snp_data['dbsnp'])
                self.conflict[snp_data['dbsnp']] = copy.deepcopy(snp_data)
            else:
                self.util.warnMe('critical', ' STRAND UNKNOWN? %s (checkStrand)' % snp_data['dbsnp'])
                self.conflict[snp_data['dbsnp']] = copy.deepcopy(snp_data)
        elif 'RefAllele_conflict=yes' in  snp_data['conflict'] and 'snp_class_conflict=in-del' not in snp_data:
            if snp_data['strand'] == "+":
                self.util.warnMe('critical', ' + STRAND AND CONFLICT? %s (checkStrand)' % snp_data['dbsnp'])
                self.conflict[snp_data['dbsnp']] = copy.deepcopy(snp_data)
            elif snp_data['strand'] == "-":
                self.util.warnMe('warning', ' SWITCH STRAND %s (checkStrand)' % snp_data['dbsnp'])
            else:
                self.util.warnMe('critical', ' STRAND UNKNOWN AND CONFLICT? %s (checkStrand)' % snp_data['dbsnp'])    
                self.conflict[snp_data['dbsnp']] = copy.deepcopy(snp_data)
        elif 'snp_class_conflict=in-del' not in snp_data:
            self.nextbio_snps[snp_data['dbsnp']]['INDEL'] = 'yes'
            self.util.warnMe('warning', ' INDEL %s (checkStrand) - %s' % (snp_data['dbsnp'], snp_data))
            self.conflict[snp_data['dbsnp']] = copy.deepcopy(snp_data)
            


    ''' compareAlleles
    '''
    def compareAlleles(self, alleles1, alleles2):
        logging.debug(' Function: compareAlleles - alleles1: %s , alleles2: %s' % (alleles1, alleles2))  
        found = 0
#         if len(alleles1) == len(alleles2):
        for a1 in alleles1:
            if a1 in alleles2:
                found+=1
        if found == len(alleles1):
            return True
        else:
            return False
#         else:
#             self.util.warnMe('warning', ' ALLELES LENGHT DIFFER (compareAlleles) - alleles1: %s , alleles2: %s' % (alleles1, alleles2))  

    ''' reverse
    '''
    def reverseAlleles(self, alleles):
        logging.debug(' Function: reverseAlleles - %s' % alleles)  
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

    ''' confirmAlleles with Ss information
    '''
    def confirmAlleles(self, alleles, data):
        logging.debug(' Function: confirmAlleles')  
    
        data = self.dict2List(data)
        for doc in data:
            if self.compareAlleles(doc['Observed'].split('/') ,alleles) and doc['orient'] =='forward':
                if  doc['strand'] == 'top' :
                    logging.info(' confirmed alleles fwd/T (confirmAlleles) - %s' % alleles)
                elif  doc['strand'] is None :
                    logging.info(' confirmed alleles fwd (confirmAlleles) - %s' % alleles)
                elif doc['strand'] == 'bottom':
                    logging.warning(' COMFIRMED ALLELE? fwd/B (confirmAlleles) - %s' % alleles)
                else:
                    self.util.warnMe('critical',' COMFIRMED ALLELE? fwd WEIRD STRAND (confirmAlleles) - %s' % alleles)
            elif  self.compareAlleles(doc['Observed'].split('/') ,self.reverseAlleles(alleles)) and doc['orient'] =='reverse':
                if  doc['strand'] == 'bottom' :
                    logging.info(' confirmed alleles rev/B (confirmAlleles) - %s' % alleles)
                elif  doc['strand'] is None :
                    logging.info(' confirmed alleles rev (confirmAlleles) - %s' % alleles)
                elif doc['strand'] == 'top':
                    logging.warning(' COMFIRMED ALLELE? rev/T (confirmAlleles) - %s' % alleles)
                else:
                    self.util.warnMe('critical',' COMFIRMED ALLELE? rev WEIRD STRAND (confirmAlleles) - %s' % alleles)
            else:
                self.util.warnMe('warning',' SOMETHINGS IS WEIRD (confirmAlleles) - alleles: %s, doc: %s' % (alleles, doc))

    ''' dict2List
    '''
    def dict2List(self,  data):
        logging.debug(' Function: dict2List') 
        if type(data) != list:
            data=[data]
        return data
 
    ''' confirmReference with orientation and /or alleles
    '''
    def confirmReference(self, ref, data, strand, key):
        logging.debug(' Function: confirmReference - key:%s, ref:%s, strand: %s, data:%s' % (key, ref,strand, data))   
        if 'orientation' in data:
            if data['orientation'] == strand:
                return ref
            else:
                self.util.warnMe('error',' ISSUE ORIENTATION MISSMATCH (confirmReference) - key:%s, ref:%s, strand: %s, data:%s' % (key, ref,strand, data))
        else:
            logging.warning(' CANNOT CONFIRM REFERENCE (confirmReference) - key:%s, ref:%s, strand: %s, data:%s' % (key, ref,strand, data))
        if ref in self.snp_info['alleles']:
            return ref
        else:
            if key =='PrimarySequence':
                logging.warning(' REF/ALLELES MISSMATCH (confirmReference) -  key:%s, ref:%s, strand: %s, alleles:%s, data:%s' % (key,ref,strand, self.snp_info['alleles'],data)) 
            else:
                self.util.warnMe('critical',' REF/ALLELES MISSMATCH (confirmReference) -  key:%s, ref:%s, strand: %s, alleles:%s, data:%s' % (key,ref,strand, self.snp_info['alleles'],data)) 
            return None
   
    ''' checkReferenceStrand
    '''
    def checkReferenceStrand(self,  data, key):
        logging.debug(' Function: checkReferenceStrand - key %s' % key)         
        ref = None
        if "orient" in data:
            if data["orient"] == 'forward':
                ref = data['refAllele'] 
                ref = self.confirmReference(ref, data, 'fwd', key)
            elif data["orient"] == 'reverse':
                ref =self.reverseAlleles(data['refAllele'])
                ref = self.confirmReference(ref, data, 'fwd', key)
            else:
                self.util.warnMe('error',' CANNOT ORIENT MYSELF (checkReferenceStrand) - key %s, orient: %s' % (data["orient"], key))
        else:
            if key =='PrimarySequence':
                logging.warning(' CANNOT ORIENT MYSELF (checkReferenceStrand) - key %s, orient: %s' % (data["orient"], key))
            else:
                self.util.warnMe('critical',' CANNOT ORIENT MYSELF (checkReferenceStrand) - key %s, orient: %s' % (data["orient"], key))
        return ref

    ''' PrimarySequence information
    define reference allele
    '''
    def getReferenceAlleleFromPrimarySequence(self,  data):     
        logging.debug(' Function: getReferenceAlleleFromPrimarySequence')    
        ### confirm with PrimarySequence
        data['PrimarySequence'] = self.dict2List(data['PrimarySequence'])
        for doc in data['PrimarySequence'] :
            logging.info(' PrimarySequence (getReferenceAlleleFromAssembly) - doc: %s' % (doc))
            if doc['orient'] is not None and doc['refAllele'] is not None:
                new_ref = self.checkReferenceStrand(doc, 'PrimarySequence') 
                if new_ref != self.snp_info['reference']:
                    logging.warning(' REFERENCE MISSNATCH PrimarySequence and Assembly (getReferenceAlleleFromAssembly) - new_ref:%s, ref:%s' % (new_ref, self.snp_info['reference']))

    ''' getReferenceAlleleFromAssembly with assembly 
    find chromosome and position
    define reference allele
    '''
    def getReferenceAlleleFromAssembly(self,  data):
        data['Assembly'] = self.dict2List(data['Assembly'])
        ### loop on assembly data
        for doc in data['Assembly']:
            ### only consider reference information
            if doc['reference'] :
                ### loop on component data -> only one data is used
                component = doc['Component']
                component = self.dict2List(component)
                nc = 0
                for comp in component:
                    if comp['physMapInt'] is not None:
                        self.snp_info['reference'] = self.checkReferenceStrand(comp, 'Assembly') 
                        self.snp_info['chromosome'] =comp['chromosome'] 
                        self.snp_info['position'] = comp['physMapInt'] 
                        comp['FxnSet'] = self.dict2List(comp['FxnSet']) 
                        for func in comp['FxnSet']:
                            if func['fxnClass'] is not None:
                                if self.snp_info['functions'] is None:
                                    self.snp_info['functions'] =[]
                                self.snp_info['functions'].append(func)
                        nc+=1
                        break
                if nc >1:# -> only component one data is used
                    self.util.warnMe('critical',' SOMETHINGS IS WEIRD (getReferenceAlleleFromAssembly) - comp: %s' % (comp))

    ''' defineAncestralAllele
    '''
    def defineAncestralAllele(self, data):
        key = 'ancestral'
        logging.debug(' Function: defineAncestralAllele - key %s' %key)   
        if self.snp_info[key] is not None:
            alleles = self.snp_info[key].split(',')
            ancestral = alleles[len(alleles)-1]
            for allele in alleles:
                if allele != ancestral:
                    logging.warning(' %s ALLELES MISSNATCH (defineAncestralAllele) - key:%s' % (key,self.snp_info[key]))
                    break
            self.snp_info[key] = ancestral
        self.checkAlleleStrand('ancestral',data)
               
    ''' checkMinorAlleleStrand
    '''
    def checkAlleleStrand(self, key, data):
        logging.debug(' Function: checkAlleleStrand - key %s' %key)                  
        ### check allele strand
        if self.snp_info[key] not in self.snp_info['alleles'] and self.snp_info[key] is not None:
            rev = self.reverseAlleles(self.snp_info[key]) 
            if rev in self.snp_info['alleles']:
                logging.info(' reversed %s allele (checkAlleleStrand)'% key)
                self.snp_info[key] = rev   
            else:
                self.util.warnMe('critical',' %s ALLELE AND ALLELES MISSNATCH (checkAlleleStrand) - alleles:%s, key:%s' % (key, self.snp_info['alleles'],self.snp_info[key]))
        else:
            logging.info(' + strand %s allele (checkAlleleStrand)'% key)
    

    ''' getMerged rs ids
    '''
    def getMerged(self, data):
        logging.debug(' Function: getMerged') 
        data['MergeHistory'] = self.dict2List(data['MergeHistory']) 
        ndoc = 0
        for doc in data['MergeHistory'] :
            if doc['rsId'] is not None:
                if ndoc == 0:
                    self.snp_info['merged']=[]
                self.snp_info['merged'].append('rs%s' %doc['rsId'])
            ndoc+=1
            
    ''' defineStrandMinorFunction
    '''
    def defineStrandMinorFunction(self, data):
        logging.debug(' Function: defineStrandMinorFunction')    
        self.snp_info ={'dbsnp':'rs%s' %data['rsId'],
                   'class':data['snpClass'],
                   'alleles':data['Sequence']['Observed'].split('/'),
                   'ancestral':data['Sequence']['ancestralAllele'],
                   'minor_allele': data['Frequency']['allele'],
                   'maf': data['Frequency']['freq'],
                   'reference':None, 'functions' : None, 'chromosome':None, 'position':None, 'merged':None
                   }
        self.confirmAlleles(self.snp_info['alleles'], data['Ss'])
        self.getReferenceAlleleFromAssembly(data)
        self.getReferenceAlleleFromPrimarySequence(data)
        self.checkAlleleStrand('minor_allele',data)
        self.defineAncestralAllele(data)
        self.getMerged(data)
        return self.snp_info
   
    ''' gatherData 
    '''
    def gatherSnpData(self, data):
        logging.debug(' Function: gatherData') 
        cooked_data = self.defineStrandMinorFunction(data)
        
        return cooked_data
   
    ''' gatherData 
    '''
    def gatherData(self):
        logging.debug(' Function: gatherData') 
        self.clean_dbsnp={}
        for snp in self.dbsnp_data:
            self.util.warnMe('info', 'start for snp %s' %snp)
            self.clean_dbsnp[snp] = self.defineStrandMinorFunction(self.dbsnp_data[snp])
            
                       
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  GatherSnpInfo main' )

        self.dbsnp_data =self.util.getTermMap('%sdbsnp_simple_info_total.json' % config.OUTPUTPATH )
        self.gatherData()

        self.util.writeOutput('%s/dbsnp_clean.json' % config.OUTPUTPATH, self.json_data)   





""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'NextBio_GatherSnpInfo'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = GatherSnpInfo()
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 