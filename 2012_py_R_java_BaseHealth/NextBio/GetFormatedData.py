#!/usr/bin/env python

"""
   Format data for Mongodb from hand-currated disease information (13 diease for 3.1)
   06/04/13: version 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Get map for SNP table snp_table_header
    2) Recover data from excel file
    3) format data for mongoDB
"""

import logging,copy,pymongo, json, xlrd
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from SnpUtils import SnpUtils
from ExcelUtils import ExcelUtils

conn = pymongo.Connection(config.MONGO, 27017)
nextbio_db = conn["nextbio"]

class GetFormatedData():
    def __init__(self,disease_name):
        self.partition = "="*30         # this separates meta data from snp information
        self.dbsnp_part  = "-"*30         # this separates column names from data in snp section
        self.disease_name = disease_name
        self.util = NextBioUtils() 
        self.snp_util = SnpUtils()
        self.excel_util = ExcelUtils()

################ general functions  ############ 
    ''' setMinMajForPopulation 
    set minor/major alleles
    check that sum of freq is 1
    tag is '' or 'hapmap_' 
    '''
    def setMinMajForPopulation(self, tag, key, info, a, b):
        logging.debug(' Function:  setMinMajForPopulation - tag: %s, key: %s, a: %s, b: %s)' % (tag, key,a, b))
        print info
        print key, a, info[a][key]
        f1 = float(info[a][key])
        f2 = float(info[b][key])
        ### find min/max
        if  f1 > f2:
            self.setMinMajAlleles(b, a, tag)
        elif  f1 < f2:
            self.setMinMajAlleles(a, b, tag)
        else:
            logging.warning(' EQUAL FREQUENCIES (setMinMajForPopulation) - bioset: %s, snp: %s, a: %s, b: %s, key: %s' % (self.bioset,self.dbsnp,a,b, key))
        ### check sum to 1
        total = (f1 + f2)
        if total != 1.0:
            self.util.warnMe('warning', ' HAPMAP FREQUENCY !=1 (setMinMajForPopulation) - bioset: %s, snp: %s, f1: %s, f2: %s' %(self.bioset,self.dbsnp,info[a][key], info[b][key]))
        return info

    ### setMinMajAlleles ###
    def setMinMajAlleles(self,  mina, maxa,tag):
        logging.debug(' Function:  setMinMajAlleles - tag: %s, min: %s, max: %s' % (tag, mina, maxa))
        self.allele_data['%smajor_allele' %tag] = maxa
        self.allele_data['%sminor_allele' %tag] = mina
        self.allele_specific_info[mina]['%sis_minor'%tag] = True
        self.allele_specific_info[maxa]['%sis_minor'%tag] = False
        
    ''' setAndCheckIfExistingValue 
     check for fit between bioset_info and next_bio meta data info'''
    def setAndCheckIfExistingValue(self, key, data ,info):
        logging.debug(' Function:  setAndCheckIfExistingValue - key: %s' % key)  

        if data[key] is None :
            data[key] = info
        elif data[key] is not None and info is None and key == 'pvalue':
            logging.warning(' I KEEP NONE VALUE (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(self.bioset,self.dbsnp,key,data[key], info))
            data[key] = info
        elif data[key] is not None:
            fit = self.util.checkSame(data[key], info)
            if fit == False:
                if key == 'minor_allele' or key =='major_allele':
                    logging.warning(' DATA DONT FIT (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(self.bioset,self.dbsnp,key,data[key], info))
                elif key =='comment':
                    logging.info(' ADDING COMMENT (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(self.bioset,self.dbsnp,key,data[key], info))
                    data[key] = '%s. NEXTBIO: %s' %(data[key], info)
                elif key == 'hapmap_pop' and info == 'Global or other':
                    data[key] = info
                    logging.warning( ' hapmap_pop DONT FIT (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(self.bioset,self.dbsnp,key,data[key], info))
                elif key == 'ci' or key == 'or_value':
                    logging.warning( ' key: %s DONT FIT (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(key,self.bioset,self.dbsnp,key,data[key], info))
                else:
                    self.util.warnMe('critical', ' DATA DONT FIT (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(self.bioset,self.dbsnp,key,data[key], info))
        else:
            self.util.warnMe('critical', ' SHOULD NEVER HAPPEN (setAndCheckIfExistingValue) - bioset: %s, snp: %s,key: %s, data1: %s, data2: %s' %(self.bioset,self.dbsnp,key,data[key], info))

########################

 
################ DBSNP check and set functions ################ 
    ''' changeAlleleStrand  '''
    def changeAlleleStrand(self, snp_doc):
        logging.debug(' Function:  changeAlleleStrand - %s' % self.dbsnp) 
        snp_doc['strand'] ='+'
        rev = self.snp_util.reverseAlleles(self.db_data['alleles'])
        for key in snp_doc:
            if snp_doc[key] in rev:
                logging.warning(' REVERSING ALLELE (changeAlleleStrand) - snp: %s, genotype: %s, dbsnp: %s' %(self.dbsnp,snp_doc[key], self.db_data['alleles']))
                snp_doc[key] = self.snp_util.reverseAlleles(snp_doc[key])
    
    ''' changeGenotypeStrand  '''
    def changeGenotypeStrand(self, snp_doc):
        logging.debug(' Function:  changeGenotypeStrand - %s' % self.dbsnp) 
        key = 'hapmap_genotype_freq'
        if snp_doc[key] is not None:
            a = snp_doc[key].split('/')[0]
            if a not in self.db_data['alleles']:
                logging.warning(' REVERSING GENOTYPE (changeGenotypeStrand) - snp: %s, genotype: %s, dbsnp: %s' %(self.dbsnp,snp_doc[key], self.db_data['alleles']))
                for allele in self.db_data['alleles']:
                    if len(allele) == 1:
                        snp_doc[key] = snp_doc[key].replace(self.snp_util.reverseAlleles(allele),allele)
        elif 'global' in snp_doc['hapmap_pop'].lower():
            logging.warning(' GLOBAL NO GENOTYPE (changeGenotypeStrand) - bioset: %s, snp: %s, genotype: %s, dbsnp: %s' %(self.bioset, self.dbsnp,snp_doc[key], self.db_data['alleles']))
        else:
            print 'hapmap_pop', snp_doc['hapmap_pop']
            self.util.warnMe('warning',' NO GENOTYPE (changeGenotypeStrand) - bioset: %s, snp: %s, genotype: %s, dbsnp: %s' %(self.bioset, self.dbsnp,snp_doc[key], self.db_data['alleles']))
            


    ''' checkDbSnpStrand + compare to dbsnp data '''
    def checkDbSnpStrand(self, snp_doc):
        logging.debug(' Function:  checkDbSnpStrand - %s' % self.dbsnp) 
        if len(self.db_data['alleles']) ==2 and self.db_data['class'] == 'snp' and snp_doc['allele'] is not None:
            snp_doc['is_indel'] = False
            snp_doc['other_allele'] = self.db_data['alleles']
            if self.snp_util.compareAlleles(snp_doc['allele'], self.db_data['alleles']):
                snp_doc['strand'] ='+'
                if self.db_data['alleles'][0] == self.snp_util.reverseAlleles(self.db_data['alleles'][0]):
                    self.util.warnMe('warning', ' CHECK CG AT SNP (checkDbSnpStrand) - snp: %s, allele: %s, dbsnp: %s' %(self.dbsnp,snp_doc['allele'], self.db_data['alleles']))
            else:
                self.changeAlleleStrand(snp_doc)
        elif snp_doc['allele'] is None:
            self.changeGenotypeStrand(snp_doc)
        else:
#             self.util.warnMe('warning', ' NOT A SNP (checkDbSnpStrand) - snp: %s, allele: %s, dbsnp: %s, class:%s' %(self.dbsnp,snp_doc['allele'], self.db_data['alleles'],self.db_data['class']))
            if self.db_data['class'] == "in-del":
                self.util.warnMe('warning', ' INDEL (checkDbSnpStrand) - snp: %s, allele: %s, dbsnp: %s, class:%s' %(self.dbsnp,snp_doc['allele'], self.db_data['alleles'],self.db_data['class']))
                snp_doc['is_indel'] = True
                snp_doc['other_allele'] = self.db_data['alleles']
            elif len(self.db_data['alleles']) >2 :
                    logging.warning( ' MULTIPLE ALLELES (checkDbSnpStrand) - snp: %s, allele: %s, dbsnp: %s, class:%s' %(self.dbsnp,snp_doc['allele'], self.db_data['alleles'],self.db_data['class']))
                    snp_doc['is_indel'] = False
                    snp_doc['other_allele'] = self.db_data['alleles']
            else:
                self.util.warnMe('warning', ' NOT INDEL? (checkDbSnpStrand) - snp: %s, allele: %s, dbsnp: %s, class:%s' %(self.dbsnp,snp_doc['allele'], self.db_data['alleles'],self.db_data['class']))


    ''' checkDbSnpPosition  + assign chr and position in proper format'''
    def checkDbSnpPosition(self, snp_doc):
        logging.debug(' Function:  checkDbSnpPosition - %s' % self.dbsnp) 
        snp_doc['chromosome'] = snp_doc['position_37'].split(':')[0]
        snp_doc['position_37'] = snp_doc['position_37'].split(':')[1]
        if  snp_doc['position_36'].split(':')[0] == snp_doc['chromosome']:
            snp_doc['position_36'] = snp_doc['position_36'].split(':')[1]
        else:
            self.util.warnMe('critical', ' CHROMOSOME DONT MATCH  (checkDbSnpPosition) - snp: %s,chr37:%s, chr36:%s, data:%s ' %(self.dbsnp,snp_doc['position_36'].split(':')[0], snp_doc['chromosome'], snp_doc))   
        if self.db_data['chromosome'] is not None and self.db_data['position'] is not None: 
            if snp_doc['chromosome'] != self.db_data['chromosome']:
                self.util.warnMe('critical', ' CHROMOSOME != DBSNP (checkDbSnpPosition) - snp: %s, chr:%s, dbchr:%s, data:%s ' %(self.dbsnp, snp_doc['chromosome'] ,self.db_data['chromosome'], snp_doc))
            elif snp_doc['position_37'] != self.db_data['position']:
                self.util.warnMe('critical', ' POSITION != DBSNP (checkDbSnpPosition) - snp: %s, chr:%s, dbchr:%s, data:%s ' %(self.dbsnp, snp_doc['position_37'] ,self.db_data['position'], snp_doc))
        else:
            self.util.warnMe('warning', ' NO DBSNP INFO CHR/POS (checkDbSnpPosition) - snp: %s, chr:%s, dbchr:%s, data:%s ' %(self.dbsnp, snp_doc['position_37'] ,self.db_data['position'], snp_doc))

    ''' setDbSnpFunction  + assign chr and position in proper format'''
    def setDbSnpFunction(self, snp_doc):
        logging.debug(' Function:  setDbSnpFunction - %s' % self.dbsnp) 
        
        if self.db_data['functions'] is not None:
            snp_doc['function'] = None
            for doc in self.db_data['functions']:
                ### check gene name
                symbol = doc['symbol'].split('-AS1')[0]
                if symbol.lower()  not in snp_doc['gene_name'].lower():
                    logging.warning( ' ADDING DBSNP GENE (setDbSnpFunction) - snp: %s, gene:%s, dbGene: %s' %(self.dbsnp,snp_doc['gene_name'], doc['symbol'] ))
                ### add function
                if doc['fxnClass'] != 'reference':
                    if snp_doc['function'] is None:
                        snp_doc['function'] = doc['fxnClass']
                    else:
                        new_index = self.function_map.index(doc['fxnClass'])
                        old_index = self.function_map.index(snp_doc['function'])
                        if new_index < old_index:
                            snp_doc['function'] = doc['fxnClass']
             
                
    ''' getSnpSpecificDbSnpInfo + compare to dbsnp data '''
    def getSnpSpecificDbSnpInfo(self, snp_doc):
        logging.debug(' Function:  getSnpSpecificDbSnpInfo - %s' % self.dbsnp) 
        if self.db_data is not None:
            ### check alleles strand
            self.checkDbSnpStrand(snp_doc)
            ### set check positions/chromosomes
            if snp_doc['chromosome'] is None:
                self.checkDbSnpPosition(snp_doc)
            #### set functions
            if 'function' not in snp_doc:
                self.setDbSnpFunction(snp_doc)                
        else:
            self.util.warnMe('critical', ' NO DBSNP (getSnpSpecificDbSnpInfo) - snp: %s' %(self.dbsnp))
################


################ functions called by formatSnp4Mongo ################
    ''' setTemplateFormat '''
    def setTemplateFormat(self, empty_doc, snp_data):
        logging.debug(' Function:  setTemplateFormat - %s' % snp_data['snp'])  
        for key in self.map:
            key2= None
            if self.map[key] is not None :
                if type(self.map[key]) is not list:
                    key2 = self.map[key].lower()
                else:
                    key2 = self.map[key] 
            if key in snp_data:
#                 print 'a', snp_data[key], empty_doc
                self.setAndCheckIfExistingValue(key,empty_doc,snp_data[key])
            elif type(key2) is list:
                for k in  key2:
                    if k.lower() in snp_data:
#                         print 'b', snp_data[k.lower()], empty_doc
                        self.setAndCheckIfExistingValue(key,empty_doc,snp_data[k.lower()])
            elif key2 in snp_data:
#                 print 'c', key, key2, snp_data[key2], empty_doc
                self.setAndCheckIfExistingValue(key,empty_doc,snp_data[key2])
        self.getSnpSpecificDbSnpInfo(empty_doc)
 
    ''' fillAlleleDocument 
    get data for each key in map
    '''
    def fillAlleleDocument(self, allele, dbsnp_doc, snp_data):
        logging.debug(' Function:  fillAlleleDocument - %s' % snp_data['snp'])  
        for key in self.map:
            key2= None
            if self.map[key] is not None:
                if type(self.map[key]) is not list:
                    key2 = self.map[key].lower()
                else:
                    key2 = self.map[key] 
            ### data for specific allele
#             print 'here', allele, key, key2, self.allele_specific_info[allele], dbsnp_doc
            if key in self.allele_specific_info[allele]:
#                 print 'a', key, key2
                self.setAndCheckIfExistingValue(key,dbsnp_doc, self.allele_specific_info[allele][key])
            elif type(key2) is list:
                for k in  key2:
                    if k.lower() in self.allele_specific_info[allele] :
#                         print 'b', key, k,self.allele_specific_info
                        self.setAndCheckIfExistingValue(key,dbsnp_doc,self.allele_specific_info[allele][k.lower()])
            elif key2 in self.allele_specific_info[allele]:
#                 print 'c', key, key2
                self.setAndCheckIfExistingValue(key,dbsnp_doc,self.allele_specific_info[allele][key2])
            ### data about snp ref/risk/protectve..
            if key in self.allele_data:
#                 print 'd', key, key2
                self.setAndCheckIfExistingValue(key,dbsnp_doc,self.allele_data[key])
            elif type(key2) is list:
                for k in  key2:
                    if k.lower() in self.allele_data:
#                         print 'e', key, k
                        self.setAndCheckIfExistingValue(key,dbsnp_doc,self.allele_data[k.lower()])
            elif key2 in self.allele_data:
#                 print 'f', key, key2
                self.setAndCheckIfExistingValue(key,dbsnp_doc,self.allele_data[key2])
            ### orphan key with noe data
            if dbsnp_doc[key] is None:
                logging.warning(' ORPHAN (createDocPerAllele) - bioset: %s, snp: %s, key: %s, key2: %s' %(self.bioset,self.dbsnp,key, key2))
        ### fix strand
        self.getSnpSpecificDbSnpInfo(dbsnp_doc)
        
    ''' checkMongoDocuments 
    check that documents for bth alleles have been added
    '''
    def checkMongoDocuments(self, bioset, snp):
        logging.debug(' Function:  checkMongoDocuments - bioset: %s, snp:%s' % (bioset,snp))
        cursor = nextbio_db['nextbio.cooked'].find({'dbsnp':snp,'bioset_id':bioset })   
        nd =0
        for doc in cursor:
            nd+=1
        if nd != 2:
            self.util.warnMe('warning', ' MONGO ISSUE (checkMongoDocuments) - snp:%s, nd:%s' % ( snp,nd))
######################## 

######################## functions called by findRefRiskAllele ####################
    ''' getAlleleSpecificNextBioData 
    return allele specific information such as is-minor/ is reference , OR, pvalue, ...
    '''
    def getAlleleSpecificNextBioData(self, snp_data):
        logging.debug(' Function:  getAlleleSpecificNextBioData %s' % snp_data['snp'])
        doc = {snp_data["allele_1"]:{'is_risk':None,'is_reference':None,'is_minor':None,'hapmap_is_minor':None, 'pvalue':snp_data['p-value']},
               snp_data["allele_2"]:{'is_risk':None,'is_reference':None,'is_minor':None,'hapmap_is_minor':None, 'pvalue':snp_data['p-value']} }
        for key in snp_data:
            if key not in doc[snp_data["allele_1"]] and 'allele_1' in key:
                for a in range(1,3):
                    tmp_key = key.replace('allele_1','allele_%s' %a)
                    new_key = key.replace('allele_1_','')
                    new_key = new_key.replace('allele_1','allele')
                    if tmp_key in snp_data:
                        doc[snp_data["allele_%s" %a]][new_key] = snp_data[tmp_key]
                    else:
                        doc[snp_data["allele_%s" %a]][new_key] = None
        return doc
    
    ''' inferFromEffectAndOR 
    define what alleles is reference vs risk/protective'''
    def inferFromEffectAndOR(self, allele_info, key, a, b, limit):
        logging.debug(' Function:  inferFromEffectAndOR - snp: %s, a: %s, b: %s, key: %s, limit: %s' % (self.dbsnp, a,b, key,limit))
        self.allele_data["reference_allele"] = b
        self.allele_data["data_allele"] = a
        allele_info[b]['is_reference'] = True
        allele_info[a]['is_reference'] = False
        allele_info[b]['pvalue'] = None
        allele_info[b][key] = 1
        if allele_info[a][key] is not None:
            try :
                if   float(allele_info[a][key]) <limit:
                    allele_info[b]['is_risk'] = True
                    allele_info[a]['is_risk'] = False
                elif float(allele_info[a][key]) >limit:
                    allele_info[a]['is_risk'] = True
                    allele_info[b]['is_risk'] = False
                else:
                    logging.warning(' NO RISK/PROTECTIVE? (inferFromEffectAndOR) -  bioset: %s, snp: %s, a: %s, b: %s, key: %s' % (self.bioset,self.dbsnp,a,b, key))
            except:
                if allele_info[a][key] == '<1' :
                    allele_info[a][key] = None
                    allele_info[b]['is_risk'] = True
                    allele_info[a]['is_risk'] = False
                    self.util.warnMe('warning',' OR <1  (inferFromEffectAndOR) -  bioset: %s, snp: %s, a: %s, b: %s, key: %s' % (self.bioset,self.dbsnp,a,b, key))
                elif allele_info[a][key] == '>1' :
                    allele_info[a][key] = None
                    allele_info[a]['is_risk'] = True
                    allele_info[b]['is_risk'] = False
                    self.util.warnMe('warning', ' OR > 1  (inferFromEffectAndOR) -  bioset: %s, snp: %s, a: %s, b: %s, key: %s' % (self.bioset,self.dbsnp,a,b, key))
                else:
                    self.util.warnMe('critical',' WHAT AM I ? (inferFromEffectAndOR) -  bioset: %s, snp: %s, a: %s, b: %s, key: %s' % (self.bioset,self.dbsnp,a,b, key))
        else:
            logging.warning(' NO OR (inferFromEffectAndOR) -  bioset: %s, snp: %s, a: %s, b: %s, key: %s' % (self.bioset,self.dbsnp,a,b, key))                   


    ''' getAlleleSpecificInferedData 
    define what alleles is reference vs risk/protective'''
    def getAlleleSpecificInferedData(self, allele_info, key, a, b):
        logging.debug(' Function:  getAlleleSpecificInferedData - snp: %s, a: %s, b: %s, key: %s' % (self.dbsnp, a,b, key))
        ## define ref/risk/protective allele
        if 'odds' in key and 'ratio' in key and "(" not in key :
#             print key
            self.inferFromEffectAndOR(allele_info,key,  a, b, 1)
        if 'effect' in key or 'z-score' in key or 'beta' in key:
            self.inferFromEffectAndOR( allele_info,key,  a, b, 0)
        ### caluclate study  frequencies
        elif 'frequency_' in key and allele_info[a][key] is not None:
            allele_info[b][key] = str(1.0-float(allele_info[a][key]))
            if 'control' in key :
                self.setMinMajForPopulation('', key, allele_info, a, b)

########################
        
################ functions called by formatSnp4Mongo ############  
    ''' getDbSnpDataFromMongo '''
    def getDbSnpDataFromMongo(self):
        logging.debug(' Function:  getDbSnpDataFromMongo-')
        db_data = nextbio_db['dbsnp'].find()
        self.dbsnp_database = {}
        for doc in db_data:
            if doc['dbsnp'] not in self.dbsnp_database:
                self.dbsnp_database[doc['dbsnp']] = doc
        
            
    ''' findRefRiskAllele 
    self.allele_data contain snp specific data like
    reference, risk minor/major alleles...
    '''
    def findRefRiskAllele(self, snp_data):
        logging.debug(' Function:  findRefRiskAllele - %s' % snp_data['snp'])
        self.allele_data={
                          "data_allele":None,"reference_allele":None,
                          'hapmap_minor_allele':None, 'hapmap_major_allele':None, 
                          'major_allele':None,'minor_allele':None}
        ### recover data specific to each allele
        self.allele_specific_info = self.getAlleleSpecificNextBioData(snp_data)
        a1 = snp_data["allele_1"]
        a2 = snp_data["allele_2"]
#         print a1, a2, self.allele_specific_info
        ### find risk/ref -> fill in the blanks
        for key in self.allele_specific_info[a1]:
            ### data is for allele_1 
            if self.allele_specific_info[a1][key] is not None and self.allele_specific_info[a2][key] is None:
                self.getAlleleSpecificInferedData(self.allele_specific_info, key, a1, a2)
            ### data is for allele_2
            elif self.allele_specific_info[a1][key] is  None and self.allele_specific_info[a2][key] is not None:
                self.getAlleleSpecificInferedData(self.allele_specific_info, key, a2, a1)
            ### data given for both 
            elif  self.allele_specific_info[a1][key] is not None and self.allele_specific_info[a2][key] is not None:
                if 'frequency' in key:
                    self.setMinMajForPopulation('hapmap_', key, self.allele_specific_info, a1, a2)
            ### data given not given 
            elif  self.allele_specific_info[a1][key] is  None and self.allele_specific_info[a2][key] is  None:
                logging.warning(' NULL DATA (findRefRiskAllele) -  bioset: %s, snp: %s, key: %s' %(self.bioset,self.dbsnp,key))  
                self.getAlleleSpecificInferedData(self.allele_specific_info, key, a1, a2) ### data come from 1
            
    ''' createDocPerAllele '''
    def createDocPerAllele(self, empty_doc, snp_data):
        logging.debug(' Function:  createDocPerAllele - %s' % snp_data['snp'])  
        ### set template data
        self.setTemplateFormat(empty_doc, snp_data)
        self.dbsnp_doc = {'reference_allele':copy.deepcopy(empty_doc),'data_allele':copy.deepcopy(empty_doc)}   
        ### clean mongo
        snp = self.dbsnp
        ### fill data for each allele
        for role in self.dbsnp_doc:
            self.fillAlleleDocument(self.allele_data[role], self.dbsnp_doc[role], snp_data)
            ### record for documentation
            self.formated_data.append(self.dbsnp_doc[role])
            if 'merged_in' in self.db_data:
                self.util.warnMe('warning', ' change merged snp %s into %s(createDocPerAllele) - bioset: %s, snp: %s' % (snp_data['snp'],self.db_data['merged_in'], self.bioset, snp_data['snp']))         
                self.dbsnp_doc[role]['dbsnp'] = self.db_data['merged_in']
                snp = self.dbsnp_doc[role]['dbsnp']
            ### add to mongo
            nextbio_db['nextbio.cooked'].insert(copy.deepcopy(self.dbsnp_doc[role]))
        self.checkMongoDocuments(self.bioset,snp)
        ### check mongo ok
        
########################



################ functions by formatData4Mongo ############ 
    ''' initTemplateFormat '''
    def initTemplateFormat(self):
        logging.debug(' Function:  initTemplateFormat - bioset: %s' % self.bioset)
        template = {'nextbio':{}}
        
        for key in self.map :
            template[key] = None
        template["disease_id"] = self.disease_name
        ### fill template with bioset data
        for key in self.nextbio_data[self.bioset]:
            if key not in self.map :
                if key != 'snps' and key != "bioset_id":
                    template['nextbio'][key] =self.nextbio_data[self.bioset][key]
            else :
                if key =='comment':
                    template['nextbio'][key] =self.nextbio_data[self.bioset][key]
                else:
                    template[key] = self.nextbio_data[self.bioset][key]
        return template


    ''' formatSnp4Mongo '''
    def formatSnp4Mongo(self, template, ns):
        logging.debug(' Function:  formatSnp4Mongo-  bioset: %s' % self.bioset)
        ### loop on snp
        for snp in self.nextbio_data[self.bioset]['snps']:
            logging.info(' bioset: %s, snp: %s (formatData4Mongo) ' % (self.bioset, snp))
            self.dbsnp = snp
            ### find dbsnp data
            self.db_data = self.dbsnp_database[snp]
#             self.getDbSnpDataFromMongo(snp)
            ### obtain self.allele_data, self.allele_specific_info
            self.findRefRiskAllele(self.nextbio_data[self.bioset]['snps'][snp])
#             print self.allele_specific_info
            ### create 2 doc for reaf and risk allele
            self.createDocPerAllele(copy.deepcopy(template), self.nextbio_data[self.bioset]['snps'][snp])
            
            ns+=1
            if ns % 100 == 0:
                print 'snp#', ns
################   
     
    ''' formatData4Mongo all bioset data '''
    def formatData4Mongo(self):
        logging.debug(' Function:  formatData4Mongo' )
        self.formated_data = []
        ns =0
        for bioset in self.nextbio_data:
#             print 'a', bioset, self.nextbio_data[bioset]
            self.bioset = bioset
            template = self.initTemplateFormat()
            self.formatSnp4Mongo(template, ns)
#             print 'b', bioset, self.formated_data
 
    ''' isBiosetPartOfCombined '''
    def isBiosetPartOfCombined(self):
        logging.debug(' Function:  isBiosetPartOfCombined' )
        for bioset in self.nextbio_data:
            if self.nextbio_data[bioset]['combined'] != 'no':
                metas = self.nextbio_data[bioset]['combined'].strip().split(',')
                for b in metas:
                    if 'in_meta' not in self.nextbio_data[b]:
                        self.nextbio_data[b]['in_meta']=[]
                    self.nextbio_data[b]['in_meta'].append(bioset)
    
    def main(self):
        logging.debug(' Function:  GetFormatedData main %s' % self.disease_name)
        ''' get Template + snp table map'''
        self.map = self.util.getTermMap("%s/maps/%s" % (config.OUTPUTPATH , "hand_curration_map.json")) 
        self.function_map = self.util.getTermMap(config.FUNCTIONOMAP)
        nextbio_db['nextbio.cooked'].remove({"disease_id":self.disease_name })
        self.getDbSnpDataFromMongo()
        ### read cooked data
        self.nextbio_data, self.bioset_info, self.snp_table_header, self.column_key= self.excel_util.readExcelData("%s/cooked_%s.xls" % (config.DATAPATH2,disease_name))
        print self.nextbio_data
        self.isBiosetPartOfCombined()
        ### format data
        self.formatData4Mongo()
        print self.formated_data
        ### record json for documentation
        self.util.writeOutput('%s_formated' % self.disease_name,self.formated_data)
        

""" MAIN function """    
if __name__ == "__main__":
    ## age_related_macular_degeneration , gallstone, venous_thrombosis, atrial_fibrillation, allergic_rhinitis,gout, kidney_stone,alcohol_dependence,melanoma,rheumatoid_arthritis
    disease_name = "allergic_rhinitis"
    logging.basicConfig(filename='%s/%s%s_%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'Format',disease_name), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )
    print'Function: __main__', disease_name
    clean = GetFormatedData(disease_name)
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 