#!/usr/bin/env python

"""
##################### functions related to .bib files ###################
    06/11/13: version 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Main function extractDataFromBib
    2) functions to extract data from bib data: getKeyValueFromBib, getValueFromBib, addProtective
    3) quality/ sanity check functions: checkValueFromBib, checkRiskPhenotypeFromBib
"""

import logging
from AlgoGeneUtil import AlgoGeneUtil
import AlgoGeneConfig as config
from MleUtil import MleUtil


class BibUtil():
    def __init__(self, path,disease_name):
        self.path = path
        self.file_name = "%s/%s_pmid.bib" % (path,disease_name )
        self.disease_name = disease_name
        self.util = AlgoGeneUtil()
        self.paper_data = {} # data for each paper
        self.keys = [] # list of possible keys
    

##################### quality/ sanity check functions ###################
    ''' Check that bib value are 'na' and not empty or 'null' '''
    def checkValueFromBib(self, bib_key ,key,value):
        logging.debug(' Function:  checkValueFromBib - bib_key: %s ,key: %s,value: %s' % (bib_key ,key,value))
        if value.lower() == 'na':
            value = None
        elif value.lower() == 'null':
            self.warnMe('warning', ' NULL VALUE - GO CHECK: bib_key:%s, key:%s, value %s' % (bib_key,key,value))
            value = None   
        elif value =="":
            self.warnMe('warning', ' EMPTY VALUE - GO CHECK: bib_key:%s, key:%s, value %s' % (bib_key,key,value))
            value = None  
        elif (key == 'ci' or key == 'ranges' or key == 'cat' or key =='odds_ratio') and ',' in value:
            value = value.replace(', ', ',').lower()
        elif (key == 'risk_id' or key == 'disease_id') :
            value = value.replace(' ', '_').lower()
        return value
    
    ''' check that risk_id and disease_id fit with bib file name '''
    def checkRiskPhenotypeFromBib(self):
        logging.debug(' Function:  checkRiskPhenotypeFromBib')
        if self.paper_data['risk_id'] != self.risk_id:
            self.util.warnMe('warning', ' WRONG RISK ID - pmid:%s, risk:%s, found %s' % (self.paper_data['pmid'],self.risk_id, self.paper_data['risk_id']))            
        if self.paper_data['disease_id'] != self.disease_name:
            self.util.warnMe('warning', ' WRONG DISEASE ID - pmid:%s, disease:%s, found %s' % (self.paper_data['pmid'],self.disease_name, self.paper_data['disease_id']))    
#####################


##################### functions to extract data from bib data ###################
    ''' extract key and values of a line of bib file '''
    def getKeyValueFromBib(self,data, line):
        logging.debug(' Function:  getKeyValueFromBib - %s' % line)
        words = line.split(" = {")
        key = words[0].lower()
        value = words[1]
        
        ### extract genophen specific keys
        if key[0] == '_':
            k_list = key.split('_')
            nk = 2
            key = k_list[2]
            for nk in range(3,len(k_list)):
                key ='%s_%s' % (key,k_list[nk]) 
        ### extract values
        data, value = self.getValueFromBib(data, value)
        value = self.checkValueFromBib(self.paper_data['bib_key'],key,value)        
        ### add key/value to data
        if key not in self.keys:
            self.keys.append(key)
        self.paper_data[key]=value
        self.addProtective(key, value)
        logging.debug(' (getKeyValueFromBib) - key: %s \t %s ' % (key, value))
        return data
    
    ''' get Value from bib '''
    def getValueFromBib(self, data, value):
        logging.debug(' Function:  getValueFromBib - value: %s' % (value))
        if "}," in value: ## key value in one line
            value = value.split('},')[0]
        elif '}}' in value: ## key in one line+ end of paper data
            value = value.split('}}')[0]
        else: # value in several lines
            l = data.readline()
            l = l.strip()
            while len(l) >0 and ("}," not in l and '}}' not in l):
                value = '%s %s' %(value,l)
                l = data.readline()
                l = l.strip()
                if "}," in l: ## get last line of value
                    v = l.split('},')[0]
                    value = '%s %s' %(value,v)
                if '}}' in l:## get last line of value = end of paper data
                    v = l.split('}}')[0]
                    value = '%s\n%s' %(value,v)
        return data, value

    ''' add if risk factor is 'protective' '''
    def addProtective(self, key, value):
        logging.debug(' Function:  addProtective - key: %s,value: %s' % (key,value))
        if key == "normal_val":
            protective = "no"
            if value is not None:
                protective = "yes"
            if 'protective' not in self.keys:
                self.keys.append('protective')
            self.paper_data['protective']=protective
#####################
     

##################### Mane function called in GetModNonModRiskFactors ###################
    ''' extract all the data of a bib file 
        .bib files must be names <disease_name>_<risk_factor>.bib
    '''
    def extractDataFromBib(self):
        logging.debug(' Function:  extractDataFromBib')
        input_data = open(self.file_name,'r') # risk factor .bib data
        line = input_data.readline() ### first line ========
        data = []
        ### get to 1st paper data
        while '@' not in line and len(line)>0:
            line = input_data.readline()
        ### loop over paper data
        while '@'  in line and len(line)>0: ### found data for a paper
            bib_key = line.strip().split('{')[1].split(',')[0]
            self.util.warnMe('info', ' -> getting data for bib_key: %s ' % bib_key )
            self.paper_data = {'bib_key':bib_key}
            line = input_data.readline()
            line = line.strip()
            ## loop on key/value data
            while " = {" in line and len(line)>0:
                input_data = self.getKeyValueFromBib(input_data, line)
                line = input_data.readline()
                line = line.strip()
#             self.checkRiskPhenotypeFromBib()
            data.append(self.paper_data)
            line = input_data.readline()
        input_data.close() 

        mle = MleUtil()
        new_disease_citations = {}
        for ref in data:
            if ref['pmid'] not in new_disease_citations:
                new_disease_citations[ref['pmid']] = {"url": "http://www.ncbi.nlm.nih.gov//pubmed?term=%s" % ref['pmid']}
                new_disease_citations[ref['pmid']]['author'] = mle.getAuthors(ref)
                new_disease_citations[ref['pmid']]['content'] = mle.getContent(ref)
        print new_disease_citations   
        self.util.writeOutput('db_citation_%s' % ( self.disease_name), new_disease_citations)
#####################

################# MAIN ###############
if __name__ == '__main__':
#     disease_name = 'atrial_fibrillation'
#     disease_name = 'gallstone'
#     disease_name = 'allergic_rhinitis'
#     disease_name = 'gout'
#     disease_name = 'venous_thrombosis'
#     disease_name = 'alcohol_dependence'
#     disease_name = 'kidney_stone'
#     disease_name = 'melanoma'
#     disease_name = 'age_related_macular_degeneration'
    disease_name = 'rheumatoid_arthritis'

    logging.basicConfig(filename='%s/%s_%s_pmid' %(config.OUTPUTPATH, config.LOG_FILENAME,disease_name), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    BibUtil = BibUtil('%s/07-03-13_10newDiseases/%s' % (config.OUTPUTPATH, disease_name),disease_name)
    BibUtil.extractDataFromBib()
