#!/usr/bin/env python

"""
     Get data from unicode-cleaned nextbio batches
     05/28/13- 1.0
     06/06/13- 1.1
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Recover all keys/tags from batches 
    2) Recover all column names for SNP information tables from all batchs
    3) correct typos for tags/keys
    4) perform analysis regarding key/tag combinations to learn how to extract information
    5)  get 13 diseases to extract data by hand
    6) change checkQuote to change " to ' inside string -> data2
 
"""

import logging, os,  io, operator, xlwt
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from GetSnpList  import GetSnpList

class GetDiseaseList():
    def __init__(self):
        self.disease = None
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.issue_list = []
        self.util = NextBioUtils()

    
    ### find a term  as simple json value in a line  ###
    ### returns the term found in the line
    def findUnicode(self, term , line, words):
        logging.debug(' Function:  findUnicode term:%s\tword: %s\tline:%s' % (term ,words, line))
        if term.lower() in line.lower():
            words.append(term.lower())
        return words
    
    ### find a disease in information ###    
    def findDisease(self, line, key):
        logging.debug(' Function:  findDisease key: %s\tline:%s' % (key, line))
        disease = None
        line = line.strip().lower()
        if 'bioset title = ' in line:
            words = line.split('bioset title = ')
            if ' asso' in  words[1]:
                words = words[1].split(' asso')
            elif 'asso' in  words[1]:
                words = words[1].split('asso')
            disease = words[0]
        elif 'tag = ' in line:
            words = line.split('tag = ')
            if '|' in words[1]:
                words = words[1].split('|')
                if 'vs' in  words[0] or 'disease' in  words[0] or  'normal' in  words[0] or 'case' in  words[0] or 'control' in  words[0] : 
                    disease = words[1]
                else:
                    disease = words[0]
        elif 'comparison = ' in line:
            words = line.split('comparison = ')
            if ' _vs_ ' in words[1]:
                words = words[1].split(' _vs_ ')
                if 'control' in  words[0] or 'without' in  words[0] or  'healthy' in  words[0]: 
                    disease = words[1].replace('patients ','').replace(' patients','').replace('patients with ','')    
                else:
                    disease = words[0].replace('patients ','').replace(' patients','').replace('patients with ','')  
        
        return disease


    ### readBiosetMetaData ###
    def readBiosetMetaData(self,input_data, line):
        logging.debug(' Function:  readBiosetMetaData - %s' % line) 
        self.util.current_bioset['batch'] = self.util.batchname   
        self.util.current_bioset['info'] ={}
        self.util.key_combination = []
        for h in self.info_header:
            self.util.current_bioset['info'][h] = ""
            if h == 'fixed':
                self.util.current_bioset['info'][h] = []
        self.disease = []
        while self.partition not in line and len(line)>0 :  
            key = self.util.getKeyValueLower(self.key_map,line)
            ### TODO recover information as you read meta data
            if key is not None: 
                if  "BIOSET TITLE".lower() in key or "TAG".lower() in key or "COMPARISON".lower() in key:
                    self.disease.append(self.findDisease(line, key))

            line = input_data.readline()
        return line

    ### readSnpTableHeader ###   
    def readSnpTableHeader(self,input_data, line):
        logging.debug(' Function:  readSnpTableHeader ')   
        if self.partition in line:
            line = input_data.readline() ### header ot table
        line = input_data.readline() ### delimiter before SNP data
        return line


    ### readFirstSnpTableRow ###   
    # recover bioset id + pmid
    def readFirstSnpTableRow(self,input_data, line):
        logging.debug(' Function:  readFirstSnpTableRow ')  
        if self.snp_part in line: 
            line = input_data.readline() ### 1st SNP info
            self.util.getPMID(line)
            for disease in self.disease:
                if disease is not None:
                    if disease not in self.disease_data and disease is not None:
                        self.disease_data[disease]={'studies':[]}
                    if self.util.bioset_id not in self.disease_data[disease]['studies']:
                        self.disease_data[disease]['studies'].append(self.util.bioset_id)
       
        return line
     
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  GetDiseaseList main' )

        ### keys expected in nextbio batch
        self.key_map = self.util.getTermMap(config.KEYMAP)
        self.info_header = self.util.getTermMap(config.INFOHEADER)
        self.disease_data={}

#         ''' '''
        self.util.batchname = "gwasbiosets_aug2013_CB4.txt"
#         ''' '''
        print 'batch name',self.util.batchname   
        logging.debug(' Function:  GetDisease main - batch %s' % self.util.batchname ) 

        input_data = io.open("%s/%s" % (config.DATAPATH,self.util.batchname),'r',encoding='utf-8-sig')

        line = input_data.readline() ### first line ========
        while len(line)>0:
            ### Starting new bioset                    
            while self.partition in line and len(line)>0:
                line = input_data.readline() ### bioset title/ID
                logging.debug(' Function:  GetData4HandCurration main - title %s' % line )
                self.util.current_bioset ={}
                ### read meta information
                line = self.readBiosetMetaData(input_data, line)
                ### starting SNP table
                line = self.readSnpTableHeader(input_data, line)
                ## read first line of table
                line = self.readFirstSnpTableRow(input_data, line)
            
            line = input_data.readline()  
    
        output = open('%s/disease_list.xls' % config.DATAPATH2, 'w')
        for d in self.disease_data:
            nd = len(self.disease_data[d]['studies'])
            output.write('%s\t%s\n' % (d,nd))
        output.close
        self.util.writeOutput('disease_list' ,self.disease_data)

        

""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'disease_list'), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = GetDiseaseList()
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 