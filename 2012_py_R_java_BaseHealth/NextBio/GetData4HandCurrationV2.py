#!/usr/bin/env python

"""
     Get data from unicode-cleaned nextbio batches
     05/28/13- 1.0
     06/06/13- 1.1
     07/14/13- 2.0
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
    7) rewrite main function to recover array for each raw and data in format as expected by fixHandCurratedData
 
"""

import csv, json
import logging,   io, operator, re
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from ExcelUtils import ExcelUtils
from FixHandCurratedData import FixHandCurratedData


class GetData4HandCurration():
    def __init__(self):
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.issue_list = []
        self.util = NextBioUtils()
        self.excelUtil = ExcelUtils()

    ### returns the term found in the line
    def findUnicode(self, term , line, words):
        logging.debug(' Function:  findUnicode term:%s\tword: %s\tline:%s' % (term ,words, line))
        if term.lower() in line.lower():
            words.append(term.lower())
        return words
    
    ### check for term rejecting a disease ###
    def checkTerm(self,disease, term , line, words,accept):
        logging.debug(' Function:  checkTerm accept:%s\tterm:%s\tword: %s\tline:%s' % (accept,term ,words, line))
        ### term = value in json
        if type(term) == unicode:
            words = self.findUnicode(term.lower() , line.lower(), words)
            if term.lower() in words and accept == True:
                self.util.current_bioset["disease"]=disease
#             elif  term.lower() in words  and accept == False:
#                 logging.info("REJECT IN WORD: disease:%s\tterm: %s\twords:%s\tline%s" %(disease,term, words, line))
#                 if "FOUND_DISEASE" in self.util.current_bioset:
#                     self.util.current_bioset.pop("FOUND_DISEASE")
        ### term in a dictionary - not just a value       
#         elif  type(term) == dict:
#             words = self.findDictionary(disease, term , line.lower(), words,accept)
#         else:
#             logging.error("NOT dict nor unicode (checkTerm): disease:%s\tterm1: %s\twords:%s\tline%s" %(disease,term, words, line))

    ### find a disease in information ###    
    def findDisease(self, line, key):
        logging.debug(' Function:  findDisease key: %s\tline:%s' % (key, line))
        if "disease" not in self.util.current_bioset:
            for d in self.disease_map:
                for term in self.disease_map[d]:
                    if "disease" not in  self.util.current_bioset:
                        line = line.strip().lower()
                        if 'bioset title = ' in line:
                            words = line.split('bioset title = ')
                            if ' asso' in  words[1]:
                                words = words[1].split(' asso')
                            elif 'asso' in  words[1]:
                                words = words[1].split('asso')
                            disease = words[0]
                            self.checkTerm(d, term,disease,[],True)
                        elif 'tag = ' in line:
                            words = line.split('tag = ')
                            disease = words[1]
                            self.checkTerm(d, term,disease,[],True)
                        elif 'comparison = ' in line:
                            words = line.split('comparison = ')
                            disease = words[1]
                            self.checkTerm(d, term,disease,[],True)

    ### write in Nextbio style batch info ###
    def getMaps(self, ):
        logging.debug(' Function:  getMaps ')
        ######### find New SNPs to add to DBSBP collection
#         self.getSnpPList.main4disease()
        #########
        
        ### maps for NLP extraction of data
        self.ethnicity_map = self.util.getTermMap(config.ETHNICITYMAP)
        self.data_type_map = self.util.getTermMap(config.DATATYPEMAP)
        self.model_map = self.util.getTermMap(config.MODELMAP)
        self.info_header = self.util.getTermMap(config.INFOHEADER)
        self.meta_map = self.util.getTermMap(config.METAMAP)
        
        ### keys expected in nextbio batch
        self.key_map = self.util.getTermMap(config.KEYMAP)

        ### get terms for and rejecting diseases of interest
        self.disease_map = self.util.getTermMap(config.DISEASEMAP)
        self.disease_data={}
        for d in  self.disease_map:
            self.disease_data[d] = {"studies":[]}

    ### checkColumnAlligned ###
    def checkColumnAlligned(self, line):
        logging.debug(' Function:  checkColumnAlligned - %s' % line)
        snp_row = {}
        hapmap_pop_index = self.util.current_bioset['columns'].index('snp_population')
        line = line.strip()
        cells = line.strip().split('\t')
        ### check that snp table is correctly aline
        if len(cells) != len(self.util.current_bioset['columns']) and len(cells[hapmap_pop_index]) !=3 and cells[hapmap_pop_index].lower() != 'Global or other'.lower():
            self.util.warnMe('warning', 'different collumns numbers (checkColumnAlligned) - batch: %s, cells: %s, header:%s' % (self.util.current_bioset['batch'],cells, self.util.current_bioset['columns']))
            line = line.replace('\t\t',"\t")
            cells = line.strip().split('\t')
        for nc in range(0,len(self.util.current_bioset['columns'])):
            if nc >= len(cells):
                cells.append('')
            if cells[nc] == '':
                cells[nc] = None
            snp_row[self.util.current_bioset['columns'][nc].lower()] =  cells[nc]
        if 'POSSIBLE_INCONSISTENCY'.lower() in snp_row:
            if  snp_row['POSSIBLE_INCONSISTENCY'.lower()] is not None:
                if 'snp_issue' not in self.util.current_bioset:
                    self.util.current_bioset['snp_issue'] = {}
                    self.util.current_bioset['info']['issue'] = (snp_row['POSSIBLE_INCONSISTENCY'.lower()])
                elif snp_row['POSSIBLE_INCONSISTENCY'.lower()] not in self.util.current_bioset['info']['issue']  :
                    self.util.current_bioset['info']['issue'] = "%s|%s" %(self.util.current_bioset['info']['issue'],snp_row['POSSIBLE_INCONSISTENCY'.lower()])
                    self.util.warnMe('warning', 'Different issue for SNP: %s twice in bioset: %s - %s - %s' %(snp_row['snp'], self.util.bioset_id, self.util.current_bioset['info']['issue'], snp_row['POSSIBLE_INCONSISTENCY'.lower()] ))
                self.util.current_bioset['snp_issue'][snp_row['snp']] ={'minor_allele':None,'major_allele':None,'issue':['_odds_ratio', '95%ci','_frequency_case','_frequency_control']}
                if 'minor_allele' in snp_row  :
                    self.util.current_bioset['snp_issue'][snp_row['snp']]['minor_allele'] = snp_row['minor_allele']
                    self.util.current_bioset['snp_issue'][snp_row['snp']]['issue'].append('minor_allele')
                if 'major_allele' in snp_row:
                    self.util.current_bioset['snp_issue'][snp_row['snp']]['major_allele'] = snp_row['major_allele']
                    self.util.current_bioset['snp_issue'][snp_row['snp']]['issue'].append('major_allele')
        return snp_row

    ### extractDataType ###    
    def extractDataType(self, line):
        logging.debug(' Function:  extractDataType - \tline:%s' % (line))
        data_type = ""
        narrow_ethnicity = ""
        if 'asso' in  line:
            words = line.split('asso')
        for d_type in self.data_type_map:
            for d_find in d_type:
                if d_find in words[1]:
                    data_type = d_type[d_find]
                    narrow_ethnicity = words[1]
                    break
            if data_type != "":
                break
            
        if data_type == "":
            self.util.warnMe('warning', ' no data type (extractDataType) - line: %s' % line)
#         return None, words[1]
        self.util.current_bioset['info']['data_type'] = data_type
        self.util.current_bioset['info']['narrow_ethnicity'] = narrow_ethnicity

    ### extractModelType ###    
    def extractModelType(self, line):
        logging.debug(' Function:  extractModelType - \tline:%s' % (line))
        model_type = ""
        for m in self.model_map:
            if m in line.lower():
                model_type = self.model_map[m]
        self.util.current_bioset['info']['model_type'] = model_type

    ### extractMetaAnalysis ###    
    def extractMetaAnalysis(self, line):
        logging.debug(' Function:  extractMetaAnalysis - \tline:%s' % (line))
        is_meta = ""
        if 'yes' in line:
            is_meta = 'yes'
        elif 'no' in line:
            is_meta = "no"
        self.util.current_bioset['info']['is_meta'] = is_meta

    ### removeComaInSampleSizes ###    
    def removeComaInSampleSizes(self, line):
        logging.debug(' Function:  removeComaInSampleSizes - \tline:%s' % (line))

        words = re.split(',|and|&|-| ',line)
        ## transform 1,000 to 1000
        digit = 0
        for word in words:
            if word.isdigit() and digit == 0:
                prev = word
                digit +=1
            elif word.isdigit() and digit > 0:
                logging.info(' (removeComaInSampleSizes)- %s,%s' % (prev,word))
                line = line.replace('%s,%s' % (prev,word), "%s%s" %(prev,word))
            else:
                digit = 0
        return line
    
    ### extractSampleSizes ###    
    def extractSampleSizes(self, line):
        logging.debug(' Function:  extractSampleSizes - \tline:%s' % (line))
        line = self.removeComaInSampleSizes(line)
        case = 0
        control = 0
        ref = re.compile('\W+')
        words = ref.split(line)
        flag = 0
        for w in range(len(words)-1,-1, -1):
            if 'case' in words[w]:
                flag = 1
            elif 'control' in words[w]:
                flag = 2
            elif flag == 1 and words[w].isdigit():
                case += int(words[w])
                flag =0
            elif flag == 2 and words[w].isdigit():
                control += int(words[w])
                flag = 0
        self.util.current_bioset['info']['case_no'] = case
        self.util.current_bioset['info']['control_no'] = control
        self.util.current_bioset['info']['sample_size'] = case+control

    ### readBiosetMetaData ###
    def readBiosetMetaData(self,input_data, line):
        logging.debug(' Function:  readBiosetMetaData - %s' % line) 
        self.util.current_bioset['batch'] = self.batchname    
        self.util.current_bioset['info'] ={}
        self.util.key_combination = []
        for h in self.info_header:
            self.util.current_bioset['info'][h] = ""
            if h == 'fixed':
                self.util.current_bioset['info'][h] = []
        
        while self.partition not in line and len(line)>0 :  
            key = self.util.getKeyValueLower(self.key_map,line)
            ### TODO recover information as you read meta data
            if key is not None: 
                if  "BIOSET TITLE".lower() in key or "TAG".lower() in key or "COMPARISON".lower() in key:
                    self.findDisease(line, key)
                ### "combined", data_type, narrow_ethnicity
                if "BIOSET TITLE".lower() in key :
                    self.extractDataType(self.util.current_bioset[key].lower())
                    if self.util.current_bioset['info']['data_type'] != 'combined':
                        self.util.current_bioset['info']['combined'] = 'no'
                    else:
                        self.util.current_bioset['info']['combined'] = self.util.current_bioset['info']['narrow_ethnicity']
                ### model_type
                elif "BIOSET SUMMARY".lower() in key or "ANALYSIS SUMMARY" .lower() in key:
                    self.extractModelType(self.util.current_bioset[key].lower())
                ### case_no control_no sample_size  
                elif "SAMPLE NUMBER".lower() in key :
                    self.extractSampleSizes(self.util.current_bioset[key].lower()) # "case_no","control_no","sample_size",
                elif "META-ANALYSIS".lower() in key :
                    self.extractMetaAnalysis(self.util.current_bioset[key].lower())#"is_meta",
            line = input_data.readline()
        return line
    
    ### readSnpTableHeader ###   
    def readSnpTableHeader(self,input_data, line):
        logging.debug(' Function:  readSnpTableHeader ')   
        if self.partition in line:
            line = input_data.readline() ### header ot table
            if "disease" in self.util.current_bioset:
                self.util.current_bioset['snps'] = {}
                self.util.current_bioset['columns'] = line.lower().strip().split('\t')
        line = input_data.readline() ### delimiter before SNP data
        return line


    ### readFirstSnpTableRow ###   
    # recover bioset id + pmid
    def readFirstSnpTableRow(self,input_data, line):
        logging.debug(' Function:  readFirstSnpTableRow ')  
        if self.snp_part in line: 
            line = input_data.readline() ### 1st SNP info
            if "disease" in self.util.current_bioset:
                self.util.getPMID(line)
                snp_row = self.checkColumnAlligned(line)
                ### hapmap_pop, broad_ethnicity
                self.util.current_bioset['info']['hapmap_pop'] = snp_row['snp_population']
                self.util.current_bioset['info']['broad_ethnicity'] = self.ethnicity_map[snp_row['snp_population']]            
        return line

    ### resolvePotentialConflict ###
    def resolvePotentialConflict(self ):
        logging.debug(' Function:  resolvePotentialConflict ')
        self.fix_data = FixHandCurratedData(self.util.current_bioset['disease'])  
        self.fix_data.fixDataForCurration(self.util.current_bioset)
        self.util.current_bioset['snps'] = self.fix_data.cooked_data[0]['SNPS']
        self.util.current_bioset['columns'] = self.fix_data.cooked_data[0]['columns']

    ### writeCvsInfo ###
    def writeCvsInfo(self ):
        logging.debug(' Function:  writeCvsInfo ')
        self.disease_data[self.util.current_bioset['disease']]['studies'].append(self.util.current_bioset['bioset_id']) 
        if len(self.disease_data[self.util.current_bioset['disease']]['studies']) == 1:
            f1 = open("%s/%s_data_%s.csv" % (config.OUTPUTPATH,self.util.current_bioset['disease'],self.output),'wb')
            f2 = open("%s/%s_info_%s.csv" % (config.OUTPUTPATH,self.util.current_bioset['disease'],self.output),'wb')
            info_file = csv.writer(f2)
            info_file.writerow(self.info_header)
        else :
            f1 = open("%s/%s_data_%s.csv" % (config.OUTPUTPATH,self.util.current_bioset['disease'],self.output),'a')
            f2 = open("%s/%s_info_%s.csv" % (config.OUTPUTPATH,self.util.current_bioset['disease'],self.output),'a')
            info_file = csv.writer(f2)
        ###  write bioset info
        row = []
        for h in self.info_header:
            row.append(self.util.current_bioset['info'][h])
        info_file.writerow(row)
        f2.close()
        return f1

    ### writeCsvData ###
    def writeCsvData(self , f1):
        logging.debug(' Function:  writeCsvData ')
        data_file = csv.writer(f1)
        data_file.writerow([self.partition])
        for h in self.meta_map:
            data_file.writerow([h,self.util.current_bioset[h.lower()]])
        ### write snp table
        data_file.writerow([self.partition])
        data_file.writerow([element.upper() for element in self.util.current_bioset['columns']]) # column names
        data_file.writerow([self.snp_part])
        for snp in self.util.current_bioset['snps']: ## snp table
            row = []
            for c in self.util.current_bioset['columns']:
                if c in self.util.current_bioset['snps'][snp]:
                    row.append(self.util.current_bioset['snps'][snp][c])
                else :
                    row.append("")
            data_file.writerow(row)
        f1.close()

    ### recordFinishedBioset ###
    def recordFinishedBioset(self ):
        logging.debug(' Function:  recordFinishedBioset ')
        if len(self.util.current_bioset)>0:
            if "disease" in self.util.current_bioset:
                self.util.warnMe('info', 'bioset_id %s - disease: %s' % (self.util.current_bioset['bioset_id'],self.util.current_bioset['disease']) )
                self.util.current_bioset['info']['bioset_id'] = self.util.current_bioset['bioset_id']
                self.util.current_bioset['info']['pmid'] = self.util.current_bioset['pmid']
                ### fix data
                self.resolvePotentialConflict()
                
                ### output data in csv
                f1 = self.writeCvsInfo()
                ### write bioset meta data
                self.writeCsvData(f1)
               

    ''' find list of keys '''
    ''' find list of columns '''
    def main(self, batchname, output_name):
        logging.debug(' Function:  GetData4HandCurration main' )

        ### get maps
        self.getMaps()
                    
        ### open input file
        self.batchname = batchname  
        self.output = output_name
        self.util.warnMe('info', 'Function:  GetDisease main - batch %s' % batchname ) 
        input_data = io.open("%s/%s" % (config.DATAPATH,batchname),'r',encoding='utf-8-sig')
        
        ### read file line by line
        line = input_data.readline() ### first line ========
        while len(line)>0:
            ### Starting new bioset                    
            while self.partition in line and len(line)>0:                       
                self.recordFinishedBioset( )
                         
                line = input_data.readline() ### bioset title/ID
                logging.debug(' Function:  GetData4HandCurration main - title %s' % line )
                self.util.current_bioset ={}
                ### read meta information
                line = self.readBiosetMetaData(input_data, line)
                
                ### starting SNP table
                line = self.readSnpTableHeader(input_data, line)
                
                ### read first line of table
                line = self.readFirstSnpTableRow(input_data, line)
            ### record snp data for each row
            if "disease" in self.util.current_bioset:
                snp_row = self.checkColumnAlligned(line)
                self.util.current_bioset['snps'][snp_row['snp']] = snp_row
            line = input_data.readline()  

        ### dont forget to record the last bioset 
        if len(line) ==0:
            self.recordFinishedBioset( )
        
        summary = open("%s/Summary_%s.json" % (config.OUTPUTPATH,self.output),'w') 
        summary.write(json.dumps(self.disease_data, indent = 4))     



""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'NextBio_get_diseases'), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = GetData4HandCurration()
    input_file = "gwasbiosets_aug2013_3items_update_v8.txt"
    output_extention = "v8"
    clean.main(input_file, output_extention)
      

    print'DONE with DMmain'
    
""" END OF main """ 