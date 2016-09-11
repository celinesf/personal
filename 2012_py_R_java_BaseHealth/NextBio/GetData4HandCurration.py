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

import csv
import logging,   io, operator, re
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from ExcelUtils import ExcelUtils
from FixHandCurratedData import FixHandCurratedData

######### find New SNPs to add to DBSBP collection
# from GetSnpList  import GetSnpList
######### 

class GetData4HandCurration():
    def __init__(self):
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.issue_list = []
        ######### find New SNPs to add to DBSBP collection
#         self.getSnpPList = GetSnpList()
        ######### 
        self.util = NextBioUtils()
        self.excelUtil = ExcelUtils()
        
        
      
    ### find a term in a dicionary of possible teamrs ###
    def findDictionary(self,disease, term , line, words,accept):
        logging.debug(' Function:  findDictionary accept:%s\tterm:%s\tword: %s\tline:%s' % (accept,term ,words, line))
        for t in term:
            if type(t) == unicode:
                words = self.findUnicode(t , line, words)
            if type(term[t]) == list and t.lower() in words:
                for t2 in term[t]:
                    if type(t2) == unicode:
                        words = self.findUnicode(t2.lower() , line.lower(), words)     
                        if t2.lower() in words  and accept == True:
                            logging.info("FOUND IN DICT: disease:%s\tterm1: %s\tterm2: %s\twords:%s\tline%s" %(disease,t, t2, words, line))
                            self.util.current_bioset["FOUND_DISEASE"]= {"disease":disease, "words":words, "line":line}    
                        elif  t2.lower() in words  and accept == False:
                            logging.info("nREJECT IN DICT: disease:%s\tterm1: %s\tterm2: %s\twords:%s\tline%s" %(disease,t, t2, words, line))
                            if "FOUND_DISEASE" in self.util.current_bioset:
                                self.util.current_bioset.pop("FOUND_DISEASE")
                            break
                    elif  type(t2) == dict:
                        words = self.findDictionary(disease, t2, line.lower(), words,accept)  
                    else: 
                        logging.error("NOT dict nor unicode (findDictionary): disease:%s\tterm1: %s\tterm2: %s\twords:%s\tline%s" %(disease,t, t2, words, line))
            elif type(term[t]) != list and t.lower() in words:
                logging.error("NOT list nor unicode (findDictionary): disease:%s\tterm1: %s\tterm2: %s\twords:%s\tline%s" %(disease,term, t, words, line))
        return words
    
    ### find a term  as simple json value in a line  ###
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
                self.util.current_bioset["FOUND_DISEASE"]={"disease":disease, "words":words, "line":line}
            elif  term.lower() in words  and accept == False:
                logging.info("REJECT IN WORD: disease:%s\tterm: %s\twords:%s\tline%s" %(disease,term, words, line))
                if "FOUND_DISEASE" in self.util.current_bioset:
                    self.util.current_bioset.pop("FOUND_DISEASE")
        ### term in a dictionary - not just a value       
        elif  type(term) == dict:
            words = self.findDictionary(disease, term , line.lower(), words,accept)
        else:
            logging.error("NOT dict nor unicode (checkTerm): disease:%s\tterm1: %s\twords:%s\tline%s" %(disease,term, words, line))
        
    ### find a disease in information ###    
    def findDisease(self, line, key):
        logging.debug(' Function:  findDisease key: %s\tline:%s' % (key, line))
        if "FOUND_DISEASE" not in self.util.current_bioset:
            for d in self.disease_map:
                for term in self.disease_map[d]:
                    if "FOUND_DISEASE" not in  self.util.current_bioset:
                        line = line.strip().lower()
                        words = line.split('bioset title = ')
                        if ' asso' in  words[1]:
                            words = words[1].split(' asso')
                        elif 'asso' in  words[1]:
                            words = words[1].split('asso')
                        disease = words[0]
                        self.checkTerm(d, term,disease,[],True)

    ### extractMetaAnalysis ###    
    def extractMetaAnalysis(self, line):
        logging.debug(' Function:  extractMetaAnalysis - \tline:%s' % (line))
        if 'yes' in line:
            return 'yes'
        elif 'no' in line:
            return "no"
        else:
            return ""

    ### extractModelType ###    
    def extractModelType(self, line):
        logging.debug(' Function:  extractModelType - \tline:%s' % (line))
        for m in self.model_map:
            if m in line.lower():
                return self.model_map[m]
        return ""

    ### extractDataType ###    
    def extractDataType(self, line):
        logging.debug(' Function:  extractDataType - \tline:%s' % (line))
        if 'asso' in  line:
            words = line.split('asso')
        for data_type in self.data_type_map:
            for d in data_type:
                if d in words[1]:
                    return data_type[d],words[1]
        self.util.warnMe('warning', ' no data type (extractDataType) - line:%s' % line)
        return None, words[1]

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
        return case, control, case+control
            

    ### write in Nextbio style batch info ###
    def writeNextBioStyle(self, filename,data):
        logging.debug(' Function:  writeNextBioStyle %s' % filename)
        header=["pmid","bioset_id", "data_type","broad_ethnicity","narrow_ethnicity","model_type","case_no","control_no","sample_size","is_meta","accepted","hapmap_pop","combined","comment","issue"]
        f1 = open("%s/%s_data4.csv" % (config.OUTPUTPATH,filename),'wb')
        f2 = open("%s/%s_info4.csv" % (config.OUTPUTPATH,filename),'wb')
        info_file = csv.writer(f2)
        info_file.writerow(header)
        data_file = csv.writer(f1)
        self.data_2_fix = {}

        data_file.writerow(["======================================================================================================================================================================================================================================="])
        for d in data:
            self.data_2_fix[d["BIOSET_ID"]] ={}
            for key in d:
                if key != 'data' and key !='column' :
                    self.data_2_fix[d["BIOSET_ID"]][key.lower()] = d[key]
                
            ### "pmid","bioset_id", 
            data_file.writerow(["BIOSET_ID" ,d["BIOSET_ID"]])
            data_file.writerow(["PMID" ,d["PMID"]])
            info_row = [d["PMID"],d["BIOSET_ID"]]                       # "pmid","bioset_id", 
            
            ### extract data_type dic, rep, comb
            if "BIOSET TITLE" in d:
                data_file.writerow(["BIOSET TITLE" ,d["BIOSET TITLE"]])
                data_type, title = self.extractDataType(d["BIOSET TITLE"].lower())
                info_row.append(data_type)                              #data_type 
            
            ### "broad_ethnicity", "narrow_ethnicity"
            info_row.extend([self.ethnicity_map[d['snp_population']],title])    #"broad_ethnicity" # "narrow_ethnicity"
            
            ### "model_type"
            if "BIOSET SUMMARY" in d:
                data_file.writerow(["BIOSET SUMMARY" ,d["BIOSET SUMMARY"]])
                model = self.extractModelType(d["BIOSET SUMMARY"].lower())
            if "ANALYSIS SUMMARY" in d:
                data_file.writerow(["ANALYSIS SUMMARY" ,d["ANALYSIS SUMMARY"]])
                model = self.extractModelType(d["ANALYSIS SUMMARY"].lower())
            info_row.append(model)                                     # "model_type"
            
            ### exctrat case_no control_no sample_size  
            if "SAMPLE NUMBER" in d:
                data_file.writerow(["SAMPLE NUMBER" ,d["SAMPLE NUMBER"]])
                info_row.extend(self.extractSampleSizes(d["SAMPLE NUMBER"].lower())) # "case_no","control_no","sample_size",
            
            if "COMPARISON" in d:
                data_file.writerow(["COMPARISON" ,d["COMPARISON"]])
            if "TAG" in d:
                data_file.writerow(["TAG" ,d["TAG"]])
            if "PLATFORM" in d:
                data_file.writerow(["PLATFORM" ,d["PLATFORM"]])
            
            ### extract is_meta yes or no                    
            if "META-ANALYSIS" in d:
                data_file.writerow(["META-ANALYSIS" ,d["META-ANALYSIS"]])
                info_row.append(self.extractMetaAnalysis(d["META-ANALYSIS"].lower()))#"is_meta",
            else:
                info_row.append("") 
            
            ### "accepted" "hapmap_pop"
            info_row.extend(["",d['snp_population']])                                        #"accepted" "hapmap_pop"
            
            ### "combined"
            if data_type != 'combined':
                info_row.append("no")              # "combined"
            else:
                info_row.append(title)
            ### add "comment","issue"]
            
            if "REFERENCE POPULATION" in d:
                data_file.writerow(["REFERENCE POPULATION" ,d["REFERENCE POPULATION"]])

            data_file.writerow(["BATCH" ,d["BATCH"]])
            data_file.writerow(["======================================================================================================================================================================================================================================="])
            
            for line in d["data"]:
                word = line.split("\t")
                data_file.writerow(word)
            data_file.writerow([""])
            info_file.writerow(info_row)
        f2.close()
        f1.close()

    ### recordFinishedBioset ###
    def recordFinishedBioset(self,data ):
        logging.debug(' Function:  recordFinishedBioset ')
        self.previous_batch = self.util.current_bioset
        self.util.warnMe('info', 'batch %s' % self.util.current_bioset['BIOSET_ID'] )
        self.util.current_bioset["data"] = data #TOADD
        self.util.current_bioset["snp_issue"] = self.snp_issue #TOADD
        if "FOUND_DISEASE" in self.util.current_bioset:
            self.disease_data[self.util.current_bioset["FOUND_DISEASE"]["disease"]]["column"][self.util.bioset_id]=(self.column_list_cap)
            self.disease_data[self.util.current_bioset["FOUND_DISEASE"]["disease"]]["studies"].append(self.util.current_bioset)


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
        
        ### keys expected in nextbio batch
        self.key_map = self.util.getTermMap(config.KEYMAP)

        ### get terms for and rejecting diseases of interest
        self.disease_map = self.util.getTermMap(config.DISEASEMAP)
        self.disease_data={}
        for d in  self.disease_map:
            self.disease_data[d] = {"studies":[],"column":{}}
     
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  GetData4HandCurration main' )

        ### get maps
        self.getMaps()
                      
        ### open input file
        self.util.batchname = "test2_CB3.txt"
        self.util.warnMe('info', 'Function:  GetDisease main - batch %s' % self.util.batchname ) 
        input_data = io.open("%s/%s" % (config.DATAPATH,self.util.batchname),'r',encoding='utf-8-sig')

        ### read input file
        line = input_data.readline() ### first line ========
        data =[]
        while len(line)>0:
            ### Starting new bioset                    
            while self.partition in line and len(line)>0:
                self.util.addNewKeyCombination()
                if self.util.issue == True:
                    self.issue_list.append({"current":self.util.current_bioset, "previous":self.previous_batch})
                self.util.issue = False
                if len(self.util.current_bioset)>0:
                    self.recordFinishedBioset(data )
#                     self.previous_batch = self.util.current_bioset
# #                     self.util.warnMe('info', 'batch %s' % self.util.current_bioset['BIOSET_ID'] )
#                     self.util.current_bioset["data"] = data #TOADD
#                     self.util.current_bioset["snp_issue"] = self.snp_issue #TOADD
#                     if "FOUND_DISEASE" in self.util.current_bioset:
#                         self.disease_data[self.util.current_bioset["FOUND_DISEASE"]["disease"]]["column"][self.util.bioset_id]=(self.column_list_cap)
#                         self.disease_data[self.util.current_bioset["FOUND_DISEASE"]["disease"]]["studies"].append(self.util.current_bioset)
                self.util.current_bioset ={}
                
                line = input_data.readline() ### bioset title/ID
                logging.debug(' Function:  GetData4HandCurration main - title %s' % line )
                while self.partition not in line and len(line)>0 :
                    key = self.util.getKeyValue(self.key_map,line)
                    if key is not None and "BIOSET TITLE" in key:
                        self.findDisease(line, key)
                    line = input_data.readline()
                ### starting SNP table
                if self.partition in line:
                    data =[]
                    self.snp_issue = {}
                    line = input_data.readline() ### header ot table
                    data.append(line)
                    self.column_list = line.lower().strip().split('\t')
                    self.column_list_cap = line.strip().split('\t')
                line = input_data.readline() ### delimiter before SNP data
                data.append(line)
                if self.snp_part in line:
                    line = input_data.readline() ### 1st SNP info
                    data.append(line)
                    self.util.getPMID(line)
            if "FOUND_DISEASE" in self.util.current_bioset:
                words = line.strip().split('\t')
                snp = words[2]
                snp_data = {}
                i = self.column_list.index('snp_population')
                ### check that snp table is correctly aline
                if len(words) != len(self.column_list) and len(words[i]) !=3 and words[i].lower() != 'Global or other'.lower():
                    self.util.warnMe('warning', 'different lenght %s %s' % (words, self.column_list))
                    line = line.replace('\t\t',"\t")
                    words = line.strip().split('\t')
                for nc in range(0,len(words)):
                    if words[nc] == '':
                        words[nc] = None
                    snp_data[self.column_list[nc].lower()] =  words[nc]
                ### record hapmap pop  
                if 'snp_population' not in self.util.current_bioset:
                    self.util.current_bioset['snps'] = {}
                    self.util.current_bioset['snp_population']= snp_data['snp_population']
                ### add new snp data in this bioset
                if snp not in self.util.current_bioset['snps']:
                    self.util.current_bioset['snps'][snp]=snp_data
                else:
                    self.util.warnMe('warning', 'found SNP: %s twice in one bioset: %s' %(snp, self.util.bioset_id))
                ### record if bioset has an issue
                if 'POSSIBLE_INCONSISTENCY'.lower() in snp_data:
                    if  snp_data['POSSIBLE_INCONSISTENCY'.lower()] != "":
                        self.snp_issue[snp]={'minor_allele':snp_data['minor_allele'],'major_allele':snp_data['major_allele'],'issue':['major_allele', 'minor_allele','_odds_ratio', '95%ci','_frequency_case','_frequency_control']}
                        if 'issue' not in self.util.current_bioset:
                            self.util.current_bioset['issue'] = snp_data['POSSIBLE_INCONSISTENCY'.lower()] 
                        elif self.util.current_bioset['issue'] != snp_data['POSSIBLE_INCONSISTENCY'.lower()] :
                            self.util.warnMe('warning', 'Different issue for SNP: %s twice in bioset: %s - %s - %s' %(snp, self.util.bioset_id, self.util.current_bioset['issue'], snp_data['POSSIBLE_INCONSISTENCY'.lower()] ))
                       
                ######### find New SNPs to add to DBSBP collection
#                 if self.util.bioset_id == words[0]:
#                     self.getSnpPList.fillSnpList(snp_data)
                ###########
            line = input_data.readline()  
            data.append(line)
        ### dont forget to record the last bioset 
        if len(line) ==0:
            self.recordFinishedBioset(data )
#             if len(self.util.current_bioset)>0:
#                 self.previous_batch = self.util.current_bioset
# #                     self.util.warnMe('info', 'batch %s' % self.util.current_bioset['BIOSET_ID'] )
#                 self.util.current_bioset["data"] = data 
#                 if "FOUND_DISEASE" in self.util.current_bioset:
#                     self.disease_data[self.util.current_bioset["FOUND_DISEASE"]["disease"]]["column"][self.util.bioset_id] = (self.column_list_cap)
#                     self.disease_data[self.util.current_bioset["FOUND_DISEASE"]["disease"]]["studies"].append(self.util.current_bioset)

            
        ######### find New SNPs to add to DBSBP collection
        ### output new SNP data
#         self.util.writeOutput("snp_List",self.getSnpPList.snp_list)
#         self.getSnpPList.writeListSnps()
        ###########
        
        ### output disease data
        for d in self.disease_data:
            nd =0 
            self.disease_data[d]['studies'].sort(key=operator.itemgetter('PMID'))
            for data in self.disease_data[d]['studies']:
                nd+=1
            self.util.warnMe('info', '%s, %s' % (d, nd))
            if nd >=1:
                self.fix_data = FixHandCurratedData(d)  
                self.fix_data.fixDataForCurration(self.data_2_fix,self.disease_data[d]['column'] )
                self.writeNextBioStyle('%s' % d,self.disease_data[d]['studies'])
                self.util.writeOutput('%s' % d,self.disease_data[d]['studies'])
             

        

""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'NextBio_get_diseases'), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = GetData4HandCurration()
    clean.main()
      

    print'DONE with DMmain'
    
""" END OF main """ 