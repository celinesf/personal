#!/usr/bin/env python

"""
    Clean up of NextBio raw batch data files
   05/09/13 - 1.1
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
"""

import logging, os, json
import nb_conf as config

class RawDataCleanUp():
    def __init__(self):
        self.batchname = ""             # batch file provided by nextbio
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.previous_batch ={}
        self.issue = False
        self.issue_list = []
        self.current_batch ={}
        self.key_list = {}
        self.key_map = None
        self.key_combination = []
        self.key_comb_list = None
        self.column_name_dict = {}
        self.pmid_list = {}
        self.meta_sp   = [":", "="]     # this separates fields whithin meta data
        self.vs        = ["_vs_", "vs"] # this separates comparison and tag fields
      
    ''' get map for typo corrections '''
    def getTermMap(self,fname):
        f = open(fname,'r')
        typos = f.read()  
        return json.loads(typos)
        
      
    ''' remove quote from key or valye 
        Issue from batch 7 
    '''
    def checkQuote(self,line):
        logging.debug(' Function:  checkQuote' )
        s = line
        if "\"" in line: 
            s = line.split("\"")
            if len(s[0]) > 0:
                return s[0]
            else:
                return s[1]
        else:
            return line
    
    '''find key before : or ='''
    def findKey(self,line):
        logging.debug(' Function:  findKey' )
        key = None
        value = None
        for sep in self.meta_sp: 
            if sep in line:
                try:
                    key, value = line.split(sep,1)
                except Exception, e:
                    print line
                    print "Exception : ", str(e)
                    logging.CRITICAL(' Function:  findKey %s' % line )
        return key, value
    
    ''' add new item to list '''
    def addToList(self,item_list, item, new_value, old_value):
        if item not in item_list:
            item_list[item] = new_value
        if old_value is not None:
            item_list[item][old_value] += 1   
        else:
            item_list[item] += 1   
        return item_list
               
    ''' get key and value from a line in meta data ''' 
    def getKeyValue(self, line):
        logging.debug(' Function:  getKeyValue' )
        
        line = line.strip()
        key, value = self.findKey(line)

        ''' record list of keys / value'''
        if key :                   
            key = key.strip().upper()
            key = self.checkQuote(key)
            
            if key in self.key_map:
                key = self.key_map[key]
            else:
                logging.critical("ERROR I COULD NOT FIND MAPPING FOR KEY %s" %key)
                print "ERROR I COULD NOT FIND MAPPING FOR KEY %s" %key
            
            ''' add new key '''
            self.key_list = self.addToList(self.key_list, key, {"count" : 0, "value" : {}, "batch":{}}, "count")
            self.current_batch = self.addToList(self.current_batch, key,{"count" : 0, "value" : {}, "batch":{}}, "count")
            
            ''' record value'''
            value = value.strip()
            value = self.checkQuote(value)
            self.key_list[key]["value"] = self.addToList(self.key_list[key]["value"], value, 0, None)  
            self.current_batch[key]["value"] = self.addToList(self.current_batch[key]["value"], value, 0, None)    
            
            ''' record batchname '''
            self.key_list[key]["batch"] = self.addToList(self.key_list[key]["batch"], self.batchname, 0, None)
            self.current_batch[key]["batch"] = self.addToList(self.current_batch[key]["batch"], self.batchname, 0, None)
            
            self.checkKeyValue(key, value)
            
            
    ''' check key and value if duplicated and came from previous bioset '''
    def  checkKeyValue(self,key,value):
        logging.debug(' Function:  checkKeyValue, %s: %s' % (key,value) )
        if key not in self.key_combination:
                self.key_combination.append(key)
        else:
            self.issue = True
            logging.warning("REPEATED KEY IN ONE BIOSET %s\t%s\t%s:%s" % (key, self.batchname, key,value))
            if self.current_batch[key]["count"] >1 and self.current_batch[key]["value"][value] == self.current_batch[key]["count"]:
                logging.warning("VALUE FOUND TWICE IN ONE BIOSET: %s" %  value)
                if value in self.previous_batch[key]["value"]:
                    logging.warning("VALUE FOUND IN PREVIOUS BIOSET: %s" %  value)
                else:
                    logging.error("VALUE FOUND NOT IN PREVIOUS BIOSET FOR SAME KEY: %s:%s\t%s" %  (key, value, self.previous_batch[key]["value"]))
            else:
                logging.error("DIFFERENT VALUE FOUND FOR SAME KEY: %s: %s\n%s" %  (key,value,self.current_batch[key]["value"]))

#                 print "I found 2 keys the same %s, %s, %s" % (key, self.batchname, line)
#                 print self.key_combination
                 
    ''' get column names for SNP data
    ''' 
    def getColumnNames(self, line):         
        logging.debug(' Function:  getColumnNames' )
        line = line.strip()
        for word in line.split('\t'):
            if word not in self.column_name_dict:
                self.column_name_dict[word] = {"count" : 0, "batch":{}}
            self.column_name_dict[word]["count"] += 1
            if self.batchname not in self.column_name_dict[word]["batch"]:
                self.column_name_dict[word]["batch"][self.batchname] = 0
            self.column_name_dict[word]["batch"][self.batchname] += 1
            
    ''' get unique PMID and they bioset IDs
    ''' 
    def getPMID(self, line):         
        logging.debug(' Function:  getPMID' )
        line = line.strip()
        data = line.split('\t')
        pmid = data[1]
        bioset_id = data[0]
        if pmid not in self.pmid_list:
            self.pmid_list[pmid] = {"count":0,"id":{}}
        self.pmid_list[pmid]["count"] += 1
        if bioset_id not in self.pmid_list[pmid]["id"]:
            self.pmid_list[pmid]["id"][bioset_id] = 0
        self.pmid_list[pmid]["id"][bioset_id] +=1    
            
    def writeOutput(self, filename,data):
        f = open("%s%s.json" % (config.OUTPUTPATH,filename), "w")
        f.write(json.dumps(data, indent=4))
        
    ''' write unique key of list '''
    def writeList(self, filename,data):
        total = 0
        f = open("%s%s.json" % (config.OUTPUTPATH,filename), "w")
        for d in data:
            f.write("%s\t%s\n" % (d, data[d]["count"])) 
            total +=  data[d]["count"]
        f.write("TOTAL\t%s\n" % (total))
    
    ''' add new key combination in list ''' 
    def addNewKeyCombination(self):
        logging.debug(' Function:  addNewKeyCombination %s' %(self.key_combination) )
        if len(self.key_combination)>0:
            ''' init combination list '''            
            if self.key_comb_list is None:
                self.key_comb_list = {}
                self.key_comb_list = self.addToList(self.key_comb_list, "0", {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
                self.key_comb_list["0"]["batch"] = self.addToList(self.key_comb_list["0"]["batch"], self.batchname, 0, None)
#                 self.key_comb_list["0"] = {"keys":self.key_combination,"count":0,"batch":{}}
            else:
                found = 0
                comb_num =""
                for comb in self.key_comb_list:
                    comb_num = comb
                    it = 0
                    for key in self.key_comb_list[comb]["keys"]:
                        if key in self.key_combination:
                            it += 1
                        else:
                            break
                    if it == len(self.key_combination) and it == len(self.key_comb_list[comb]['keys']):
                        found = 1
                        self.key_comb_list = self.addToList(self.key_comb_list, comb, {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
                        self.key_comb_list[comb]["batch"] = self.addToList(self.key_comb_list[comb]["batch"], self.batchname, 0, None)
                        break
                if found == 0:
                    comb_num = str(len(self.key_comb_list))
                    self.key_comb_list = self.addToList(self.key_comb_list, comb_num, {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
                    self.key_comb_list[comb_num]["batch"] = self.addToList(self.key_comb_list[comb_num]["batch"], self.batchname, 0, None)
        self.key_combination = []
     
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  RawDataCleanUp main' )
        ''' read batch files'''
        batches_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if i.startswith("GWAS") and i.endswith('.txt')] 
        self.key_map = self.getTermMap(config.KEYMAP)

        for batch in batches_list:
            self.batchname = batch
            print self.batchname
#             if self.batchname == 'GWAS_standard1.txt':
            input_data = open("%s%s" % (config.DATAPATH,self.batchname),'r')

            line = input_data.readline() ### first line ========
            while len(line)>0:
                ### Starting new bioset                    
                while self.partition in line and len(line)>0:
                    self.addNewKeyCombination()
                    if self.issue == True:
                        self.issue_list.append({"current":self.current_batch, "previous":self.previous_batch})
                    self.issue = False
                    if len(self.current_batch)>0:
                        self.previous_batch = self.current_batch
                    self.current_batch ={}
                    
                    line = input_data.readline() ### bioset title/ID
                    logging.debug(' Function:  RawDataCleanUp main - title %s' % line )
                    while self.partition not in line and len(line)>0 :
                        self.getKeyValue(line)
                        line = input_data.readline()
                    ### starting SNP table
                    if self.partition in line:
                        line = input_data.readline() ### header ot table
                        self.getColumnNames(line)
                    line = input_data.readline() ### delimiter before SNP data
                    if self.snp_part in line:
                        line = input_data.readline() ### 1st SNP info
                        self.getPMID(line)
                line = input_data.readline()  

        self.writeOutput('Issues',self.issue_list)
        self.writeOutput('key',self.key_list)
        self.writeOutput('key_comb',self.key_comb_list)
#         self.writeOutput('column',self.column_name_dict)
#         self.writeOutput('pmid',self.pmid_list)
        
        self.writeList('key_list',self.key_list)
#         self.writeList('column_name_dict',self.column_name_dict)
#         self.writeList('pmid_list',self.pmid_list)



""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'NextBio_Clean_up'), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = RawDataCleanUp()
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 