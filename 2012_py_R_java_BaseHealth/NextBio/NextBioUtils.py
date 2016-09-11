#!/usr/bin/env python

"""
     Utility functions to obtain Nextbio format batch data
     06/12/12- 1.0
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
    6) change checkQuote to change " to ' inside string
 
"""

import logging,  json
import NextBioConfig as config

class NextBioUtils():
    def __init__(self):
        self.batchname = ""             # batch file provided by nextbio
        self.previous_batch ={}
        self.issue = False
        self.current_bioset ={}
        self.key_list = {}
        self.key_combination = []
        self.column_list = []
        self.key_comb_list = None
        self.column_name_dict = {}
        self.pmid_list = {}
        self.meta_sp   = [":", "="]     # this separates fields whithin meta data
        self.vs        = ["_vs_", "vs"] # this separates comparison and tag fields
        self.bioset_id = None


################# general function ################
    ''' keyToUpperCase '''
    def keyToUpperCase(self, data):
        logging.debug(' Function:  keyToUpperCase ')
        new_data = {}
        for key in data:
            new_key = key.upper()
            new_data[new_key] = data[key]
            if type(new_data[new_key]) == dict:
                new_data[new_key] = self.keyToUpperCase(new_data[new_key] )
        return new_data
    
    ''' checkSame
    compare two value'''
    def checkSame(self, data1, data2):
        logging.debug(' Function:  checkSame - data1: %s, data2: %s' %(data1, data2))
        if data1 != data2 :
            if data1 != '' and  data2 != '' and  data1 is not None and  data2 is not None  :
                logging.warning(' DATA DONT FIT (checkSame) - data1: %s, data2: %s' %(data1, data2))
                return False
        return True
    
    ''' get map for typo corrections '''
    def getTermMap(self,fname):
        logging.debug(' Function:  getTermMap' )
        f = open(fname,'r')
        typos = f.read()  
        f.close()
        return json.loads(typos)

######################### OUT FUNCTIONS ##############
    def writeOutput(self, filename,data):
        logging.debug(' Function:  writeOutput %s:' % filename)
        f = open("%s/%s.json" % (config.OUTPUTPATH,filename), "w")
        f.write(json.dumps(data, indent=4))
        f.close()
        
    ''' write unique key of list '''
    def writeList(self, filename,data):
        logging.debug(' Function:  writeList %s:' % filename)
        total = 0
        f = open("%s%s.json" % (config.OUTPUTPATH,filename), "w")
        for d in data:
            f.write("%s\t%s\n" % (d, data[d]["count"])) 
            total +=  data[d]["count"]
        f.write("TOTAL\t%s\n" % (total))

    ''' warnMe '''   
    def warnMe(self, flag, comment):
        if flag == 'info':
            logging.info(comment)
            print comment
        elif flag== 'warning':
            logging.warning(comment)
            print "%s -%s" % (flag.upper(),comment)
        elif flag== 'critical':
            logging.critical(comment)
            print "%s -%s" % (flag.upper(),comment)
        elif flag== 'error':
            logging.error(comment)
            print "%s -%s" % (flag.upper(),comment)
        else:
            logging.info(comment)
            print "%s -%s" % (flag.upper(),comment)     



################ ORIGINAL BATCH CLEAN UP ####################      
    ''' add new item to list '''
    def addToList(self,item_list, item, new_value, old_value):
        logging.debug(' Function:  addToList %s' % item )
        if item not in item_list:
            item_list[item] = new_value
        if old_value is not None:
            item_list[item][old_value] += 1   
        else:
            item_list[item] += 1   
        return item_list
    
    ############### check functions ###############
    ''' check key and value if duplicated and came from previous bioset '''
    def  checkKeyValue(self,key,value):
        logging.debug(' Function:  checkKeyValue, %s: %s' % (key,value) )
        if key not in self.key_combination:
                self.key_combination.append(key)
        else:
            self.issue = True
            logging.warning("REPEATED KEY IN ONE BIOSET (checkKeyValue) %s\t%s\t%s:%s" % (key, self.batchname, key,value))
            if self.current_bioset[key]["count"] >1 and self.current_bioset[key]["value"][value] == self.current_bioset[key]["count"]:
                logging.warning("VALUE FOUND TWICE IN ONE BIOSET (checkKeyValue): %s" %  value)
                if value in self.previous_batch[key]["value"]:
                    logging.warning("VALUE FOUND IN PREVIOUS BIOSET (checkKeyValue): %s" %  value)
                else:
                    logging.error("VALUE FOUND NOT IN PREVIOUS BIOSET FOR SAME KEY (checkKeyValue): %s:%s\t%s" %  (key, value, self.previous_batch[key]["value"]))
            else:
                self.warnMe("error", "DIFFERENT VALUE FOUND FOR SAME KEY (checkKeyValue): %s: %s\n%s" %  (key,value,self.current_bioset[key]["value"]))

    ''' remove quote from key or value 
        Issue from batch 7 
    '''
    def checkQuote(self,line):
        logging.debug(' Function:  checkQuote' )
        s=line.encode("utf-8")
        if "\"" in s: 
            s = line.split("\"")
            if len(s) == 2:
                if len(s[1]) == 0 :
                    logging.warning(' REMOVE QUOTE s[1]=0 (checkQuote) - %s - %s\n%s' % (self.batchname,line, s))
                    s =   s[0]
                elif len(s[0])== 0 :
                    logging.warning(' REMOVE QUOTE s[0]=0 (checkQuote) - %s - %s\n%s' % (self.batchname,line, s))
                    s =   s[1]  
                else:
                    self.warnMe('error', ' ERROR QUOTE (checkQuote)  - %s - %s\n%s' % (self.batchname, line,s))
                    s= line
            else:
                s1 = None
                for w in s:
                    if s1 is None:  s1=w
                    else:  s1 = "%s'%s" % (s1,w)
                s = s1
                logging.warning(' CHANGED QUOTE (checkQuote) - %s - %s\n%s' % (self.batchname, line,s))
        return s
    ############### 
    
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
                    self.warnMe('critical',' NO KEY VALUE? (findKey) %s' % line )
        return key, value
    
    ''' get key and value from a line in meta data and count # occurence in batch''' 
    def getKeyValueCount(self,key_map,  line):
        logging.debug(' Function:  getKeyValue' )
        line = line.strip()
        key, value = self.findKey(line)
        ''' record list of keys / value'''
        self.current_bioset["BATCH"]=self.batchname
        if key :                   
            key = key.strip().upper()
            key = self.checkQuote(key)            
            if key in key_map:
                key = key_map[key]
            else:
                self.warnMe('critical', "I COULD NOT FIND MAPPING FOR KEY %s (getKeyValueCount)- %s" % (key, self.current_bioset))
            ''' add new key '''
            self.key_list = self.addToList(self.key_list, key, {"count" : 0, "value" : {}, "batch":{}}, "count")
            ''' record value'''
            value = value.strip()
            value = self.checkQuote(value)
            self.key_list[key]["value"] = self.addToList(self.key_list[key]["value"], value, {"count":0,"batch":{}}, "count")
            self.key_list[key]["value"][value]["batch"] = self.addToList(self.key_list[key]["value"][value]["batch"], self.batchname, 0, None)  
            self.current_bioset[key]["value"] = self.addToList(self.current_bioset[key]["value"], value, 0, None)    
            self.current_bioset[key]=value
            ''' record batchname '''
            self.key_list[key]["batch"] = self.addToList(self.key_list[key]["batch"], self.batchname, 0, None)
                        
            self.checkKeyValue(key, value)
        return key

    ''' get key and value from a line in meta data ''' 
    def getKeyValueLower(self, key_map, line):
        logging.debug(' Function:  getKeyValue' )
        line = line.strip()
        key, value = self.findKey(line)
        ''' record list of keys / value'''
        self.current_bioset["BATCH".lower()]=self.batchname
        if key :                   
            key = key.strip().upper()
            key = self.checkQuote(key)
            if key in key_map:
                key = key_map[key].lower()
            else:
                self.warnMe('critical', "I COULD NOT FIND MAPPING FOR KEY %s (getKeyValue)- %s" % (key, self.current_bioset))
            ''' add new key '''
            self.key_list = self.addToList(self.key_list, key, {"count" : 0, "value" : {}, "batch":{}}, "count")
            ''' record value'''
            value = value.strip()
            value = self.checkQuote(value)
            self.key_list[key]["value"] = self.addToList(self.key_list[key]["value"], value, 0, None)
            self.current_bioset[key]=value
            ''' record batchname '''
            self.key_list[key]["batch"] = self.addToList(self.key_list[key]["batch"], self.batchname, 0, None)
                        
            self.checkKeyValue(key, value)
        return key
               
    ''' get key and value from a line in meta data ''' 
    def getKeyValue(self, key_map, line):
        logging.debug(' Function:  getKeyValue' )
        line = line.strip()
        key, value = self.findKey(line)
        ''' record list of keys / value'''
        self.current_bioset["BATCH"]=self.batchname
        if key :                   
            key = key.strip().upper()
            key = self.checkQuote(key)
            if key in key_map:
                key = key_map[key]
            else:
                self.warnMe('critical', "I COULD NOT FIND MAPPING FOR KEY %s (getKeyValue)- %s" % (key, self.current_bioset))
            ''' add new key '''
            self.key_list = self.addToList(self.key_list, key, {"count" : 0, "value" : {}, "batch":{}}, "count")
            ''' record value'''
            value = value.strip()
            value = self.checkQuote(value)
            self.key_list[key]["value"] = self.addToList(self.key_list[key]["value"], value, 0, None)
            self.current_bioset[key]=value
            ''' record batchname '''
            self.key_list[key]["batch"] = self.addToList(self.key_list[key]["batch"], self.batchname, 0, None)
                        
            self.checkKeyValue(key, value)
        return key
           
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
            
    ''' get column names for SNP data
    ''' 
    def getColumnList(self, line):         
        logging.debug(' Function:  getColumnNames' )
        line = line.strip()
        self.column_list = []
        for word in line.split('\t'):
            if word not in self.column_list:
                self.column_list.append(word)
            else:
                word = '%s2' % word
                self.column_list.append(word)
                self.warnMe('warning', ' WORD REPEATED %s (getColumnList) %s' % (word, line))
            
    ''' get unique PMID and they bioset IDs
    ''' 
    def getPMID(self, line):         
        logging.debug(' Function:  getPMID' )
        line = line.strip()
        data = line.split('\t')
        pmid = data[1]
#         self.current_bioset["PMID"] = pmid
#         self.current_bioset["BIOSET_ID"] = data[0]
        bioset_id = data[0]
        if pmid not in self.pmid_list:
            self.pmid_list[pmid] = {"count":0,"id":{}}
        self.pmid_list[pmid]["count"] += 1
        if bioset_id not in self.pmid_list[pmid]["id"]:
            self.pmid_list[pmid]["id"][bioset_id] = 0
        self.pmid_list[pmid]["id"][bioset_id] +=1   
        self.bioset_id =  bioset_id
        self.current_bioset["BIOSET_ID".lower()]=bioset_id
        self.current_bioset["PMID".lower()]=pmid    

    ''' add new key combination in list ''' 
    def addNewKeyCombination(self):
        logging.debug(' Function:  addNewKeyCombination %s' %(self.key_combination) )
        if len(self.key_combination)>0:
            ''' init combination list '''            
            if self.key_comb_list is None:
                self.key_comb_list = {}
                self.key_comb_list = self.addToList(self.key_comb_list, "0", {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
                self.key_comb_list["0"]["batch"] = self.addToList(self.key_comb_list["0"]["batch"], self.batchname, 0, None)
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

