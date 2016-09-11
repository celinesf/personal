#!/usr/bin/env python

"""
    Clean up of NextBio raw batch data files
   05/06/13 - 1.0 
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Recover all keys/tags from batches 
    2) Recover all column names for SNP information tables from all batchs
"""

import logging, os, json
import nb_conf as config

class RawDataCleanUp():
    def __init__(self):
        self.batchname = ""             # batch file provided by nextbio
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.key_list = {}
        self.column_name_dict = {}
        self.pmid_list = {}
        self.meta_sp   = [":", "="]     # this separates fields whithin meta data
        self.vs        = ["_vs_", "vs"] # this separates comparison and tag fields
                  
    ''' get key and value from a line in meta data
    ''' 
    def getKeyValue(self, line):
        logging.debug(' Function:  getKeyValue' )
        key = None
        value = None
        line = line.strip()
        for sep in self.meta_sp:
            if sep in line:
                try:
                    key, value = line.split(sep,1)
                except Exception, e:
                    print line
                    print "Exception : ", str(e)
        ### record list of keys
        if key :                   
            key = key.strip()
            value = value.strip()
            if key not in self.key_list:
#                 self.key_list[key] = {"count" : 0, "value" : {}, "batch":{}}
                self.key_list[key] = {"count" : 0,  "batch":{}}
            self.key_list[key]["count"]+= 1
            
#             if value not in self.key_list[key]["value"]:
#                 self.key_list[key]["value"][value] = 0
#             self.key_list[key]["value"][value] += 1     
            
            if self.batchname not in self.key_list[key]["batch"]:
                self.key_list[key]["batch"][self.batchname] = 0
            self.key_list[key]["batch"][self.batchname] += 1     
                 
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
        
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  RawDataCleanUp main' )
        ''' read batch files'''
        batches_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if i.startswith("GWAS") and i.endswith('.txt')] 

        for batch in batches_list:
            self.batchname = batch
            print self.batchname
#             if self.batchname == 'GWAS_standard1.txt':
            input_data = open("%s%s" % (config.DATAPATH,self.batchname),'r')

            line = input_data.readline() ### first line ========
            while len(line)>0:
                ### Starting new bioset
                while self.partition in line and len(line)>0:
                    line = input_data.readline() ### bioset title/ID
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

        self.writeOutput('key_',self.key_list)
        self.writeOutput('column',self.column_name_dict)
        self.writeOutput('pmid',self.pmid_list)
        
        self.writeList('key_list',self.key_list)
        self.writeList('column_name_dict',self.column_name_dict)
        self.writeList('pmid_list',self.pmid_list)



""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'NextBio_Clean_up'), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = RawDataCleanUp()
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 