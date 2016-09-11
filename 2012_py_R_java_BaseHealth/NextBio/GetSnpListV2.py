#!/usr/bin/env python

"""
     Get SNP List from all the batches
     08/15/13- 2.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Recover all the SNPs in NextBio batches in order to get info from dbsnp 
    2) read data from csv files output by GetData4HandCurrationV2
"""

import logging, os, pymongo, csv
import NextBioConfig as config
from NextBioUtils import NextBioUtils

conn = pymongo.Connection(config.MONGO, 27017)
nextbio_db = conn["nextbio"]

class GetSnpList():
    def __init__(self,filename):
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.snp_list = {}
        self.filename = filename
        self.util = NextBioUtils()

    ''' csvDataFiles '''
    def csvDataFiles(self):
        logging.debug(' Function:  csvDataFiles' )
        ### get batch file names
 
        ### get old list of snp
        db_data = nextbio_db['dbsnp'].find({},{'dbsnp':1,'_id':0}) # find all
        old_list =[]
        for d in db_data :
            old_list.append(d['dbsnp'])

        ### output file of snps
        new_list = []
        ### open input file
        disease = self.filename.split('_data')[0]
        nfile = 0
        output = open("%s/snp_list2_%s_%s.txt"% (config.OUTPUTPATH,disease,nfile),'w')
        csvfile = open('%s/%s.csv' % (config.DATAPATH,self.filename), 'rb') 
        spamreader = csv.reader(csvfile)
        ### read file
        nsnp = 0
        nrow = 0
        for row in spamreader:
            if len(row)>2:
                nrow +=1
                if row[2] not in old_list and row[2] != 'SNP' and row[2] not in new_list:
                    new_list.append(row[2])
                    output.write("%s\n" % row[2] )
                    nsnp += 1
                    if nsnp % 5000 == 0:
                        output.close()
                        nfile += 1
                        output = open("%s/snp_list_%s_%s.txt"% (config.OUTPUTPATH,disease,nfile),'w')
            if nrow % 5000 == 0 and nrow !=0:
                self.util.warnMe('warning', 'nrow %s' % nrow)
        output.close() 



""" MAIN function """    
if __name__ == "__main__":
    
    filename = "attention_disorder_data_v8"
    logging.basicConfig(filename='%s/%s_%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,filename,'_GetSnpList'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )
  
    clean = GetSnpList(filename)
    clean.csvDataFiles()
    
    print'DONE with %s' % filename
    
""" END OF main """ 