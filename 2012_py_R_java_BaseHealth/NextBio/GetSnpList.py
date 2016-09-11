#!/usr/bin/env python

"""
     Get SNP List from all the batches
     06/12/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Recover all the SNPs in NextBio batches in order to get info from dbsnp 
"""

import logging, os, sys,pymongo,  io, copy, xlrd, csv
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from ExcelUtils import ExcelUtils

conn = pymongo.Connection(config.MONGO, 27017)
nextbio_db = conn["nextbio"]

class GetSnpList():
    def __init__(self):
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section
        self.snp_list = {}
        self.util = NextBioUtils()
        self.in_util = ExcelUtils()
    
    ### writeListSnps ###  
    def writeListSnps(self):
        logging.debug(' Function:  writeListSnps' )
        nfile = 0
        nsnp = 0
        output = open("%s/snp_list_%s.txt"% (config.OUTPUTPATH,nfile),'w')
        output_count = open("%s/snp_list_count.txt"% (config.OUTPUTPATH),'w')
        db_data = nextbio_db['dbsnp'].find({},{'dbsnp':1,'_id':0}) # find all
        self.util.warnMe('info', ' total snps - %s' %len(self.snp_list))
        for d in db_data:
            if d['dbsnp'] in self.snp_list:
                logging.info( 'found doc for snp %s' % d['dbsnp'] ) 
                self.snp_list.pop(d['dbsnp'])
        self.util.warnMe('info', ' new snps - %s' %len(self.snp_list))
        for snp in self.snp_list:
#             db_data = nextbio_db['dbsnp'].find({'dbsnp':snp})
            doc = None
#             for doc in db_data:
#                 logging.info( 'found doc for snp %s' % snp) 
            if doc is None:
                logging.info(' NEW  snp %s' % snp)
                output.write("%s\n" % snp)
                output_count.write("%s\t%s\n" % (snp,len(self.snp_list[snp]['id'])))
                nsnp += 1
                if nsnp % 5000 == 0:
                    output.close()
                    nfile += 1
                    output = open("%s/snp_list_%s.txt"% (config.OUTPUTPATH,nfile),'w')
        output.close() 
        output_count.close()
   
    ### get list of columns names to recover per table ### 
    def getSnpInfoNames(self):
        logging.debug(' Function:  getSnpInfoNames' )
        self.column_map =[]
        for key in self.snp_info_template:
            self.column_map.append(key)
            if type(self.snp_info_template[key]) == dict : # SNP_POPULATION dictionary
                for k in self.snp_info_template[key]['<pop>']:
                    self.column_map.append(k)
 
    ### addNewPopulationFrequencies ###  
    def addNewPopulationFrequencies(self, record, new_data):
        logging.debug(' Function:  addNewPopulationFrequencies' )
        key = 'snp_population'
        pop = new_data[key]
        hap_data = copy.deepcopy(self.snp_info_template[key]['<pop>'])
        for k in hap_data:
#             print 'addNewPopulationFrequencies', k
            if k in new_data:
                hap_data[k] = new_data[k]
                new_data.pop(k)
        record[pop]=hap_data 
        return record, new_data
    
    ### addNewSnp ###  
    def addNewSnp(self, data):
        logging.debug(' Function:  addNewSnp' )
        new_snp = copy.deepcopy(self.snp_info_template)
        for key in new_snp:
            key2 = key.lower()
#             print 'addNewSnp',new_snp, key, data
            if type(new_snp[key]) == dict : # SNP_POPULATION dictionary
                new_snp[key]={} # re-init SNP_POPULATION
                new_snp[key], data = self.addNewPopulationFrequencies(new_snp[key], data)
            elif type(new_snp[key]) == list: #  list of bioset ID
                new_snp[key].append(data[key2])
            else:
                new_snp[key] = data[key2]
        p = 0
        for k in new_snp['snp_population']:
            self.pop = k
            p+=1
        if p==0 or p>1:
            self.util.warnMe('warning', 'here %s %s' %(p, new_snp['snp_population']))
        self.snp_list[new_snp['snp']] = copy.deepcopy(new_snp)

    ### checkGene ###  
    def checkGene(self, gene1, gene2):
        logging.debug(' Function:  checkGene : gene1:%s, gene2 %s' %( gene1, gene2) )
        flag = True
        gene1= gene1.split("|")
        gene2=gene2.split("|")
        for gene in gene1:
            if gene not in gene2:
                flag = False
                break
        for gene in gene2:
            if gene not in gene1:
                flag = False
                break
        return flag

     ### reverseGenotypeFrequencies ###  
    def reverseGenotypeFrequencies(self, key, data):
        logging.debug(' Function:  reverseGenotypeFrequencies - key:%s, data:%s' %(key, data ) )

        if 'G' in data['allele_1'] or 'G' in data['allele_2'] and 'G' not in data[key] and 'C' in data[key] :
            data[key] = data[key].replace('C',"G")
        elif 'C' in data['allele_1'] or 'C' in data['allele_2'] and 'C' not in data[key]and 'G' in data[key] :
            data[key] = data[key].replace('G',"C")
        if 'T' in data['allele_1'] or 'T' in data['allele_2'] and 'T' not in data[key]and 'A' in data[key] :
            data[key] = data[key].replace('A',"T")
        elif 'A' in data['allele_1'] or 'A' in data['allele_2'] and 'A' not in data[key]and 'T' in data[key] :
            data[key] = data[key].replace('T',"A")

    ### flipHeterozygous ###  
    def flipHeterozygous(self, key, data):    
        logging.debug(' Function:  flipHeterozygous - key:%s, data:%s' %(key, data ) )
        if 'C/T' in data[key]:
            data[key] = data[key].replace('C/T','T/C')
        elif 'T/C' in data[key]:
            data[key] = data[key].replace('T/C','C/T')
        elif 'A/T' in data[key]:
            data[key] = data[key].replace('A/T','T/A')
        elif 'T/A' in data[key]:
            data[key] = data[key].replace('T/A','A/T')
        elif 'A/C' in data[key]:
            data[key] = data[key].replace('A/C','C/A')
        elif 'C/A' in data[key]:
            data[key] = data[key].replace('C/A','A/C')
        elif 'A/G' in data[key]:
            data[key] = data[key].replace('A/G','G/A')
        elif 'G/A' in data[key]:
            data[key] = data[key].replace('G/A','A/G')
        elif 'C/G' in data[key]:
            data[key] = data[key].replace('C/G','G/C')
        elif 'G/C' in data[key]:
            data[key] = data[key].replace('G/C','C/G')
        elif 'T/G' in data[key]:
            data[key] = data[key].replace('T/G','G/T')
        elif 'G/T' in data[key]:
            data[key] = data[key].replace('G/T','T/G')

        ### checkGenotypingFrequencies ###  
    def checkGenotypingFrequencies(self, data1, data2, num):
        logging.debug(' Function:  checkGenotypingFrequencies -num:%s, data1:%s, data2:%s' %(num, data1['genotype_frequency'], data2['genotype_frequency'] ) )
        key = 'genotype_frequency'
        flag = False
        if num == 0:
            self.reverseGenotypeFrequencies(key, data1)
        elif num == 1:
            self.reverseGenotypeFrequencies(key, data2)
        elif num == 2:
            self.flipHeterozygous(key, data2)
            
        if data1[key] == data2[key]:
            flag = True
        elif data1[key].split('|')[0] != data2[key].split('|')[0]:
            flag = self.checkGenotypingFrequencies(data1, data2, 1)
        elif data1[key].split('|')[2] != data2[key].split('|')[2]:
            flag = self.checkGenotypingFrequencies(data1, data2, 2)
        
        return flag


    ### checkKeyValue ###  
    def checkKeyValue(self, key1, old_data, key2,new_data, flag_rev):
        logging.debug(' Function:  checkKeyValue flag: %s, key1: %s, old:%s, key2: %s, new:%s' %( flag_rev,key1,old_data[key1],key2,new_data[key2] ) )
        flag_fit = False
        if flag_rev == True:
            if "allele_1" in key2 :
                key2 = key2.replace('allele_1','allele_2')
            elif "allele_2" in key2 :
                key2 = key2.replace('allele_2','allele_1')
            elif "genotype_frequency" in key2 :
                ### flip allele1/allele2 freqs
                geno = new_data[key2].split('|')
                new_data[key2] = '%s|%s|%s|%s|%s|%s' % (geno[4],geno[5],geno[2],geno[3],geno[0],geno[1]) 
                
        if old_data[key1] != new_data[key2] and new_data[key2] !="" and old_data[key1] !="":
            if flag_rev == False:
                flag_fit, flag_rev = self.checkKeyValue(key1,old_data,key2, new_data, True)
            else:   
                if key1 ==  'genotype_frequency':
                    ### flip allele1/allele2 freqs
                    geno = new_data[key2].split('|')
                    new_data[key2] = '%s|%s|%s|%s|%s|%s' % (geno[4],geno[5],geno[2],geno[3],geno[0],geno[1])
                    flag_fit = self.checkGenotypingFrequencies(old_data,new_data, 0)
                    print flag_fit
                    flag_rev = False               
                elif key1 =="gene":
                    flag_fit = self.checkGene(old_data[key1],new_data[key1])
                    flag_rev = False
                else:
                    flag_fit = False
                    self.util.warnMe('critical', ' VALUES DONT FIT - (checkKeyValue) snp:%s, pop:%s, flag: %s, key1: %s, old:%s, key2: %s, new:%s' %(self.snp, self.pop,flag_rev,key1,old_data[key1],key2,new_data[key2] ) )
                    print old_data,new_data
        else:
            if "genotype_frequency" in key2 :
                flag_rev = False
            logging.debug(' VALUES FIT - (checkKeyValue) flag: %s, key1: %s, old:%s, key2: %s, new:%s' %(flag_rev,key1,old_data[key1],key2,new_data[key2] ) )
            flag_fit = True
        return flag_fit, flag_rev
    
    ### checkPopulationFrequenciesFit ###  
    def checkPopulationFrequenciesFit(self, old_data, new_data, flag_rev):
        logging.debug(' Function:  checkPopulationFrequenciesFit flag: %s, old:%s, new:%s' %(flag_rev,old_data,new_data ) )      
        flag_fit = False
        for pop in new_data:
            for key in old_data[pop]:
#                 print 'checkPopulationFrequenciesFit', key
                flag_fit, flag_rev = self.checkKeyValue(key,old_data[pop],key, new_data[pop], flag_rev)
                if flag_rev == False:
                    break
            if flag_rev == False:
                    break
        return flag_fit, flag_rev
  
    ### checkExistingSnp ###  
    def checkExistingSnp(self, data):
        logging.debug(' Function:  checkExistingSnp' )
        self.snp = data['snp']
        old_snp = copy.deepcopy(self.snp_list[data['snp']])
        
        flag_rev = False
        flag_fit = False
        for key in old_snp:
#             print 'checkExistingSnp', key
            if type(old_snp[key]) == dict : # SNP_POPULATION dictionary
                if data[key] in old_snp[key]: # existing hapmap pop
                    data[key],data = self.addNewPopulationFrequencies({}, data)
                    flag_fit , flag_rev = self.checkPopulationFrequenciesFit(old_snp[key],data[key], flag_rev )
                else:
                    old_snp[key],data = self.addNewPopulationFrequencies(old_snp[key], data)
            elif type(old_snp[key]) == list: #  list of bioset ID
                if data[key] not in old_snp[key]:
                    old_snp[key].append(data[key])
#                 else:
#                     self.util.warnMe("critical", ' BIOSET IS REPLICATED (checkExistingSnp) - snp:%s old: %s, new: %s' % (data['snp'],old_snp[key],data[key]))
            else:
                flag_fit, flag_rev = self.checkKeyValue(key,old_snp,key, data, flag_rev)
            if flag_fit == False:
                print "NO FIT",self.snp, key, flag_rev
                print "old", old_snp
                print 'new', data
                break
        self.snp_list[old_snp['snp']] = copy.deepcopy(old_snp)
            

    ### fill up SNP list ###  
    def fillSnpList(self, data):
        logging.debug(' Function:  fillSnpList %s' ,data['snp'])
        if data['snp'] not in self.snp_list:
            self.addNewSnp(data)
        else:
            self.checkExistingSnp(data)
    

 ### recover the data for snp in a bioset ###  
    def recoverSnpData(self):
        logging.debug(' Function:  recoverSnpData' )  
        for bioset in self.in_util.nextbio_data:
            for snp in self.in_util.nextbio_data[bioset]['snps']:
                self.fillSnpList(self.in_util.nextbio_data[bioset]['snps'][snp])


    ''' main4disease '''
    def main4disease(self):
        logging.debug(' Function:  main4disease' )
        ### get batch file names
        batches_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if i.endswith('.xls')] 
        
        ### keys expected in nextbio batch
        self.key_map = self.util.getTermMap(config.KEYMAP)
        ### get column name to record info from
        tmp =  self.util.getTermMap(config.SNPINFOMAP)
        self.snp_info_template = {}
        for key in tmp:
            if tmp[key] is not None and type(tmp[key]) is str:
                self.snp_info_template[key.lower()] = tmp[key].lower()
            elif type(tmp[key]) is dict:
                self.snp_info_template[key.lower()] = {'<pop>':{}}
                for k in tmp[key]['<pop>']:
                    self.snp_info_template[key.lower()]['<pop>'][k.lower()] = tmp[key]['<pop>'][k]
            else : self.snp_info_template [key.lower()] = tmp[key]

        self.getSnpInfoNames()
        
#         print batches_list
#         for batch in batches_list:
#             self.util.batchname = batch
#             print batch
#             ''' get data for this disease '''
#             self.in_util.readExcelData('%s/%s' % (config.DATAPATH, batch))
#             print self.in_util.nextbio_data
#             self.recoverSnpData()
#         # outputData   
#         self.util.writeOutput("snp_List",self.snp_list)
#         self.writeListSnps()
     
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  GetSnpList main' )
        ### get batch file names
        batches_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if i.startswith("GWAS") and i.endswith('.txt')] 
        
        ### keys expected in nextbio batch
        self.key_map = self.util.getTermMap(config.KEYMAP)
        ### get column name to record info from
        self.snp_info_template = self.util.getTermMap(config.SNPINFOMAP)
        self.getSnpInfoNames()
        
        print batches_list
        for batch in batches_list:
            self.util.batchname = batch
#         ''' '''
        self.util.batchname = "GWAS_standard1_DPD_CB3.txt"
#         ''' '''
        print self.util.batchname   
        logging.debug(' Function:  GetSnpList main - batch %s' % self.util.batchname ) 

        input_data = io.open("%s/%s" % (config.DATAPATH,self.util.batchname),'r',encoding='utf-8-sig')

        line = input_data.readline() ### first line ========
        nl = 0
        while len(line)>0:
            ### Starting new bioset                    
            while self.partition in line and len(line)>0:
                line = input_data.readline() ### bioset title/ID
                logging.info(' (GetSnpList main) - title %s' % line )
                while self.partition not in line and len(line)>0 :
                    line = input_data.readline()
                ### starting SNP table
                if self.partition in line:
                    line = input_data.readline() ### snp_table_header ot table
                    self.util.getColumnList(line)
                line = input_data.readline() ### delimiter before SNP data
                if self.snp_part in line:
                    line = input_data.readline() ### 1st SNP info
                    self.recoverSnpData(line, self.util.column_list)
            line = input_data.readline()  
            if self.partition not in line and len(line)>0:
                self.recoverSnpData(line, self.util.column_list)
            nl +=1
            if nl % 100000 ==0:
                print nl
                sys.stdout.flush()
    
        # outputData   
        self.util.writeOutput("snp_List",self.snp_list)
        self.writeListSnps()


    ''' csvDataFiles '''
    def csvDataFiles(self):
        logging.debug(' Function:  csvDataFiles' )
        ### get batch file names
        batches_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if  i.endswith('_data_v7.csv')] 
        
        ### get old list of snp
        db_data = nextbio_db['dbsnp'].find({},{'dbsnp':1,'_id':0}) # find all
        self.old_list =[]
        for d in db_data :
            self.old_list.append(d['dbsnp'])

        
        output = open("%s/snp_list_%s.txt"% (config.OUTPUTPATH,nfile),'w')
        new_list = []
        for batch in batches_list:
            self.util.batchname = batch
            print batch
            csvfile = open('%s/%s' % (config.DATAPATH,batch), 'rb') 
            spamreader = csv.reader(csvfile)
            for row in spamreader:
                if len(row)>2:
                    if row[2] not in self.old_list:
                        new_list.append(row[2])
                        
         output.write("%s\n" % snp)
                output_count.write("%s\t%s\n" % (snp,len(self.snp_list[snp]['id'])))
                nsnp += 1
                if nsnp % 5000 == 0:
                    output.close()
                    nfile += 1
                    output = open("%s/snp_list_%s.txt"% (config.OUTPUTPATH,nfile),'w')
        output.close() 
        output_count.close()





""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'NextBio_GetSnpList'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = GetSnpList()
    clean.csvDataFiles()
#     clean.main4disease()
#     clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 