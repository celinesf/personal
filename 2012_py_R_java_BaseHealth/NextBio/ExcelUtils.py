#!/usr/bin/env python

"""
     Utility functions to obtain read and write on excel files
     06/19/13- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 



import logging, xlrd, copy
from NextBioUtils import NextBioUtils

class ExcelUtils():
    def __init__(self):
        self.util = NextBioUtils() 
        self.partition = "="*30         # this separates meta data from snp information
        self.snp_part  = "-"*30         # this separates column names from data in snp section

############## read excel data ###################
    ''' goto next row while getting key/value on a row
    assumes data from exel file '''
    def getKeyValueRow(self,data, rownum):
        logging.debug(' Function:  getKeyValueRow' )
        rownum+=1
#         print rownum, data[0][rownum]
        key = data[0][rownum].lower()
        if len(key.split(" "))>1:
            key = "%s_%s" % (key.split(" ")[0],key.split(" ")[1])
        value  = data[1][rownum]
        return rownum, key, value
    
    ''' extract data from a row from SNP table '''
    def getSnpRow(self, header, row, data):
        logging.debug(' Function:  getSnpRow' )
        csnp = header.index("SNP" )
        dbsnp = row[csnp]
        data[dbsnp]={}
        self.line = None
        for num in range(len(header)):
            if len(header[num])>0:
                if self.line is None:
                    self.line =row[num].replace('\n','').replace(' ',"")
                else:
#                     print self.line
                    self.line = '%s\t%s' % (self.line,row[num].replace('\n','').replace(' ',""))
                if header[num].lower().replace('\n','').replace(' ',"") == 'comment':
                    data[dbsnp][header[num].lower().replace('\n','').replace(' ',"")] = row[num].replace('\n','').replace(' ',"_")
                elif header[num].lower().replace('\n','').replace(' ',"") == 'snp_population':
                    data[dbsnp][header[num].lower().replace('\n','').replace(' ',"")] = row[num].replace('\n','')
                else:
                    data[dbsnp][header[num].lower().replace('\n','').replace(' ',"")] = row[num].replace('\n','').replace(' ',"")
                if data[dbsnp][header[num].lower().replace('\n','').replace(' ',"")] == "":
                    data[dbsnp][header[num].lower().replace('\n','').replace(' ',"")] = None
        return data    

    ''' getBiosetMetaData'''
    def getBiosetMetaData(self, key,value,  excel_data,row_num):
        logging.debug(' Function:  getBiosetMetaData - key: %s, nrow: %s' %(key,row_num))
        # skip empty line before new bioset
        row_num, key, value = self.getKeyValueRow(excel_data,row_num)
        if len(key) == 0: 
            row_num, key, value = self.getKeyValueRow(excel_data,row_num)
        self.bioset_id= value # bioset #
#         print self.bioset_id, self.bioset_info
        ### only record accepted biosets by science team
        if self.bioset_info[self.bioset_id]["accepted"].lower() == 'yes': 
            self.nextbio_data[self.bioset_id]=copy.deepcopy(self.bioset_info[self.bioset_id] )
        ### bioset meta data
        while self.partition not in key and len(key)>0: #
            if excel_data[2][row_num] == "" :
                if self.bioset_info[self.bioset_id]["accepted"].lower() == 'yes':## only record accepted biosets by science team
                    self.nextbio_data[self.bioset_id][key.lower()]=value
                    if key not in self.bioset_info[self.bioset_id]:
                        self.nextbio_data[self.bioset_id][key.lower()]=value
                    else:
                        self.util.checkSame(self.bioset_info[self.bioset_id][key], value)
            else : print "issue A", excel_data[2][row_num]
            row_num, key, value = self.getKeyValueRow(excel_data,row_num)
        return row_num, key, value 

    ''' getBiosetSnpData'''
    def getBiosetSnpData(self, key,value,  excel_data,row_num, sheet):
        logging.debug(' Function:  getBiosetSnpData - key: %s, nrow: %s' %(key,row_num))
        
        ### startingB SNP table
        if self.partition in key:
            row_num+=1
            self.column_key[self.bioset_id] = sheet.row_values(row_num)
            if self.bioset_info[self.bioset_id]["accepted"].lower() == 'yes':## only record accepted biosets by science team
                self.nextbio_data[self.bioset_id]["snps"] = {}
        ### snp partition
        row_num, key, value = self.getKeyValueRow(excel_data,row_num)
        ### get all snp data
        if self.snp_part in key:
           
            row_num, key, value = self.getKeyValueRow(excel_data,row_num)
            while  key == self.bioset_id:
                if self.bioset_info[self.bioset_id]["accepted"].lower() == 'yes':## only record accepted biosets by science team
                    self.nextbio_data[self.bioset_id]["snps"] = self.getSnpRow(self.column_key[self.bioset_id], sheet.row_values(row_num),self.nextbio_data[self.bioset_id]["snps"])
                row_num, key, value = self.getKeyValueRow(excel_data,row_num)
        return row_num, key, value 
    
    
    ''' get bioset data including SNP table from NextBio xsl format'''
    def readNextBioDataSheet(self, sheet):
        logging.debug(' Function:  readNextBioDataSheet' )
        self.nextbio_data = {}
        excel_data = []
        ''' get all column'''
        for cnum in range(sheet.ncols ):
            excel_data.append(sheet.col_values(cnum))
        ''' get bioset data per data_row'''
        row_num =  -1
        while row_num+1 < len(excel_data[0]):
            row_num, key, value = self.getKeyValueRow(excel_data,row_num) ### first bioset partision line
            ### start new bioset
            while self.partition in key and len(key)>0 and row_num+1 <len(excel_data[0]): 
                ### get bioset meta data
                row_num, key, value = self.getBiosetMetaData(key,value ,excel_data,row_num)
                ### get snp data
                row_num, key, value = self.getBiosetSnpData(key,value ,excel_data,row_num, sheet)
                 
    
    ''' get data from informative sheets in execl file currated by science team'''
    def readExcelData(self,fname):
        logging.debug(' Function:  readExcelData' )
        self.column_key ={}
        infile = xlrd.open_workbook(fname)
        ''' get bioset info '''
        self.bioset_info,  self.header = self.readBiosetInfoSheet(infile.sheet_by_name("info"))

        ''' get bioset+ SNP data '''
        self.readNextBioDataSheet(infile.sheet_by_name("data")) 
        return self.nextbio_data, self.bioset_info,  self.header, self.column_key 

    
    ''' get information as specified by science team'''
    def readBiosetInfoSheet(self, sheet):
        logging.debug(' Function:  readBiosetInfoSheet' )
        header =[]
        data = {}
        data_tmp = {}
        for cnum in range(sheet.ncols ):
            tmp = sheet.col_values(cnum) 
            col_name = tmp[0]
            header.append(col_name)
            tmp.remove(col_name)
            data_tmp[col_name] = tmp
        for bnum in range (len(data_tmp['bioset_id'])):
            bioset = data_tmp['bioset_id'][bnum]
            print bioset
            data[bioset]={}
            for tag in data_tmp:
                data[bioset][tag] = data_tmp[tag][bnum]
        return data, header
    
#################### write to excel ##############   
          
    ### writeExcelRow
    def writeExcelRow(self,sheet,row,key,value):
        sheet.write(row,0,key )
        if value is not None:
            sheet.write(row,1, value )
        row += 1
        return row
    
    ### writeExcelSnpTableData ###
    def writeExcelSnpTableData(self, sheet,row, header, data):
        logging.debug(' Function:  writeExcelSnpTableData')
        for snp in data:
            c = 0
            for h in header:
                if h in data[snp] and data[snp][h] is not None:
                    sheet.write(row,c,data[snp][h])
                else:
                    sheet.write(row,c,'')
                c+=1
            row += 1
        sheet.write(row,0,"=======================================================================================================================================================================================================================================" )
        row += 1
        sheet.write(row,0,"")
        row += 1
        return row, sheet
    
    ### writeExcelSnpTableHeader ###
    def writeExcelSnpTableHeader(self, sheet,row, header):
        logging.debug(' Function:  writeExcelSnpTableHeader')
        c = 0
        new_header = []
        for h in header:
            if h != "":
                h = h.replace('/n','').strip()
                sheet.write(row,c,h)
                new_header.append(h)
                c+=1
        row += 1
        sheet.write(row,0,"------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
        row += 1
        return row, sheet, new_header
    
    ### writeInfoSheet ###
    def writeInfoSheet(self, sheet,row, header, data):
        logging.debug(' Function:  writeInfoSheet')
        if row == 0:
            c=0
            for key in header:
                sheet.write(row,c,key)
                c+=1
            row +=1 
        if row >0:
            c=0
            for key in header:
                sheet.write(row,c,data[key])
                c+=1
            row +=1 
        return row, sheet
    
    ### write in Nextbio bioset key balye information data###
    def writeNextBioBiosetInfo(self,  sheet, row,data):
        logging.debug(' Function:  writeNextBioBiosetInfo ')
        row= self.writeExcelRow(sheet,row,"BIOSET_ID" ,data["BIOSET_ID"])
        row= self.writeExcelRow(sheet,row,"PMID" ,data["PMID"])
        if "BIOSET TITLE" in data:
            row= self.writeExcelRow(sheet,row,"BIOSET TITLE" ,data["BIOSET TITLE"])
        if "SAMPLE NUMBER" in data:
            row= self.writeExcelRow(sheet,row,"SAMPLE NUMBER" ,data["SAMPLE NUMBER"])
        if "COMPARISON" in data:
            row= self.writeExcelRow(sheet,row,"COMPARISON" ,data["COMPARISON"])
        if "TAG" in data:
            row= self.writeExcelRow(sheet,row,"TAG" ,data["TAG"])
        if "BIOSET SUMMARY" in data:
            row= self.writeExcelRow(sheet,row,"BIOSET SUMMARY" ,data["BIOSET SUMMARY"])
        if "ANALYSIS SUMMARY" in data:
            row= self.writeExcelRow(sheet,row,"ANALYSIS SUMMARY" ,data["ANALYSIS SUMMARY"])
        
        if "PLATFORM" in data:
            row= self.writeExcelRow(sheet,row,"PLATFORM" ,data["PLATFORM"])
        if "META-ANALYSIS" in data:
            row= self.writeExcelRow(sheet,row,"META-ANALYSIS" ,data["META-ANALYSIS"])
        if "REFERENCE POPULATION" in data:
            row= self.writeExcelRow(sheet,row,"REFERENCE POPULATION" ,data["REFERENCE POPULATION"])
        row= self.writeExcelRow(sheet,row,"BATCH" ,data["BATCH"])
        row= self.writeExcelRow(sheet,row,"ACCEPTED" ,None)
        row= self.writeExcelRow(sheet,row,"COMMENT" ,None)
        sheet.write(row,0,"=======================================================================================================================================================================================================================================" )
        row+= 1
        return row, sheet


        
#         self.batchname = ""             # batch file provided by nextbio
#         self.previous_batch ={}
#         self.issue = False
#         self.current_batch ={}
#         self.key_list = {}
#         self.key_combination = []
#         self.column_list = []
#         self.key_comb_list = None
#         self.column_name_dict = {}
#         self.pmid_list = {}
#         self.meta_sp   = [":", "="]     # this separates fields whithin meta data
#         self.vs        = ["_vs_", "vs"] # this separates comparison and tag fields

#     ''' checkSame
#     compare two value'''
#     def checkSame(self, data1, data2):
#         logging.debug(' Function:  checkSame - data1: %s, data2: %s' %(data1, data2))
#         if data1 != data2 :
#             if data1 != '' and  data2 != '' and  data1 is not None and  data2 is not None  :
#                 logging.warning(' DATA DONT FIT (checkSame) - data1: %s, data2: %s' %(data1, data2))
#                 return False
#         return True
# 
#     ### writeExcelRow
#     def writeExcelRow(self,sheet,row,key,value):
#         sheet.write(row,0,key )
#         if value is not None:
#             sheet.write(row,1, value )
#         row += 1
#         return row
# 
#     ''' keyToUpperCase '''
#     def keyToUpperCase(self, data):
#         logging.debug(' Function:  keyToUpperCase ')
#         new_data = {}
#         for key in data:
#             new_key = key.upper()
#             new_data[new_key] = data[key]
#             if type(new_data[new_key]) == dict:
#                 new_data[new_key] = self.keyToUpperCase(new_data[new_key] )
#         return new_data
# 
# #     ''' isSameAllele '''
# #     def isSameAllele(self, allele1, allele2):
# #         logging.debug(' Function:  isSameAllele - allele1:%s, allele2:%s' % (allele1, allele2))
# #         if allele1 == allele2:
# #             return 1
# #         else:
# #             return 0
# 
# #     ''' reverse
# #     '''
# #     def reverseAlleles(self, alleles):
# #         logging.debug(' Function: reverseAlleles - %s' % alleles)  
# #         reversed_allele = []
# #         for index in range(0,len(alleles)):
# #             allele= alleles[index ]
# #             if allele == 'G':
# #                 reversed_allele.append('C')
# #             elif allele == 'C':
# #                 reversed_allele.append('G')
# #             elif allele == 'A':
# #                 reversed_allele.append('T')
# #             elif allele == 'T':
# #                 reversed_allele.append('A')
# #             elif allele == '-':
# #                 reversed_allele.append('-')
# #             else :
# #                 self.warnMe('warning', ' CANT REVERSE ALLELE %s (reverseAlleles) - ' % (allele,alleles))  
# #         if len(reversed_allele) == 1:
# #             reversed_allele = reversed_allele[0]
# #         return reversed_allele

   
   
    
#     ''' warnMe '''   
#     def warnMe(self, flag, comment):
#         if flag == 'info':
#             logging.info(comment)
#             print comment
#         elif flag== 'warning':
#             logging.warning(comment)
#             print "%s -%s" % (flag.upper(),comment)
#         elif flag== 'critical':
#             logging.critical(comment)
#             print "%s -%s" % (flag.upper(),comment)
#         elif flag== 'error':
#             logging.error(comment)
#             print "%s -%s" % (flag.upper(),comment)
#         else:
#             logging.info(comment)
#             print "%s -%s" % (flag.upper(),comment)
#       
#     ''' get map for typo corrections '''
#     def getTermMap(self,fname):
#         logging.debug(' Function:  getTermMap' )
#         f = open(fname,'r')
#         typos = f.read()  
#         f.close()
#         return json.loads(typos)
#         
#       
#     ''' remove quote from key or value 
#         Issue from batch 7 
#     '''
#     def checkQuote(self,line):
#         logging.debug(' Function:  checkQuote' )
#         s=line.encode("utf-8")
#         if "\"" in s: 
#             s = line.split("\"")
#             if len(s) == 2:
#                 if len(s[1]) == 0 :
#                     logging.warning(' REMOVE QUOTE s[1]=0 (checkQuote) - %s - %s\n%s' % (self.batchname,line, s))
#                     s =   s[0]
#                 elif len(s[0])== 0 :
#                     logging.warning(' REMOVE QUOTE s[0]=0 (checkQuote) - %s - %s\n%s' % (self.batchname,line, s))
#                     s =   s[1]  
#                 else:
#                     self.warnMe('error', ' ERROR QUOTE (checkQuote)  - %s - %s\n%s' % (self.batchname, line,s))
#                     s= line
#             else:
#                 s1 = None
#                 for w in s:
#                     if s1 is None:
#                         s1=w
#                     else:
#                         s1 = "%s'%s" % (s1,w)
#                 s = s1
#                 logging.warning(' CHANGED QUOTE (checkQuote) - %s - %s\n%s' % (self.batchname, line,s))
#         return s
#     
#     '''find key before : or ='''
#     def findKey(self,line):
#         logging.debug(' Function:  findKey' )
#         key = None
#         value = None
#         for sep in self.meta_sp: 
#             if sep in line:
#                 try:
#                     key, value = line.split(sep,1)
#                 except Exception, e:
#                     self.warnMe('critical',' NO KEY VALUE? (findKey) %s' % line )
#         return key, value
    
#     ''' add new item to list '''
#     def addToList(self,item_list, item, new_value, old_value):
#         logging.debug(' Function:  addToList %s' % item )
#         if item not in item_list:
#             item_list[item] = new_value
#         if old_value is not None:
#             item_list[item][old_value] += 1   
#         else:
#             item_list[item] += 1   
#         return item_list


#     ''' get key and value from a line in meta data and count # occurence in batch''' 
#     def getKeyValueCount(self,key_map,  line):
#         logging.debug(' Function:  getKeyValue' )
#         line = line.strip()
#         key, value = self.findKey(line)
# 
#         ''' record list of keys / value'''
#         self.current_batch["BATCH"]=self.batchname
#         if key :                   
#             key = key.strip().upper()
#             key = self.checkQuote(key)
#             
#             if key in key_map:
#                 key = key_map[key]
#             else:
#                 self.warnMe('critical', "I COULD NOT FIND MAPPING FOR KEY %s (getKeyValueCount)- %s" % (key, self.current_batch))
#             
#             ''' add new key '''
#             self.key_list = self.addToList(self.key_list, key, {"count" : 0, "value" : {}, "batch":{}}, "count")
# 
#             ''' record value'''
#             value = value.strip()
#             value = self.checkQuote(value)
#             self.key_list[key]["value"] = self.addToList(self.key_list[key]["value"], value, {"count":0,"batch":{}}, "count")
#             self.key_list[key]["value"][value]["batch"] = self.addToList(self.key_list[key]["value"][value]["batch"], self.batchname, 0, None)  
#         
#             self.current_batch[key]["value"] = self.addToList(self.current_batch[key]["value"], value, 0, None)    
#             self.current_batch[key]=value
#             
#             ''' record batchname '''
#             self.key_list[key]["batch"] = self.addToList(self.key_list[key]["batch"], self.batchname, 0, None)
#                         
#             self.checkKeyValue(key, value)
#         return key
               
#     ''' get key and value from a line in meta data ''' 
#     def getKeyValue(self, key_map, line):
#         logging.debug(' Function:  getKeyValue' )
#         line = line.strip()
#         key, value = self.findKey(line)
# 
#         ''' record list of keys / value'''
#         self.current_batch["BATCH"]=self.batchname
#         if key :                   
#             key = key.strip().upper()
#             key = self.checkQuote(key)
#             
#             if key in key_map:
#                 key = key_map[key]
#             else:
#                 self.warnMe('critical', "I COULD NOT FIND MAPPING FOR KEY %s (getKeyValue)- %s" % (key, self.current_batch))
#             
#             ''' add new key '''
#             self.key_list = self.addToList(self.key_list, key, {"count" : 0, "value" : {}, "batch":{}}, "count")
# 
#             ''' record value'''
#             value = value.strip()
#             value = self.checkQuote(value)
#             self.key_list[key]["value"] = self.addToList(self.key_list[key]["value"], value, 0, None)
#         
#             self.current_batch[key]=value
#             
#             ''' record batchname '''
#             self.key_list[key]["batch"] = self.addToList(self.key_list[key]["batch"], self.batchname, 0, None)
#                         
#             self.checkKeyValue(key, value)
#         return key
            
            
#     ''' check key and value if duplicated and came from previous bioset '''
#     def  checkKeyValue(self,key,value):
#         logging.debug(' Function:  checkKeyValue, %s: %s' % (key,value) )
#         if key not in self.key_combination:
#                 self.key_combination.append(key)
#         else:
#             self.issue = True
#             logging.warning("REPEATED KEY IN ONE BIOSET (checkKeyValue) %s\t%s\t%s:%s" % (key, self.batchname, key,value))
#             if self.current_batch[key]["count"] >1 and self.current_batch[key]["value"][value] == self.current_batch[key]["count"]:
#                 logging.warning("VALUE FOUND TWICE IN ONE BIOSET (checkKeyValue): %s" %  value)
#                 if value in self.previous_batch[key]["value"]:
#                     logging.warning("VALUE FOUND IN PREVIOUS BIOSET (checkKeyValue): %s" %  value)
#                 else:
#                     logging.error("VALUE FOUND NOT IN PREVIOUS BIOSET FOR SAME KEY (checkKeyValue): %s:%s\t%s" %  (key, value, self.previous_batch[key]["value"]))
#             else:
#                 self.warnMe("error", "DIFFERENT VALUE FOUND FOR SAME KEY (checkKeyValue): %s: %s\n%s" %  (key,value,self.current_batch[key]["value"]))
#                  
#     ''' get column names for SNP data
#     ''' 
#     def getColumnNames(self, line):         
#         logging.debug(' Function:  getColumnNames' )
#         line = line.strip()
#         for word in line.split('\t'):
#             if word not in self.column_name_dict:
#                 self.column_name_dict[word] = {"count" : 0, "batch":{}}
#             self.column_name_dict[word]["count"] += 1
#             if self.batchname not in self.column_name_dict[word]["batch"]:
#                 self.column_name_dict[word]["batch"][self.batchname] = 0
#             self.column_name_dict[word]["batch"][self.batchname] += 1
#             
#     ''' get column names for SNP data
#     ''' 
#     def getColumnList(self, line):         
#         logging.debug(' Function:  getColumnNames' )
#         line = line.strip()
#         self.column_list = []
#         for word in line.split('\t'):
#             if word not in self.column_list:
#                 self.column_list.append(word)
#             else:
#                 word = '%s2' % word
#                 self.column_list.append(word)
#                 self.warnMe('warning', ' WORD REPEATED %s (getColumnList) %s' % (word, line))
            
#     ''' get unique PMID and they bioset IDs
#     ''' 
#     def getPMID(self, line):         
#         logging.debug(' Function:  getPMID' )
#         line = line.strip()
#         data = line.split('\t')
#         pmid = data[1]
#         self.current_batch["PMID"] = pmid
#         self.current_batch["BIOSET_ID"] = data[0]
#         bioset_id = data[0]
#         if pmid not in self.pmid_list:
#             self.pmid_list[pmid] = {"count":0,"id":{}}
#         self.pmid_list[pmid]["count"] += 1
#         if bioset_id not in self.pmid_list[pmid]["id"]:
#             self.pmid_list[pmid]["id"][bioset_id] = 0
#         self.pmid_list[pmid]["id"][bioset_id] +=1    
#             
#     def writeOutput(self, filename,data):
#         logging.debug(' Function:  writeOutput %s:' % filename)
#         f = open("%s/%s.json" % (config.OUTPUTPATH,filename), "w")
#         f.write(json.dumps(data, indent=4))
#         f.close()
        
#     ''' write unique key of list '''
#     def writeList(self, filename,data):
#         logging.debug(' Function:  writeList %s:' % filename)
#         total = 0
#         f = open("%s%s.json" % (config.OUTPUTPATH,filename), "w")
#         for d in data:
#             f.write("%s\t%s\n" % (d, data[d]["count"])) 
#             total +=  data[d]["count"]
#         f.write("TOTAL\t%s\n" % (total))
#     
#     ''' add new key combination in list ''' 
#     def addNewKeyCombination(self):
#         logging.debug(' Function:  addNewKeyCombination %s' %(self.key_combination) )
#         if len(self.key_combination)>0:
#             ''' init combination list '''            
#             if self.key_comb_list is None:
#                 self.key_comb_list = {}
#                 self.key_comb_list = self.addToList(self.key_comb_list, "0", {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
#                 self.key_comb_list["0"]["batch"] = self.addToList(self.key_comb_list["0"]["batch"], self.batchname, 0, None)
# #                 self.key_comb_list["0"] = {"keys":self.key_combination,"count":0,"batch":{}}
#             else:
#                 found = 0
#                 comb_num =""
#                 for comb in self.key_comb_list:
#                     comb_num = comb
#                     it = 0
#                     for key in self.key_comb_list[comb]["keys"]:
#                         if key in self.key_combination:
#                             it += 1
#                         else:
#                             break
#                     if it == len(self.key_combination) and it == len(self.key_comb_list[comb]['keys']):
#                         found = 1
#                         self.key_comb_list = self.addToList(self.key_comb_list, comb, {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
#                         self.key_comb_list[comb]["batch"] = self.addToList(self.key_comb_list[comb]["batch"], self.batchname, 0, None)
#                         break
#                 if found == 0:
#                     comb_num = str(len(self.key_comb_list))
#                     self.key_comb_list = self.addToList(self.key_comb_list, comb_num, {"keys":self.key_combination,"count" : 0, "batch":{}},"count")
#                     self.key_comb_list[comb_num]["batch"] = self.addToList(self.key_comb_list[comb_num]["batch"], self.batchname, 0, None)
#         self.key_combination = []
     
    
######################### READ EXCEL DATA FILE #####################

#  