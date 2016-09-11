#!/usr/bin/env python

"""
     Simplify dbsnp json structure
     06/15/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) get dbsp in json format
    2_ simplify json structure
"""

import logging, json, copy,re
import DbSnpConfig as config
from DbSnpUtils import DbSnpUtils
from GatherSnpInfo import GatherSnpInfo


class SimpleDbsnpInfo():
    def __init__(self):
        self.util = DbSnpUtils()
        self.dbsnp_map =self.util.getTermMap(config.XMLMAP )
        self.gatherInfo = GatherSnpInfo()
        
    ''' getDataDict 
        loop over key/objects in  a dictionary
    '''
    def getDataDict(self,data,snp_info):
        logging.debug(' Function: getDataDict ')
        for key in data:
            self.chooseDataType(key, data, snp_info)

    ''' getDataList 
        loop over key/objects in  a list
        each new object is added to the information list
    '''
    def getDataList(self,data,template_doc):
        logging.debug(' Function: getDataList ')
        list_info = []
        for index in range(0,len(data)):
            new_doc=copy.deepcopy(template_doc)
            self.chooseDataType(index, data, new_doc)
            if new_doc not in list_info:
                list_info.append(new_doc)
        return list_info


    ''' chooseDataType '''
    def chooseDataType(self,key, data,snp_info):
        logging.debug(' Function: chooseDataType ') 
        if key not in self.key_list:
            self.key_list.append(key)
        else:
            if key in snp_info:
                if snp_info[key] is not None and type(snp_info[key]) != dict and type(snp_info[key]) != list :
                    self.util.warnMe('warning', ' REPATED KEY (chooseDataType) - key:%s, value:%s, data:%s' %( key ,snp_info[key], data[key]))
                    return None
        if type(data[key]) == dict :
            if key in snp_info:
                logging.info("DICT key FOUND (chooseDataType) - %s" % key)
                self.getDataDict(data[key],snp_info[key])
            else:
                logging.info("DICT key NOT found (chooseDataType) - %s" % key)
                self.getDataDict(data[key],snp_info)
        elif type(data[key]) == list:
            if key in snp_info:
                logging.info("LIST key FOUND (chooseDataType) - %s" % key)
                if type(snp_info[key]) == list:
                    snp_info[key] = self.getDataList(data[key],copy.deepcopy(snp_info[key][0]) )
                else:
                    snp_info[key] = self.getDataList(data[key],copy.deepcopy(snp_info[key]) )
            else :
                logging.info("LIST key NOT found (chooseDataType) - %s" % key)
        else:
            if key in snp_info :
                logging.info("UNICODE key FOUND (chooseDataType) - %s" % key)
                snp_info[key] = data[key]
            else:
                logging.warning("UNICODE key NOT found (chooseDataType) - %s" % key)

    ''' fillKeyValue '''
    def fillKeyValue(self,data,snp_info):
        logging.debug(' Function: fillKeyValue ') 
        for key in data:
            self.chooseDataType(key, data, snp_info)


    ''' simplifyStructure '''
    def simplifyStructure(self, data):
        logging.debug(' Function: simplifyStructure') 
        snp_info = copy.deepcopy(self.dbsnp_map)
        self.key_list =[]
        self.fillKeyValue(data,snp_info)
        dbsnp_info = self.gatherInfo.gatherSnpData(copy.deepcopy(snp_info))
        return dbsnp_info


    ''' getDdSnpList '''
    def getDdSnpList(self):
        logging.debug(' Function: getDdSnpList') 
        self.dbsnp_info = {}
        dbsnp_dic = self.util.getTermMap(config.OUTPUTJSON )
        ns = 0
        for key in dbsnp_dic:
            snp_info = copy.deepcopy(self.dbsnp_map)
            self.fillKeyValue(dbsnp_dic[key],snp_info)
            self.dbsnp_info[key] = copy.deepcopy(snp_info)
            ns+=1
            if ns %100 == 0:
                print ns

                       
    ''' find list of keys '''
    ''' find list of columns '''
    def main(self):
        logging.debug(' Function:  SimpleDbsnpInfo main' )

        self.getDdSnpList()
        self.util.writeOutput('%s/dbsnp_simple.json' % config.OUTPUTPATH, self.json_data)



""" MAIN function """    
if __name__ == "__main__":
    
    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'DBSNP_SimpleDbsnpInfo'), filemode='w',
                        level=logging.INFO,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = SimpleDbsnpInfo()
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 