#!/usr/bin/python
'''
Created on October 9, 2012
@author: celine

update 2.1 genomics documents (in json format)
- add real and relative position
- add gene name
- update/clean effect
'''
import json
import pymongo
import copy



# tiger_dev
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"

#tiger_qa
#mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.229.194/admin"

#mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.140.237/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genophen30"]
GENETICS= CW["genetics"]
DIR = '/Users/celine/Box Documents/Celine/ScriptInOut/mongodb-update'


def getHapListFromGeneticsRep():
    data = GENOPHEN30['genetics.rep'].find()
    hap_list = {}
    for disease in data:
        for doc in disease['snps']:
            if 'rs' not in doc['snp']:
                if doc['snp'] not in hap_list:
                    hap_list[doc['snp']] = {'disease':disease, 'data':[]}
                hap_list[doc['snp']]['data'].append(doc)
    return hap_list
                    
def getHapSnpFromIhap(hap_list ):
    snp_list = {}
    for hap in hap_list:
        data = GENETICS['idhap'].find({'dbhap':hap})
        for doc in data:
            if doc['dbsnp'] not in snp_list:
                snp_list[doc['dbsnp']]= { 'dbhap':[],"minor_allele":doc["minor_allele"], "major_allele":doc["major_allele"]}
            snp_list[doc['dbsnp']]['dbhap'].append(hap)
            if doc["minor_allele"] != snp_list[doc['dbsnp']]['minor_allele']:
                print 'issue A'
            if doc["major_allele"] != snp_list[doc['dbsnp']]['major_allele']:
                print 'issue B'
    return snp_list

def checkAllelesDbsnp(snp_list):
    for snp in snp_list:
        print snp
        data = GENOPHEN30['dbsnp'].find({'dbsnp':snp})
        for doc in data:
            if snp_list[snp]["minor_allele"] not in doc['alleles'] or  snp_list[snp]["major_allele"] not in doc['alleles']:
                print 'issue wront alleles' , snp ,doc['alleles'] , snp_list[snp]

             
# hap_list = getHapListFromGeneticsRep()
# # print hap_list
# # print json.dumps(hap_list,indent=4)
# snp_list = getHapSnpFromIhap(hap_list)
# # print 'B', json.dumps(snp_list,indent=4)
# checkAllelesDbsnp(snp_list)

rs706793

