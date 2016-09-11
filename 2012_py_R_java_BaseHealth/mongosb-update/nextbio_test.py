#!/usr/bin/python
'''
Created on April 23, 2013
@author: celine

'''
import pymongo, logging
import json


#ares_qa
#mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.229.194/admin"

#### tiger_dev
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genetics"]



def findUnicKey():
    search = GENOPHEN30.nextbio_raw.find({});
    keys = {}
    for doc in search:
        for key in doc:
            if key not in keys:
                print 'new key', key
                keys[key] = 0
            keys[key] += 1   
             
    print json.dumps(keys, indent=4)
            

        
def findUnicComparison():
    print 'findUnicTag'
    search = GENOPHEN30.nextbio_raw.find({},{ "comparisons":1,"comparison":1,"tag":1,"tags":1,"_id":0,"bioset title":1,"bioset summary":1,"analysis summary":1,"batch_name":1});

    unic_comp = {}
    for doc in  search:
        comp_flag  = doc['comparison']
        if "comparisons" in doc:
            print json.dumps(doc)
        tag_flag = None
        if "tag" in doc:
            tag_flag = doc['tag']
        elif "tags" in doc:
            tag_flag = doc['tags']
        else:
            print  'not tags ',json.dumps(doc)
        if doc['comparison'] not in unic_comp:
            unic_comp[comp_flag] = {"total_count":0,"tags":{}}
        unic_comp[comp_flag]["total_count"]+=1
        if tag_flag  is not None:
            if tag_flag not in unic_comp[comp_flag]['tags']:
                unic_comp[comp_flag]['tags'][tag_flag] = {"tag_count":0,"summary":{}}
            unic_comp[comp_flag]['tags'][tag_flag]["tag_count"] += 1    
        title = None
        if 'bioset title' in doc:
            title =  doc['bioset title']
        else:
            print 'not title ',json.dumps(doc)
        summary = None
        if "bioset summary" in doc:
            summary = doc["bioset summary"]
            
        elif "analysis summary" in doc:
            print 'ANALYSIS ',json.dumps(doc)
            summary = doc["analysis summary"]
        else:
            print 'not summary ',json.dumps(doc)

        if tag_flag is not None:
            if summary not in unic_comp[comp_flag]['tags'][tag_flag]['summary']:
                unic_comp[comp_flag]['tags'][tag_flag]['summary'][summary]={"count":0,"title":{}}
            unic_comp[comp_flag]['tags'][tag_flag]['summary'][summary]['count'] += 1
            unic_comp[comp_flag]['tags'][tag_flag]['summary'][summary]['title'][title] = '' 
        else :
            unic_comp[comp_flag]['title_RECOVER'] = title
        
    print json.dumps(unic_comp, indent=4)

        
################ MAIN ###############
if __name__ == "__main__":



    findUnicKey()
    print 'DONE'

