#!/usr/bin/python
'''
Created on January 28, 2012
@author: celine

Check the gid match for dev and QA in 
Sasha.authentication
profile
riskoutput

'''
import pymongo, logging
import json
import copy


class GidChecks():
    """ Init FUNCTION """ 
    def __init__(self,conn ):
        logging.debug(' Function: __init__' )
        self.CW = pymongo.Connection(conn, 27017);
        self.GENOPHEN30 = self.CW["genophen30"]
#        mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@107.21.17.73/crm"
#        self.CW_CRM = pymongo.Connection(mongo_cnxn_string, 47954); 
#        self.CRM = CW_CRM["crm"]
    
    def update_profile(self,sasha):
        logging.debug(' Function: update_profile' )
        gid = sasha["static"]["gid"]
        find_profile= self.GENOPHEN30.profile.find({"gid":gid});    
        size = 0
        for p in find_profile:
            size += 1
        if size == 1:
            logging.info(' update_profile -- I found gid %s' % gid)
        elif size == 0 :
            logging.warning(' update_profile -- I could not find gid %s (username: %s)' % (gid,sasha["static"]["username"]) )
            print "size=", size, gid,sasha["static"]["username"]
            find_template= self.GENOPHEN30.profile.find({"gid":"0"},{"_id":0}); 
            for t in find_template:
               t['gid'] = gid
#               self.GENOPHEN30.profile.insert(t)
               
    def update_riskoutput(self,sasha):
        logging.debug(' Function: update_riskoutput' )
        gid = sasha["static"]["gid"]
        find_riskoutput= self.GENOPHEN30.riskoutput.find({"gid":gid});    
        size = 0
        for r in find_riskoutput:
            size += 1
        if size == 1:
            logging.info(' update_riskoutput -- I found gid %s' % gid)
        elif size == 0 :
            logging.warning(' update_riskoutput -- I could not find gid %s (username: %s)' % (gid,sasha["static"]["username"]) )
            print "size=", size, gid,sasha["static"]["username"]
            find_template= self.GENOPHEN30.riskoutput.find({"gid":"0"},{"_id":0}); 
            for t in find_template:
               t['gid'] = gid
               print gid, sasha["static"]["username"]
#               print t
#               self.GENOPHEN30.riskoutput.insert(t)
           
        
    def check_gid(self):
        logging.debug(' Function: check_gid' )
        find_sasha= self.GENOPHEN30.sasha.authentication.find({},{"static.username":1,"static.gid":1,"_id":0});        
#        for member in find_sasha:
#            self.update_profile(member)
#            self.update_riskoutput(member)


        
################ MAIN ###############
if __name__ == "__main__":

    logging.basicConfig(filename='log' , filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')

    ######## dev
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.142.17/admin"
    ######## qa
#    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.140.237/admin"

    gid_checks = GidChecks(mongo_cnxn_string)
    gid_checks.check_gid()

    print 'DONE'

