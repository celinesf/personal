#!/usr/bin/env python

######################## STOP - BEFORE STARTING ######################
# a.    CalculateModNmod.R
#     i.    in getContinuouslRiskFactors
#         1.    add any weird risk factors (i.e, BMI values are 10*BMI)
# b.    ConfigModNmod.py
#     i.    Configuration file for GetModNonModRiskFactors
#         1.    LOG_FILENAME 
#         2.    OUTPUTPATH 
#         3.    ROUTPUTPATH = to OUTPUTPATH but white space as '\ '
#         4.    DATAPATH : path to .bib files
#         5.    MONGO : connection string
# c.    main:
#     i.    disease_name : key specifying disease in database
#     ii.    version 
#     iii.     update_mongo = True/False -> do would want to pdate DB


"""
    Prepare algorithm documents for Mod Non-Mod risk factors for a disease
    06/09/13: version 1.0
    06/11/13: version 1.1: functionalized
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 
"""
    1) Input is <disease_name>_<risk_factor>.bib - DONE
    2)  Data recorded in Mongodb 
    3) Genophen yes extracted / run R on them
    4) output data into Algorithm, riskouput/ general/attributes, reference   
    5) clear up divide with different classes and functions
"""

import logging, os, pymongo
import copy
import subprocess
import datetime
import ConfigModNmod as config
from MleUtil import MleUtil
from BibUtil import BibUtil
from InputOutputUtil import InputOutputUtil

########### connect to mongo  ###############
# mongo_cnxn = "mongodb://adm1n:adm1npw0rd@54.245.239.40:27017/admin" ### tiger
# # mongo_cnxn = "mongodb://adm1n:adm1npw0rd@54.244.24.253:27017/admin" ### mars
conn = pymongo.Connection(config.MONGO, 27017)
comprehensive_db = conn["comprehensive"]
genophen_db = conn["genophen30"]

DATE = datetime.datetime.now().strftime("%m-%d-%Y")

class GetModNonModRiskFactors():
    def __init__(self,disease_name, version, update_mongo, paf):
        self.disease_name = disease_name
        self.version = version
        self.update_mongo = update_mongo
        self.paf = paf
        self.risk_factor_data = [] # record data for all papers/ all risk factors
        self.util = InputOutputUtil(config.OUTPUTPATH,self.disease_name )

##################### function to run R ###################
    ''' run R script to get ModNmod data '''
    def calculateModNMod(self):
        logging.debug(' Function:  calculateModNMod')
        r_output = '%s/%s_Rconsole.txt' % (config.ROUTPUTPATH ,self.disease_name)
        script_filename = 'CalculateModNmod.R'
        param_filename = self.disease_name
        process = subprocess.Popen(
                ["R --vanilla --args %s %s %s/ < %s > %s" % (param_filename, config.MONGOHOST,config.ROUTPUTPATH, script_filename,r_output) ], shell=True) 
        print "R --vanilla --args %s %s %s/ < %s > %s" % (param_filename, config.MONGOHOST,config.ROUTPUTPATH, script_filename,r_output) 
        process.wait()
       
        
    ''' getRiskFactorAlgorithm '''
    def getRiskFactorAlgorithm(self):
        logging.debug(' Function:  getRiskFactorAlgorithm')
        r_data = open('%s/%s_data_modNmod.xls' % (config.OUTPUTPATH ,self.disease_name),'r')
        line = r_data.readline().strip()
        header = line.split('\t')
        ### indexes of valuable information
        risk_id = header.index('risk_id')
        value_as_int = header.index('value_as_int')
        log_or = header.index('log_or')
        normal = header.index('normal')

        risk_factors=[]
        ### loop on risk factor values
        while len(line)>0:
            line = r_data.readline().strip()
            if len(line)>0:
                words = line.replace(' ','').split('\t')# remove extra white space from R data.frame
                if words[value_as_int] == 'NA': # deal with missing data Na to -999
                    words[value_as_int] = int(-999)
                ### add new risk document
                risk_doc = {"factor_value": float(words[log_or]), "attr_name": words[risk_id],"attr_value": int(words[value_as_int])}
                risk_factors.append(risk_doc)
        r_data.close()
        return risk_factors
#####################


#####################  mongo util functions #####################  
    ''' getRecordExistingMongoDocument '''
    def getRecordExistingMongoDocument(self, db, collection, query, fields, name):
        logging.debug(' Function:  getRecordExistingMongoDocument- name %s' % (name)) 
        data = db[collection].find(query,fields) 
        for doc in data:
            ## get document date
            if 'date_update' in doc:
                date = doc['date_update']
            else :
                date = doc['date']
            ### output json
            self.util.writeJsonFile("%s_%s_%s" %(name, date,doc['version']), doc)
            return doc

    ''' recordNewMongoDocument in json and mongo'''
    def recordNewMongoDocument(self, filename,new_doc,flag, db, collection,query ):
        logging.debug(' Function:  recordNewMongoDocument- filename %s, flag: %s' % (filename, flag))  
        if flag :
            ### output json
            self.util.writeJsonFile(filename, new_doc)
            ### update algorithms
            if self.update_mongo:
                self.util.warnMe('info', ' mongo %s/%s was updated' %(collection, query))
                db[collection].update(query,new_doc) 
            else:
                self.util.warnMe('warning', ' DID NOT UPDATE mongo %s/%s' %(collection, query))
        else :
            self.util.warnMe('info', ' No changes in %s/%s' %(collection, query))   
            ### output json
            self.util.writeJsonFile(filename, new_doc)
    
    ''' initNewMongoDocument'''
    def initNewMongoDocument(self, document, version, doc_name ):
        logging.debug(' Function:  initNewMongoDocument %s' % doc_name)
        if "name" in document:
            document['name']=doc_name
        if "date_update" in document:
            document['date_update']= DATE
        else:
            document['date']=DATE
        document['version'] = version 
        filename = "%s_%s_%s" %(doc_name, DATE,document['version'])  
        return document, filename, False
#####################


#####################  functions related to updateMongoGeneralAttributes #####################     
    ''' getTemplate '''
    def getTemplate(self, template, doc):
        logging.debug(' Function:  getTemplate')
        if 'whatif_map' in doc :
            if doc['whatif_map'][0]['continuous'] ==  False and template['modifiable_categories'] is None :
                template['modifiable_categories'] = doc
            elif doc['whatif_map'][0]['continuous'] == True and template['modifiable_continuous'] is None :
                template['modifiable_continuous'] = doc
        elif 'whatif_map'  not in doc and template['non_modifiable'] is None :
            template['non_modifiable'] = doc
        return template

    ''' updateWhatIf '''
    def updateWhatIf(self, new_doc,risk_info ):
        logging.debug(' Function:  updateWhatIf %s' % risk_info['risk_id'])
        ### change basic information
        what_if = new_doc['whatif_map'][0]
        what_if['d'] = risk_info['risk_display']
        what_if['group_name'] = risk_info['risk_display']
        what_if['group_name_key'] = risk_info['risk_id']
        what_if['v'] = risk_info['risk_id']
        what_if['unit'] = None
        ### categorical risk factor
        if what_if['continuous'] == False: 
            what_if['categories'] = []
            for cat in risk_info['cat'].split(','):
                what_if['categories'].append(cat)
        ### continuous risk factpr
        else :                             
             what_if['categories'] =[]
             what_if['min'] = risk_info['min_val']
             what_if['max'] = risk_info['max_val']
        self.util.warnMe('warning'," CHECK WHAT IF OF NEW ATTRIBUTE %s" % risk_info['risk_id'])   
        new_doc['whatif_map']=[] 
        new_doc['whatif_map'].append(what_if)
        return new_doc

    ''' getExistingAttributes '''
    def getExistingAttributes(self, attributes_doc):
        logging.debug(' Function:  getExistingAttributes')
        list = []
        temp = {'non_modifiable' : None, 'modifiable_continuous' : None , 'modifiable_categories' : None }
        new_rf = copy.deepcopy(self.risk_genophen)
        ### find which risk factor is new
        for doc in attributes_doc['list']:
            ### get templates
            temp = self.getTemplate(temp, doc)
            
            ### Add existing attributes
            attribute = doc
            for risk in self.risk_genophen:
                if risk == doc['v'] :
                    ### add disease normal value to existing protective factor
                    if new_rf[risk]['protective'] == 'yes': 
                        self.util.warnMe('warning', ' PROTECTIVE EXISTING RISK FACTOR - %s: %s' % ( risk,new_rf[risk]['normal_val']) )
                        attribute['norm'][self.disease_name] = float(new_rf[risk]['normal_val'])
                        flag_changes = True
                    new_rf.pop(risk) ## remove this existing risk factor -> no need to add it
            list.append(attribute)
        return list, temp, new_rf

    ''' getNewAttributes '''
    def getNewAttributes(self, attributes_list,rf_list, template, flag) :
        logging.debug(' Function:  getNewAttributes')
        for new_risk in rf_list:
            if new_risk  != 'lifetime_prevalence' and new_risk  != 'incidence' :
                self.util.warnMe('warning'," NEW ATTRIBUTE (updateMongoGeneralAttributes) - GO CHECK %s in %s " % (new_risk, self.disease_name))   
                if rf_list[new_risk]['risk_type'] == 'modifiable' and rf_list[new_risk]['is_cat'] == 'yes':
                    attribute = copy.deepcopy(template['modifiable_categories'])
                    attribute = self.updateWhatIf(attribute,rf_list[new_risk] )
                    ### TODO udate whatif
                elif rf_list[new_risk]['risk_type'] == 'modifiable' and rf_list[new_risk]['is_cat'] == 'no':
                    attribute = copy.deepcopy(template['modifiable_continuous'])
                    attribute = self.updateWhatIf(attribute,rf_list[new_risk] )
                elif rf_list[new_risk]['risk_type'] == 'non_modifiable':
                    attribute = copy.deepcopy(template['non_modifiable'])
                else:
                    self.util.warnMe('warning', ' Cannot figure out my template (updateMongoGeneralAttributes) - %s' % rf_list[new_risk])
                attribute['v'] = new_risk
                attribute['d'] = rf_list[new_risk]['risk_display']
                if rf_list[new_risk]['protective'] == 'yes': # if protective 
                    self.util.warnMe('warning',' PROTECTIVE NEW RISK FACTOR (updateMongoGeneralAttributes) - %s: %s' % ( new_risk,rf_list[new_risk]['normal_val']) )
                    attribute['norm'][self.disease_name] = float(rf_list[new_risk]['normal_val'])
                attributes_list.append(attribute)
                flag = True
        return attributes_list, flag
#####################


#####################  functions related to updateMongoGeneralCitation #####################     
    ''' getNewDataCitations '''
    def getNewDataCitations(self):
        logging.debug(' Function:  getNewDataCitations')
        references = comprehensive_db['mod.n.mod'].find({'disease_id':self.disease_name,'genophen':'yes'},
                                                      {"_id":0,  'pmid':1,'author':1,'year':1,'title':1,
                                                       'journal':1, 'volume':1,'number':1,'pages':1,
                                                       'date-added':1}) 
        mle = MleUtil(config.DATAPATH, self.disease_name)
        new_disease_citations = {}
        for ref in references:
            if ref['pmid'] not in new_disease_citations:
                if ref['pmid']  != None:
                    new_disease_citations[ref['pmid']] = {"url": "http://www.ncbi.nlm.nih.gov//pubmed?term=%s" % ref['pmid']}
                    new_disease_citations[ref['pmid']]['author'] = mle.getAuthors(ref)
                    new_disease_citations[ref['pmid']]['content'] = mle.getContent(ref)
                else:
                    ref['pmid'] = ref['title']
                    new_disease_citations[ref['pmid']] = {"url": ref['title']}
                    new_disease_citations[ref['pmid']]['author'] =""
                    new_disease_citations[ref['pmid']]['content'] = ""
                
        return new_disease_citations        
##################### 


#####################  functions related to updateMongoRiskOutput ##################### 
    ''' getNewRiskOutput'''
    def getNewRiskOutput(self, diseases_list, flag):
        logging.debug(' Function:  getNewRiskOutput')
        new_diseases = []
        for doc in diseases_list:
            ### add risk factor in new disease
            if doc['disease'] == self.disease_name:
                doc["lfp"] = float(self.risk_genophen['lifetime_prevalence']['odds_ratio'])
                doc['risk_factors'] = {'genetics':{'f': 0, 'm': 0, 'v': 0}}
                doc['genetic_only'] = 'no'
                doc['paf'] = self.paf
                print doc
                for risk in self.risk_genophen:
                    risk_template ={'f': 0, 'm': 0, 'v': 0}
                    if risk != "ethnicity" and risk != "age"  and risk != "gender" and risk != "lifetime_prevalence" and risk != "incidence" and risk not in doc['risk_factors']:
                        if self.risk_genophen[risk]['risk_type'] == 'modifiable':
                            risk_template['m']=1                            
                        flag = True
                        doc['risk_factors'][risk] = risk_template
            new_diseases.append(doc) 
        return new_diseases, flag
#####################
    
##################### functions related to .bib files ###################
    ''' extract all the data of a bib file '''
    def extractDataFromBib(self):
        logging.debug(' Function:  extractDataFromBib')
        risk_data_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if i.startswith("%s_" % self.disease_name) and i.endswith('.bib')] 
        for fname in risk_data_list:
            risk_id = fname.split('%s_' % self.disease_name)[1].split('.bib')[0]
            self.util.warnMe('info', ' Risk factor: %s' % risk_id )
            bib = BibUtil(config.DATAPATH ,fname,self.disease_name, risk_id)
            self.risk_factor_data = bib.extractDataFromBib( self.risk_factor_data)
#####################


##################### function to update mongo to new data ###################   
    ''' add data to mongo comprehensive db / mod.nonmod collection '''
    def updateMongoComprehensive(self):
        logging.debug(' Function:  updateMongoComprehensive')
        ### output json
        self.util.writeJsonFile("comprehensive", self.risk_factor_data)
        ### clean up mod.n.mod for this disease
        count = comprehensive_db['mod.n.mod'].find({"disease_id":self.disease_name}).count()
        if  count>0:
            self.util.warnMe('info', ' deleting comprehensive/mod.n.mod data for %s - ndoc: %s' % (self.disease_name, count))    
            comprehensive_db['mod.n.mod'].remove({"disease_id":self.disease_name})
        ### add data in mod.n.mod for this disease
        self.risk_genophen = {}
        for doc in self.risk_factor_data:
            if doc['genophen']== 'yes':
                self.risk_genophen[doc['risk_id']] = doc
            doc["date_update"] = DATE
            doc["version"] = self.version
            comprehensive_db['mod.n.mod'].insert(doc)
        count = comprehensive_db['mod.n.mod'].find({"disease_id":self.disease_name}).count()
        self.util.warnMe('info', ' updated comprehensive/mod.n.mod data for %s - ndoc: %s' % (self.disease_name, count))  


    ''' add data to mongo algorithms collection '''
    def updateMongoAlgorithms(self):
        logging.debug(' Function:  updateMongoAlgorithms')
        ### get cooked data
        risk_factors = self.getRiskFactorAlgorithm()
        ### update disease_doc general information
        disease_doc={"date_update":DATE, "version": self.version, "disease" : self.disease_name, "risk_factors":risk_factors}
        ### write json + update/or not mongo
        self.recordNewMongoDocument("algorithms",disease_doc,True, genophen_db, 'algorithms',{"disease":self.disease_name} )

    
    ''' add data to mongo General/attributes collection '''
    def updateMongoGeneralAttributes(self):
        logging.debug(' Function:  updateMongoGeneralAttributes')
        
        ### get and record existing document
        attributes_doc = self.getRecordExistingMongoDocument(genophen_db,'general',{"name":"attributes"},{"_id":0}, 'attributes')

        ### initialize new attributes
        new_attributes, filename, flag_changes = self.initNewMongoDocument({}, 
                                                                  str(int(float(attributes_doc['version']))+0.1), 
                                                                  'attributes' )
        ### find which risk factor is new
        new_attributes['list'], template, new_rf_list = self.getExistingAttributes(attributes_doc) 
        
        ### add new risk factors  
        new_attributes['list'],flag_changes = self.getNewAttributes(new_attributes['list'], new_rf_list, template,  flag_changes)

        ### write json + update/or not mongo
        self.recordNewMongoDocument(filename,new_attributes,flag_changes, genophen_db, 'general',{"name":"attributes"} )

    
    ''' add data to mongo General/citation collection '''
    def updateMongoGeneralCitation(self):
        logging.debug(' Function:  updateMongoGeneralCitation')

        ### get list of citations for new disease
        new_disease_citations = self.getNewDataCitations()
        
        ### get and record existing document
        citations_doc = self.getRecordExistingMongoDocument(genophen_db,'general',{"name":"citation"},{"_id":0}, 'citation')

        ### initialize new citation doc
        citations_doc, filename, flag_changes = self.initNewMongoDocument(citations_doc, 
                                                                  str(int(float(citations_doc['version']))+0.1), 
                                                                  'citation' )
        ### add new disease in citation
        if self.disease_name not in  citations_doc['data']:
            citations_doc['data'][self.disease_name]={'medical':[],'genetics':[]}
        ### update citations if empty or asked for update
        if len(citations_doc['data'][self.disease_name]['medical']) == 0 or self.update_mongo: 
            citations_doc['data'][self.disease_name]['medical'] = []
            for pmid in new_disease_citations:
                citations_doc['data'][self.disease_name]['medical'].append(new_disease_citations[pmid])
        else:
            self.util.warnMe('info', ' No changes (updateMongoGeneralCitation)')
        ### write json + update/or not mongo
        self.recordNewMongoDocument(filename,citations_doc,True,genophen_db, 'general',{"name":"citation"} )
   

   
    ''' updateMongoRiskOutput template gid:0'''
    def updateMongoRiskOutput(self):
        logging.debug(' Function:  updateMongoRiskOutput')
         ### get and record existing document
        template = self.getRecordExistingMongoDocument(genophen_db,'riskoutput',{'gid':'0'},{"_id":0}, 'riskoutput')

        ### initialize new riskoutput doc
        template, filename, flag_changes = self.initNewMongoDocument(template, 
                                                                  str(int(float(template['version']))+0.1), 
                                                                  'riskoutput' )
        ### update list of diseases in riskoutput
        template['diseases'], flag_changes = self.getNewRiskOutput(template['diseases'], flag_changes)

        ### write json + update/or not mongo
        self.recordNewMongoDocument(filename,template,flag_changes,genophen_db, 'riskoutput',{"gid":"0"} )
         
    
    ''' main for class GetModNonModRiskFactors'''
    def main(self):
        logging.debug(' Function:  GetModNonModRiskFactors main %s' % self.disease_name)
        ### get data from .bib into mongo
        self.extractDataFromBib()
        ### update comprehensize
        self.updateMongoComprehensive()
        
        ### run R ###
        self.calculateModNMod()

        ### get R output into algorithm, general /"name": "attributes"riskoutput/"gid":"0"
        self.updateMongoAlgorithms()
        self.updateMongoGeneralAttributes()
        self.updateMongoGeneralCitation()
        self.updateMongoRiskOutput()
        logging.debug(' DONE with Function:  GetModNonModRiskFactors main %s' % self.disease_name)
        
        
""" MAIN function """    
if __name__ == "__main__":
    
    disease_name = "lung_cancer" 
    version = '3.0'
    update_mongo =False
    paf=   [
            0.0, 
            0.3, 
            0.6, 
            0.8, 
            0.85, 
            0.9, 
            0.99
        ]
    
    print '%s/%s%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'ModNmod_',disease_name)
    logging.basicConfig(filename='%s/%s%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'ModNmod_',disease_name), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    disease_risk = GetModNonModRiskFactors(disease_name, version, update_mongo, paf)
    disease_risk.main()
    
    print'DONE with DMmain'
""" END OF main """ 