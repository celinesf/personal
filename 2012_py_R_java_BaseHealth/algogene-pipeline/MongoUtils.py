#!/usr/bin/env python



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

import logging, json, copy
import datetime

import AlgoGeneConfig as config
from AlgoGeneUtil import AlgoGeneUtil
# from MleUtil import MleUtil
# from BibUtil import BibUtil
DATE = datetime.datetime.now().strftime("%m-%d-%Y")

class MongoUtils():
    def __init__(self, disease_name, genetic_data,version, update_mongo):
        self.disease_name = disease_name
        self.version = version
        self.update_mongo = update_mongo
        self.genetic_data = genetic_data
        self.ethmap = {"caucasian":0, "native_american":0, "hispanic":0, "asian":0,"african_american":0,'pacific_islander':0,'mixed':0}   
        self.util = AlgoGeneUtil()

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
            self.util.writeOutput("%s_%s_%s" %(name, date,doc['version']), doc)
            return doc

    ''' recordNewMongoDocument in json and mongo'''
    def recordNewMongoDocument(self, filename,new_doc,flag, db, collection,query ):
        logging.debug(' Function:  recordNewMongoDocument- filename %s, flag: %s' % (filename, flag))  
        if flag :
            ### output json
            self.util.writeOutput('%s_%s_old' %(self.disease_name,filename), new_doc)
            ### update algorithms
            if self.update_mongo:
                cursor = db[collection].find(query)
                nd = 0
                for doc in cursor:
                    nd+=1
                if nd ==1: 
                    self.util.warnMe('info', ' mongo %s/%s was updated (recordNewMongoDocument)' %(collection, query))
                    db[collection].update(query,new_doc)
                else: 
                    self.util.warnMe('info', ' mongo %s/%s was inserted (recordNewMongoDocument)' %(collection, query))
                    db[collection].insert(new_doc)
                ## check ok
                cursor = db[collection].find(query)
                nd = 0
                for doc in cursor:
                    nd+=1
                if nd !=1:
                    self.util.warnMe('error','NO DOCUMENT (recordNewMongoDocument)' %(collection, query))
            else:
                self.util.warnMe('warning', ' DID NOT UPDATE mongo %s/%s (recordNewMongoDocument)' %(collection, query))
        else :
            self.util.warnMe('info', ' No changes in %s/%s (recordNewMongoDocument)' %(collection, query))   
        ### output json
        self.util.writeOutput('%s_%s_new' %(self.disease_name,filename), new_doc)
    
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

#####################  functions related to updateMongoGeneralCitation #####################     
#     ''' getNewDataCitations '''
#     def getNewDataCitations(self):
#         logging.debug(' Function:  getNewDataCitations')
# 
#         mle = MleUtil(config.DATAPATH, self.disease_name)
#         new_disease_citations = {}
#         for ref in references:
#             if ref['pmid'] not in new_disease_citations:
#                 new_disease_citations[ref['pmid']] = {"url": "http://www.ncbi.nlm.nih.gov//pubmed?term=%s" % ref['pmid']}
#                 new_disease_citations[ref['pmid']]['author'] = mle.getAuthors(ref)
#                 new_disease_citations[ref['pmid']]['content'] = mle.getContent(ref)
#         return new_disease_citations        
##################### 


#####################  functions related to updateMongoRiskOutput ##################### 
    ''' getNewRiskOutput'''
    def getNewRiskOutput(self, diseases_list, flag):
        logging.debug(' Function:  getNewRiskOutput')
        new_diseases = []
        for doc in diseases_list:
            ### add risk factor in new disease
            if doc['disease'] == self.disease_name:
                risk_template ={'f': 0, 'm': 0, 'v': 0}
                doc['genetic_only'] = 'no'
                for risk in self.risk_genophen:
                    if risk != "incidence" and risk != "lifetime_prevalence" and risk != "intercepts" and risk not in doc['risk_factors']:
                        flag = True
                        doc['risk_factors'][risk] = risk_template
            new_diseases.append(doc) 
        return new_diseases, flag
#####################
    
##################### functions related to .bib files ###################
#     ''' extract all the data of a bib file '''
#     def extractDataFromBib(self):
#         logging.debug(' Function:  extractDataFromBib')
#         risk_data_list = [i for i in os.listdir('%s' % (config.DATAPATH)) if i.startswith("%s_" % self.disease_name) and i.endswith('.bib')] 
#         for fname in risk_data_list:
#             risk_id = fname.split('%s_' % self.disease_name)[1].split('.bib')[0]
#             self.util.warnMe('info', ' Risk factor: %s' % risk_id )
#             bib = BibUtil(config.DATAPATH ,fname,self.disease_name, risk_id)
#             self.risk_factor_data = bib.extractDataFromBib( self.risk_factor_data)
#####################
    
    ''' add data to mongo General/citation collection '''
    def updateMongoGeneralCitation(self):
        logging.debug(' Function:  updateMongoGeneralCitation')

        ### get list of citations for new disease
        new_disease_citations = self.getNewDataCitations()
        
        ### get and record existing document
        citations_doc = self.getRecordExistingMongoDocument(config.GENOPHEN_DB,'general',{"name":"citation"},{"_id":0}, 'citation')

        ### initialize new citation doc
        citations_doc, filename, flag_changes = self.initNewMongoDocument(citations_doc, 
                                                                  str(int(float(citations_doc['version']))+1.0), 
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
        self.recordNewMongoDocument(filename,citations_doc,True,config.GENOPHEN_DB, 'general',{"name":"citation"} )
   

    ''' getSnpListEthnicities'''
        ### count snp per ethnciity + add NA to studied in + list snp/ethncity per snps
    def getSnpListEthnicities(self):
        logging.debug(' Function:  getSnpListEthnicities')
        self.cooked_data = []
        self.snp_list ={}
        for doc in self.genetic_data:
            doc['gene'] = doc['gene'].replace('_',',')
            doc['gene'] = doc['gene'].replace('hCG,','hCG_')
            if doc['snp'] not in self.snp_list:
                self.snp_list[doc['snp']] = {}
            if doc['study_in'] not in self.snp_list[doc['snp']]:
                self.ethmap[doc['study_in']] += 1
                new_pop = 0
                self.snp_list[doc['snp']][doc['study_in']] = []
                self.snp_list[doc['snp']][doc['study_in']].append(doc)
            else:
                self.snp_list[doc['snp']][doc['study_in']].append(doc)
                new_pop +=1
            if new_pop == 2:
                new_doc = copy.deepcopy(doc)
                new_doc['gt'] = 'NA'
                new_doc['or'] = 0
                new_doc['oldor'] = -1
                new_doc['expor'] = -1
                self.snp_list[doc['snp']][doc['study_in']].append(new_doc)
                self.cooked_data.extend(self.snp_list[doc['snp']][doc['study_in']])
            elif new_pop >2:
                print 'issue'
        

    ''' findDominantEthnicity'''
    def findDominantEthnicity(self):
        logging.debug(' Function:  findDominantEthnicity')
        dominant = ""
        max_snp  = 0
        for eth in self.ethmap:
            if self.ethmap[eth] > max_snp:
                dominant = eth
        return dominant

    ''' getGeneticDataFormat'''
    def getGeneticDataFormat(self):
        logging.debug(' Function:  getGeneticDataFormat')
        self.getSnpListEthnicities()
        dominant = self.findDominantEthnicity()
        
        print 'dominant ethnicity',  dominant
        print 'snp_list', self.snp_list
        ### fiil missing ethnicies with  dominant data
        for snp in self.snp_list:
            if dominant in self.snp_list[snp]:
                for eth in self.ethmap:
                    if self.ethmap[eth] == 0:
                        for doc in self.snp_list[snp][dominant]:
                            new_doc = copy.deepcopy(doc)
                            new_doc['ethnicity'] = eth
                            self.cooked_data.append(new_doc)


    ''' updateGenetics '''
    def updateGenetics(self):
        logging.debug(' Function:  updateGenetics')
        ### get cooked data
        self.getGeneticDataFormat()

        ### update disease_genetics_doc general information
        self.disease_genetics_doc={"date":DATE, "version": self.version, "disease" : self.disease_name, "snps":self.cooked_data}
        ### write json + update/or not mongo
        self.recordNewMongoDocument("genetics.rep",self.disease_genetics_doc,True, config.GENOPHEN_DB, 'genetics.rep',{"disease":self.disease_name} )

    ''' main for class MongoUtils'''
    def main(self):
        logging.debug(' Function:  MongoUtils main %s' % self.disease_name)
        ### get data from .bib into mongo
        
        self.updateGenetics()
        
#         self.extractDataFromBib()
#         self.updateMongoGeneralCitation()
#         self.updateMongoRiskOutput()
        logging.debug(' DONE with Function:  MongoUtils main %s' % self.disease_name)
        
        
""" MAIN function """    
if __name__ == "__main__":
    data = config.GENOPHEN_DB['general'].find({'name':'gene_function'},{'_id':0})
    gene_list = []
    for doc in data:
        for gene in doc:
            if gene not in gene_list:
                gene_list.append(gene)
            else: 
                print 'found duplicate', gene
    print gene_list
    data = config.GENOPHEN_DB['genetics.rep'].find({},{'_id':0})
    new_genes = []
    for disease in data:
        for doc in  disease['snps']:
            if '_' in doc['gene'] and 'hcg_' not in doc['gene'].lower():
                print doc
            genes = doc['gene'].split(',')
            print disease['disease'], doc['gene']
            for g in genes:
                g = g.replace('.','_')
                if g not in gene_list and g not in new_genes:
                    new_genes.append(g)
    print new_genes
    f = open("%s/%s.json" % (config.OUTPUTPATH,'new_genes'), "w")
    f.write(json.dumps(new_genes, indent=4))
    f.close()
       
        
       
    print'DONE with DMmain'
""" END OF main """ 