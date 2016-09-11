#!/usr/bin/python

'''
Created on October, 17 2012
@author: gwas_data
Version: 0.1

Class AggregateVariantData
former cauto with improvements
    - mongo db instead of sql
    - logs
    - record variation data has dictionary (no more output files)
'''

import logging
import copy
import AlgoGeneConfig as config
from AlgoGeneUtil import AlgoGeneUtil

###### AggregateVariantData class
class AggregateVariantData(object):
    ### __init__ ###
    ### Constructor
    def __init__(self, disease_name):
        logging.debug(' Function: __init__ of AggregateVariantData class -- disease: %s' % disease_name)
        self.version = '2.0'
        self.disease = disease_name
        self.ethmap =["ASW","CEU","JPT","CHB","CHD","GIH","LWK","MEX","MKK","TSI","YRI",'Global or other']
        self.unique_list = {}
        self.variant_data = {}   
        self.util = AlgoGeneUtil()
    ### END __init__
 

    ### snpHapmap ###
    ###  to see whether snp exists in hapmap and has frequency for the ethnicity
    def snpHapmap(self, snp, ethn):
        logging.debug(' Function: snpHapmap -- snp: %s, ethn: %s' % (snp,ethn))
        finddoccursor = config.NEXTBIO_DB.nextbio.cooked.find(
                {'disease_id':self.disease,'dbsnp':snp,'broad_ethnicity':ethn },
                {'_id':0})
        found = False
        for doc in finddoccursor:
            if doc['hapmap_pop'] in self.ethmap:
                if doc['hapmap_allele_freq'] is None:
                    self.util.warnMe('warning',' snp: %s, ethn: %s (snpHapmap) -- I did not find HapMap frequencies %s' % (snp,ethn,doc))
                    found = False
                    break
                else:
                    found = True
                    logging.info(' snp: %s, ethn: %s (snpHapmap) -- I got HapMap frequencies' % (snp,ethn))   
            else:                   
                self.util.warnMe('warning',' snp: %s, ethn: %s (snpHapmap) -- ETHNICITY NOT FOUND %s' % (snp,ethn,doc))
                found = False
                break
        return found
#     ### END OF snp_hapmap

    ### isBiAllelicSnp ###
    def isBiAllelicSnp(self, data):
        logging.debug(' Function: isBiAllelicSnp -- bioset: %s, snp: %s' % (self.bioset, self.snp))
        if data['is_indel'] == True:
            return False
        elif len( data['other_allele']) != 2:
            return False
        else:
            return True

    ### getSampleModelData ###
    ### fill bioset general data + record model raw data
    def getSampleModelData(self):
        logging.debug(' Function: getSampleModelData -- bioset: %s, snp: %s' % (self.bioset, self.snp))
        model_data = {} 
        self.sample_id = '%s_%s' %(self.snp,self.bioset)
        study_type = 0
        row_num = 0 
       
        finddoccursor = config.NEXTBIO_DB.nextbio.cooked.find(
        {'disease_id':self.disease,'dbsnp':self.snp,'bioset_id':self.bioset },
        {'_id':0})
        ### Loop on data
        for or_data in finddoccursor:
            ### define study_type
            if row_num == 0:
                study_type = self.defineStudyType(or_data)
                is_snp = self.isBiAllelicSnp(or_data)
                sample_data = {"sample_id" : self.sample_id, 
                      "sample_size" : float(or_data['sample_size']), 
                      "study_type":study_type, 
                      "is_snp":is_snp,
                      "meta_sample":or_data['in_meta'],
                      "association": {}}          
            model_data = self.defineModelType(or_data,model_data)
            row_num += 1
        return sample_data, model_data

    ### getCi ###
    def getCi(self, data):
        logging.debug(' Function: getCi -- bioset: %s, snp: %s' % (self.bioset, self.snp))
        if data['ci']:
            return  [ float(data['ci'].split("-")[0]) ,   float(data['ci'].split("-")[1]) ]
        else:
            return  [-1, -1]  

    ### getOR ###
    def getOR(self, data):
        logging.debug(' Function: getOR -- bioset: %s, snp: %s' % (self.bioset, self.snp))
        if str(data['or_value']) == '1.0':
            return 1
        elif data['or_value'] is None:
            return -1 
        else:
            return float(data['or_value'])
    ### getEffectSize ###
    def getEffectSize(self, data):
        if str(data['effect_size']) == 0:
            return 1
        elif data['effect_size'] is None:
            return -1 
        else:
            return float(data['effect_size'])

    ### getPvalue ###
    def getPvalue(self, data):
        logging.debug(' Function: getPvalue -- bioset: %s, snp: %s' % (self.bioset, self.snp))
        if data['pvalue'] == None:
            return -1
        else:
            return float(data['pvalue'])   

    ### fillOrData ###
    def fillOrData(self, sample_data, model_data, template_eth):
        logging.debug(' Function: fillOrData -- bioset: %s, snp: %s' % (self.bioset, self.snp))
        ### loop on models
        for model_id in model_data:                    
            model = '%s_%s_%s' %(self.snp,self.bioset,model_id)          
            rows = [row for row in model_data[model_id]['data']]
            sample_data['association'][model] ={"model_id":model, "control_freq": model_data[model_id]['control_freq'], "OR_data":{} }           
            template_eth["sample"][self.sample_id] = sample_data
            ## loop on OR rows
            for row in rows:
                ## get CI
                ci = self.getCi(row)
                ## reference has OR=1 (integer)     
                OR = self.getOR(row) 
#                 if OR is None:
                          
                ## pvalue as -1 if NULL   
                pvalue =self.getPvalue(row) 
                ## make sure the is allele
                if  row['allele']:
                    allele_genotype = row['allele']
                    template_eth['sample'][self.sample_id]["association"][model]["OR_data"][allele_genotype ] = {
                                        "OR_id" : str(allele_genotype),
                                        "odd_ratio":OR, "CI" : ci,
                                        "pvalue" : pvalue}
                else:
                    self.util.warnMe('critical', ' NO ALLELE (fillOrData) - bioset: %s, snp: %s, data: %s' % (self.bioset, self.snp, row))
            self.finalize(self.snp, self.ethnicity, template_eth)
            return template_eth
       
    ### getSnpPaperSample ###
    ###  for a given snp, pubmedid find the sample  
    def getSnpPaperSample(self, snp):
        logging.debug(' Function: getSnpPaperSample -- snp: %s' % snp)
        
        ### init output data structure
        self.variant_data[snp] =  {'variant_id':snp,  "nextbio_score":-1,  "nextbio_pvalue":-1, "broad_ethnicity" : { }}
        ### loop on ethnicities
        for eth in self.unique_list[snp]:
            ### init template for this ethnicity
            eth_dict = {"ethnicity_id" : eth, 'hapmap_freq':self.snpHapmap(snp, eth), "sample" : {} }
            self.ethnicity = eth
            for bioset in self.unique_list[snp][eth]:
                self.bioset = bioset
                sample_data, model_data = self.getSampleModelData()
                eth_dict = self.fillOrData(sample_data, model_data, copy.deepcopy(eth_dict))

    ### END OF getSnpPaperSample
 
    ### defineModelType ###
    ### Called by getSnpPaperSample                       
    def  defineModelType(self,or_data,model_data):
        logging.debug(' Function: defineModelType -- bioset:%s, snp:%s' % (self.bioset,self.snp))
        if len(or_data['allele']) == 1:
            model_data = self.addAllelicGenotypic(model_data,'allelic', or_data)
        elif len(or_data['allele']) == 2:
            model_data = self.addAllelicGenotypic(model_data,'genotypic', or_data)
        else:
            logging.warning(' bioset:%s, snp:%s (defineModelType) -- I GOT NOR Allelic NOR Genotypic -- allele:%s' % (self.bioset,self.snp,or_data['allele']))      
        return model_data
    ### END OF defineModelType
  
    ''' defineStudyType 
    define meta 4, gwas discovery 1, replication 2
    ### Called by getSnpPaperSample           '''            
    def  defineStudyType(self,or_data):
        logging.debug(' Function: defineStudyType -- bioset:%s, snp:%s' % (self.bioset,self.snp))
        if or_data['is_meta'].lower() == 'yes':
            return 4
        elif or_data['data_type'].lower() == 'discovery':
            return 1
        elif or_data['data_type'].lower() == 'replication':
            return 2
        elif or_data['data_type'].lower() == 'combined':
            return 3
        else:
            logging.warning(' bioset:%s, snp:%s (getSnpPaperSample) -- CANNOT DEFINE STUDY TYPE' %(self.bioset,self.snp))
            return None        
    ### END OF defineStudyType
       
    ### initializeData ###
    ### Called by addAllelicGenotypic                       
    def  initializeData(self,data,model):
        logging.debug(' Function: initializeData -- bioset:%s, snp:%s, model: %s' % (self.bioset,self.snp, model))
        data[model] = {}
        data[model]['control_freq'] = False
        data[model]['data'] = []
        return data
    ### END OF initializeData
       
    '''
    ### addAllelicGenotypic ###
    - add data to allelic or genotypic model
    - if new model, initialize data
    ### called by getSnpPaperSample '''
    def addAllelicGenotypic(self, data,model,orrow):
        logging.debug(' Function: addAllelicGenotypic -- bioset:%s, snp:%s, model: %s' % (self.bioset,self.snp, model))
        try:
            data[model]['data'].append(orrow)
            if (orrow['frequency_pop'] or orrow['frequency_control']) and data[model]['control_freq'] == False:
                logging.info(' model: %s, allele_genotype: %s (addAllelicGenotypic) -- found control frequencies' % (model,orrow['allele']))
                data[model]['control_freq'] = True
            else: 
                logging.warning(' bioset:%s, snp:%s, model: %s, allele_genotype: %s (addAllelicGenotypic) -- CANNOT FIND control frequencies' % (self.bioset,self.snp, model,orrow['allele']))                  
        except:
            data = self.initializeData(data,model)
            data = self.addAllelicGenotypic( data,model,orrow)
        return data
    ### END OF addAllelicGenotypic
       
    ### finalize ###
    ###    
    def finalize(self, snp, ethnicity, eth_dic):
        logging.debug(' Function: finalize -- snp: %s, ethnicity: %s' % (snp, ethnicity))
        self.variant_data[snp]["broad_ethnicity"][ethnicity] = eth_dic
    ### END OF finalize
       
    ### getGwasData ###
    ###     
    def getGwasData(self):
        logging.debug(' Function: getGwasData' )
        ns = 0
        for snp in self.unique_list:
            self.snp = snp
            self.getSnpPaperSample(str(snp))
            ns += 1
            if ns %50 ==0:
                print ns
    ### END OF getGwasData
    

    ### getUniqueLists ###
    ###  obtain unique SNPS , ethnicity per snp, bioset per ethnicity per snp, of the disease
    def getUniqueLists(self):
        logging.debug(' Function: getUniqueLists' )
        finddoccursor = config.NEXTBIO_DB.nextbio.cooked.find({'disease_id':self.disease},
                                                              {'dbsnp':1,'bioset_id':1,'broad_ethnicity':1,'_id':0})
        for doc in finddoccursor:
            if doc['dbsnp'] not in self.unique_list:
                self.unique_list[doc['dbsnp']] = {}
            if doc['broad_ethnicity'] not in self.unique_list[doc['dbsnp'] ] :
                self.unique_list[doc['dbsnp']][doc['broad_ethnicity']] = []
            if doc['bioset_id'] not in self.unique_list[doc['dbsnp'] ][doc['broad_ethnicity']] :
                self.unique_list[doc['dbsnp']][doc['broad_ethnicity']].append(doc['bioset_id'] )
    ### END OF getUniqueSnps
       
    ################# MAIN ###############
    def main(self):
        logging.debug(' **************** Function: main AggregateVariantData -- VERSION: %s ****************' % self.version )
        self.getUniqueLists()
        self.getGwasData()
        
    ### END OF main

