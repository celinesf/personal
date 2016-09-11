#!/usr/bin/python

'''
Created on October, 17 2012
@author: celine
Version: 0.1

Base/concrete Class VariantSelector
selectVariant is Equivalent to the main function of the formally called Genophen Variant Selection.
'''

import logging
import AlgoGeneConfig as config

### instance from abstract class 
from GeneticEngine import AbstractVariantSelector

### variant scoring class functions
import ScoreVariants as SC
ScoreVariants = SC.ScoreVariants()

### LD variant selection class
import SelectVariantPerLdBlock as LD
SelectVariantPerLdBlock = LD.SelectVariantPerLdBlock()

### shared functions and variables
import VariantSelectorSharedFunctions as SF
STAT = SF.STATISTICS
SharedFunctions = SF.VariantSelectorSharedFunctions()

from AlgoGeneUtil import AlgoGeneUtil

class VariantSelector(AbstractVariantSelector):
    
    def __init__(self,disease_name,  variant_data,ld_block_data):
        logging.debug(' Function: VariantSelector, disease_name %s' % (disease_name))
        '''
        self.disease   = disease_name
        self.version = '0.1'
        self.variant_data = variant_data
        self.ld_block_data = ld_block_data
        self.selected_variant = {}
        '''
        self.util = AlgoGeneUtil()
        ### call abstract VariantSelector
        super(VariantSelector, self).__init__(disease_name, variant_data,ld_block_data)
    ### END OF init
   
    ### reformatDataForScoring ###
    # filters out input for scoreVariants and selectVariantPerLdBlock
    # reformat into dictionaries (i.e., no array)
    # called by selectVariant()
    ### 
    def reformatDataForScoring(self,indata):
        logging.debug(' Function: reformatDataForScoring')
    
        ### Initialization 
        input_data = {} # gwas data -> scoreVariants
        input_info = {} # variant specific scores -> selectVariantPerLdBlock   
        
        STAT['variant_scoring'] = {}
        self.initializeSummaryStatistics(STAT['variant_scoring'], 0)
    
        ### Loop on variant documents/dictionaries
        for rs in indata:
            input_data[rs] = {}
            input_info[rs] = {}
            STAT['variant_scoring']['sum_variant'] += 1
            
            for key in indata[rs]:
                value = indata[rs][key]
                ## Records variant ID 
                if key == 'variant_id': #key.split("_")[1] == 'id':
                    input_info[rs][key] = value
                
                ## Records nextbio score and pvalue
#                if key.split("_")[0] == 'nextbio':
                if  'nextbio' in key:
                    input_info[rs][key] = value
                
                ##  records document ethnicity  
#                elif key.split("_")[1] == 'ethnicity':
                elif  'ethnicity' in key:
                    input_data[rs][key] = value
                    input_info[rs][key] = {}
                    for key2 in indata[rs][key]:
                        STAT['variant_scoring']['total_ethnicity'] += 1
                        
                        ## update list of ethnicities
                        if key2 not in STAT['variant_scoring']['different_ethnicity']:
                            STAT['variant_scoring']['different_ethnicity'].append(key2)
    
                            ### INIT summary per ethnicity
                            STAT['variant_scoring'][key2] = {}
                            self.initializeSummaryStatistics(STAT['variant_scoring'][key2], 1 )
                        
                        ## record hapfreq in info dict for LD selection  
                        hapfreq = SharedFunctions.convertToBoolean(indata[rs][key][key2]['hapmap_freq'])
                        input_info[rs][key][key2] = hapfreq
                        
        return input_data, input_info
    ### END OF reformatDataForScoring
    
    ### initializeSummaryStatistics ###
    # Initialize the summary statistics dictionary
    ### 
    def initializeSummaryStatistics(self,stat, num):
        logging.debug(' Function: initializeSummaryStatistics')
        if num == 0:
            stat['different_ethnicity'] =[] 
            stat['total_ethnicity'] = 0 
            stat['total_sign_ethnicity'] = 0 
            stat['total_nsign_ethnicity'] = 0 
            stat['total_num_ethnicity'] = 0
            stat['sum_variant'] = 0 
            stat['sum_risk_allele'] = 0 
            stat['max_sample_size'] = 0
            stat['WARNING'] = {}
            stat['ERROR'] = {}
    
        stat['total_variant'] = 0 
        stat['total_risk_allele'] = 0 
        
        ## sample info
        stat['total_sample'] = 0 # IN
        stat['sum_sample_size'] = 0
        stat['sum_sample_type'] = 0
        stat['total_sign_sample'] = 0 
        stat['total_nsign_sample'] = 0 
        stat['total_num_sample'] = 0 
        
        ## model info
        stat['total_model'] = 0 
        stat['total_sign_model'] = 0 
        stat['total_nsign_model'] = 0 
        stat['total_OK_OR_model'] = 0 
        stat['total_NOTOK_OR_model'] = 0 
        stat['total_OK_model'] = 0 
        stat['total_NOTOK_model'] = 0 
    ### END OF initializeSummaryStatistics
    
    ### calculateSummaryStatistics ###
    # clean dictionaries for output purposes
    ### 
    def calculateSummaryStatistics(self):
        logging.debug(' Function: calculateSummaryStatistics')
        for et in STAT['variant_scoring']['different_ethnicity']:
            for key in STAT['variant_scoring'][et]:
                STAT['variant_scoring'][key] += STAT['variant_scoring'][et][key]
            if et in STAT['LD_selection']:
                for key in STAT['LD_selection'][et]:
                    STAT['LD_selection'][key] += STAT['LD_selection'][et][key]
    ### END OF calculateSummaryStatistics
          

        
    ### MAIN of Genophen Variant Selector
    def selectVariant(self):
        super(VariantSelector, self).selectVariant()
        logging.debug(' ************** Function: VariantSelector.selectVariant -- VERSION: %s, disease: %s **************' % (self.version,self.disease ))
#         baseline = '%s/%s' %(config.OUTPUTPATH, self.disease)
       
        ### filters the input data into input for scoreVariants and selectVariantPerLdBlock
        (variant_data, variant_info) = self.reformatDataForScoring(self.variant_data)
# 
#         ### score variants
        variant_score = ScoreVariants.scoreVariants(variant_data)
#         variant_score = self.util.fecth_input_data('%s/%s_score.json'%(config.OUTPUTPATH,self.disease),'score' )
            
        ### select best variants for each LD block
        self.selected_variant, self.genophen_score = SelectVariantPerLdBlock.selectVariantPerLdBlock(variant_score,variant_info,self.ld_block_data)
        ### output summary statistics
        self.calculateSummaryStatistics()
        self.util.writeOutput('%s_summary' % self.disease, STAT)
#         SharedFunctions.output_data(STAT, baseline+'_summary.json','summaries' )  
    
#        ########## TEST AND DEBUG ONLY #########
        ### intermediate output of variant scoring
        self.util.writeOutput('%s_score' % self.disease, variant_score)
        self.util.writeOutput('%s_selection' % self.disease, self.selected_variant)
        self.util.writeOutput('%s_working' % self.disease, SF.WORKING)
#         SharedFunctions.output_data(variant_score, baseline+'_score.json','variant score' )   
# 
#         ### out best variant per LD block per ethnicity
#         SharedFunctions.output_data(self.selected_variant, baseline+'_selection.json','LD variant selection' )   
#         SharedFunctions.output_data(SF.WORKING, baseline+'_working.json','working dictionary' ) 
 
    ### END OF selectVariant __main__

