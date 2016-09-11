#!/usr/bin/env python
"""
   Genetic Engine Dashboard
   09/15/12 - 0.0.1 : 
"""


"""
   Algorithm & Flow
   
   Genetic engine consists of these major phases
   I)   variant aggregation : all the information is gathered for the genetic variation of interest
   II)  Data aggregation    : all the variants are aggregated for the disease of interest
   III) Data selection      : good variants are selected for the disease of interest (QC & selection routine)
   IV)  Adjusted calculation: selected variatns go throug adjusted odds ration calculation
   V)   Ethnicity decision  : Genophen covers about 6 different ethinicities at the time of this writting. what is the default ethnicity for the disease of interest
   
   and log properly
"""

import logging
import abc
import json
import sys
import AlgoGeneConfig as config
from MongoUtils import MongoUtils
from AlgoGeneUtil import AlgoGeneUtil

'''
GeneticEngine class
- calls the concreate classes
- generate the selected_variant needed for each step
'''
class GeneticEngine():
    ### __init__ ###
    ### Constructor    
#     def __init__(self,disease_name, VariantDataAggregator, VariantSelector, haplotype_selector, haplotype_or_calculator, AbstractVariantOrAdjuster): #, adj_ethn):
    def __init__(self,disease_name,lifetime_prevalence, VariantDataAggregator, VariantSelector,VariantOrAdjuster): #, adj_ethn):
        logging.debug(' Function: __init__ of GeneticEngine, disease_name %s' % (disease_name))
        self.disease = disease_name
        self.lifetime_prevalence = lifetime_prevalence
        self.VariantDataAggregator = VariantDataAggregator
        self.VariantSelector = VariantSelector
        self.VariantOrAdjuster = VariantOrAdjuster
        logging.debug(' Function: __init__ of GeneticEngine, disease_name %s' % (disease_name))
        self.util = AlgoGeneUtil()

    ### modelGeneticDisease ###
    # pipeline for genetic selected_variant for a disease
    ###   
    def modelGeneticDisease(self):
        logging.debug(' Function: modelGeneticDisease of GeneticEngine, disease_name %s' % (self.disease))
        ### aggregate variant and LD block selected_variant
        aggregated_data = self.VariantDataAggregator(self.disease)
   
        ### output variant OR selected_variant
        self.util.writeOutput('%s_variants' % self.disease, aggregated_data.aggregated_variant_data)
        ### output LD/gene block selected_variant
        self.util.writeOutput('%s_LD' % self.disease, aggregated_data.aggregated_ld_block_data)  
              
        ## QC, select best variant for each ld block
#         aggregated_variant_data = self.util.fecth_input_data('%s/%s_variants.json'%(config.OUTPUTPATH,self.disease),'variants' )
#         aggregated_ld_block_data = self.util.fecth_input_data('%s/%s_LD.json'%(config.OUTPUTPATH,self.disease),'LD' )
#         selected_variant = self.VariantSelector(self.disease,aggregated_variant_data, aggregated_ld_block_data)
     
        selected_variant = self.VariantSelector(self.disease,aggregated_data.aggregated_variant_data, aggregated_data.aggregated_ld_block_data)
         
        ### calculated adjusted OR
#         selected_variant = self.util.fecth_input_data('%s/%s_selection.json'%(config.OUTPUTPATH,self.disease),'selection' )
#         print selected_variant
#         variant_ajusted_or = self.VariantOrAdjuster(self.disease,self.lifetime_prevalence, selected_variant) 
 
        variant_ajusted_or = self.VariantOrAdjuster(self.disease,self.lifetime_prevalence, selected_variant.selected_variant) 
 
        self.util.writeOutput('%s_adjOR' % self.disease, variant_ajusted_or.variant_adjusted_or)
        self.util.writeOutput('%s_pmid' % self.disease, variant_ajusted_or.pmid)
    
#         variant_adjusted_or = self.util.fecth_input_data('%s/%s_adjOR.json'%(config.OUTPUTPATH,self.disease),'adjOR' )
#         mongoUtils = MongoUtils(self.disease, variant_adjusted_or,'1.0' ,True)
  
        mongoUtils = MongoUtils(self.disease, variant_ajusted_or.variant_adjusted_or,'2.0' ,False)
        mongoUtils.main()
          
########################### HAPLOTYPE ONLY TO ADD LATER ###################
#         ### replaces selected_variant with haplotypes in applicable ld_blocks
#         selected_variant_with_haplotype = self.haplotype_selector(self.disease,selected_variant.selected_variant)
#         print selected_variant.selected_variant
#      
#         ### replaces selected_variant with haplotypes in applicable ld_blocks
#         haplotype_adjusted_or = self.haplotype_or_calculator(self.disease,selected_variant_with_haplotype.selected_variant_with_haplotype)
#         print json.dumps(haplotype_adjusted_or.haplotype_adjusted_or)

        ### update mongodb genetics.rep
#         mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.142.17/admin"
#         CW = pymongo.Connection(mongo_cnxn_string, 27017);
#         GENOPHEN30 = CW["genetics"]
#         GENOPHEN30.genuimap.remove({'disease':self.disease})
#         for hap in haplotype_adjusted_or.haplotype_adjusted_or:
#             GENOPHEN30.genuimap.insert(hap)


        print "I AM DONE with disease" , self.disease
    ### END OF GeneticEngine


''' - aggregates the genetic selected_variant for all the varint studied in a disease
    - aggregates selected_variant specific to each gene/LD block
'''
class AbstractVariantDataAggregator():
    __metaclass__ = abc.ABCMeta
    ### __init__ ###
    ### Constructor    
    def __init__(self,disease_name):
        logging.debug(' Function: __init__ of AbstractVariantDataAggregator, disease_name %s' % (disease_name))
        self.disease   = disease_name
        self.version = '2.0'
        
        ### initiation
        self.aggregated_variant_data = {}
        self.aggregated_ld_block_data  = {}
        
        ### call AbstractVariantDataAggregator.aggerateVariantData
        print "SARTING variant aggregation disease" , disease_name
        self.aggerateVariantData()
        print "I AM DONE variant aggregation disease" , disease_name
         
        ### call AbstractVariantDataAggregator.aggerate_ld_block_datat
        print "SARTING generating LD blocks for disease" , disease_name
        self.aggerateLdBlockData()
        print "I AM DONE generating LD blocks for disease" , disease_name
    ### END OF __init__
        
    @abc.abstractmethod
    def aggerateVariantData(self):
        logging.info(' Function: AbstractVariantDataAggregator.aggerateVariantData -- VERSION: %s, disease: %s, selected_variant: %s' % (self.version,self.disease,self.aggregated_variant_data ))
        ''' aggregates the genetic selected_variant for all the variant studied in a disease'''
 
    @abc.abstractmethod
    def aggerateLdBlockData(self):
        logging.info(' Function: AbstractVariantDataAggregator.aggerateVariantData -- VERSION: %s, disease: %s' % (self.version,self.disease ))
        ''' aggregates snps within each gene block for a disease'''


''' select best variant per LD block and best OR per variant '''
class AbstractVariantSelector():
    __metaclass__ = abc.ABCMeta
    ### __init__ ###
    ### Constructor 
    def __init__(self, disease_name, variant_data,ld_block_data):
        logging.debug(' Function: AbstractVariantSelector, disease_name %s' % (disease_name))
        self.disease   = disease_name
        self.version = '2.0'
        self.variant_data = variant_data
        self.ld_block_data = ld_block_data
        self.selected_variant = {}

        ### call AbstractVariantSelector.selectVariant
        print "STARTING WITH SELECTING VARIANT FOR disease" , disease_name
        self.selectVariant()
        print "I AM DONE WITH SELECTING VARIANT FOR disease" , disease_name
    ### END OF __init__
    
    @abc.abstractmethod
    def selectVariant(self):
        logging.info(' Function: AbstractVariantSelector.selectVariant -- VERSION: %s, disease: %s' % (self.version,self.disease ))
        ''' Select the best variant for each LD block for disease_name '''
    
''' calculate adjusted OR '''    
class AbstractVariantOrAdjuster():
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, disease_name,lifetime_prevalence, selected_variant):
        self.disease = disease_name
        self.lifetime_prevalence = lifetime_prevalence
        self.selected_variant     = selected_variant
        self.variant_adjusted_or = []
        
        # function calls
        self.calculateAdjustedOr()

    @abc.abstractmethod
    def calculateAdjustedOr(self):
        ''' calculate adjusted OR for variants '''
        


  
        

