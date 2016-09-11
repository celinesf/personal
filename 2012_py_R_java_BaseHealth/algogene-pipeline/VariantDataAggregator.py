#!/usr/bin/python

'''
Created on October, 17 2012
@author: celine
Version: 0.1

Base/concrete Class VariantDataAggregator -> aggerateVariantData
1) former cauto with improvements
    - mongo db instead of sql
    - logs
    - record variation data has dictionary (no more output files)
2) former LDVariantSelectionInputGenerator -> aggerateLdBlockData
    - sql to mongo query
    - record ld_blok information data has dictionary (no more output files)
'''

import logging

### instance from abstract class 
from GeneticEngine import AbstractVariantDataAggregator
from AggregateVariantData import AggregateVariantData
from AggregateLdBlockData import AggregateLdBlockData

'''
VariantDataAggregator
'''
class VariantDataAggregator(AbstractVariantDataAggregator):
    ### __init__ ###
    ### Constructor
    def __init__(self,disease_name):
        logging.debug(' Function: __init__ of VariantDataAggregator class, disease_name %s' % (disease_name))
        '''
        self.disease   = disease_name
        self.version = '0.1'
        self.aggregated_variant_data = {}
        self.aggregated_ld_block_data  = {}
        '''
        ### call abstract VariantDataAggregator
        super(VariantDataAggregator, self).__init__(disease_name)
    ### END OF __init__
    
    ### aggerateVariantData
    def aggerateVariantData(self):
        super(VariantDataAggregator, self).aggerateVariantData()
        logging.debug(' Function: VariantDataAggregator.aggerateVariantData -- VERSION: %s, disease: %s' % (self.version,self.disease ))

        variant_data = AggregateVariantData(self.disease)
        variant_data.main()
        self.aggregated_variant_data = variant_data.variant_data 
    ### END OF aggerateVariantData
     
    ### aggerateLdBlockData
    def aggerateLdBlockData(self):
        super(VariantDataAggregator, self).aggerateLdBlockData()
        logging.debug(' Function: VariantDataAggregator.aggerateLdBlockData -- VERSION: %s, disease: %s' % (self.version,self.disease ))

        ld_block_data = AggregateLdBlockData(self.disease)
        ld_block_data.main()
        self.aggregated_ld_block_data = ld_block_data.gene_block 
    ### END OF aggerateVariantData   
    


    
