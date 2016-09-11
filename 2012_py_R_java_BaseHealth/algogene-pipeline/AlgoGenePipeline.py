#!/usr/bin/python

'''
Created on October, 17 2012
@author: celine

Algorithm & Flow

Genetic engine consists of these major phases
I)   variant aggregation : all the information is gathered for the genetic variation of interest
II)  Data aggregation    : all the variants are aggregated for the disease of interest
III) Data selection      : good variants are selected for the disease of interest (QC & selection routine)
IV)  Adjusted calculation: selected variatns go throug adjusted odds ration calculation
V)   Ethnicity decision  : Genophen covers about 6 different ethinicities at the time of this writting. what is the default ethnicity for the disease of interest

and log properly

main function to call
1) data aggregator class
2) variant selector class
'''


import sys
import logging
import AlgoGeneConfig as config
from GeneticEngine import GeneticEngine
from VariantDataAggregator import VariantDataAggregator
from VariantSelector import VariantSelector
from VariantOrAdjuster import VariantOrAdjuster


################# MAIN ###############
if __name__ == '__main__':

#     disease_name = 'atrial_fibrillation'
#     lifetime_prevalence = 0.166 #
  
#     disease_name = 'gallstone'
#     lifetime_prevalence = 0.1

#     disease_name = 'allergic_rhinitis'
#     lifetime_prevalence = 0.15
    
#     disease_name = 'gout'
#     lifetime_prevalence = 0.026
     
#     disease_name = 'venous_thrombosis'
#     lifetime_prevalence = 0.05 #
    
#     disease_name = 'kidney_stone'
#     lifetime_prevalence = 0.05
    
#     disease_name = 'alcohol_dependence'
#     lifetime_prevalence = 0.125 #

#     disease_name = 'melanoma'
#     lifetime_prevalence = 0.02 #

#     disease_name = 'age_related_macular_degeneration'
#     lifetime_prevalence = 0.06 #
     
#     disease_name = 'rheumatoid_arthritis'
#     lifetime_prevalence = 0.05 #
      
    
    if len(sys.argv) > 1:
        disease_name = sys.argv[1]
    if len(sys.argv) > 2:
        lifetime_prevalence = sys.argv[2]

    logging.basicConfig(filename='%s/%s_%s' %(config.OUTPUTPATH, config.LOG_FILENAME,disease_name), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')

    GeneticEngine = GeneticEngine(disease_name,lifetime_prevalence,VariantDataAggregator,VariantSelector,VariantOrAdjuster)
    GeneticEngine.modelGeneticDisease()
    
    print "I AM DONE with the genetic engine for disease" , disease_name
### END OF MAIN
    
    
    