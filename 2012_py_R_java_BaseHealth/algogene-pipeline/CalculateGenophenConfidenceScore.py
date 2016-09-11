#!/usr/bin/python

'''
Created on October, 29 2012
@author: celine
Version: 0.1

calculate genophen confidence score for snp a LD block
'''

import logging
import math
import scipy.stats

'''
haplotype_aggregator
'''
class CalculateGenophenConfidenceScore():
    ### __init__ ###
    ### Constructor
    def __init__(self):
        ''' '''
        self.genophen_score = {}
        #logging.debug(' Function: __init__ of CalculateGenophenConfidenceScore class' )
    ### END OF __init__

    ### CalculateGenophenConfidenceScore
    ###
    def calculateGenophenConfidenceScore(self,selected_variant):
        logging.debug(' ************ Function: calculateGenophenConfidenceScore ************** %s' % selected_variant )
        snp = {'snp': selected_variant["LD_best_variant"]['variant_id'], 'score': 0}
        ### replication
        rep = selected_variant['sign_sample'] - selected_variant['nsign_sample']
        key = 'LD_best_variant' 
        error = selected_variant[key]['ethnicity_score']['best_sample']['best_model']['max_error']
        OR = selected_variant[key]['ethnicity_score']['best_sample']['best_model']['min_OR']
        pvalue = selected_variant[key]['ethnicity_score']['best_sample']['best_model']['max_pvalue']

        ### errorsc
        if error != -1:
            sc = math.log(OR)/(error/(2 * 1.96)) # CI ~ 4*std
        elif pvalue != -1 :
            if pvalue == 0.00 : pvalue = 0.001
            Z = abs(scipy.stats.norm.ppf(pvalue/2))
            sample_size =  selected_variant[key]['ethnicity_score']['best_sample']['sample_size']      
            sc = math.log(OR)/(OR * math.sqrt(sample_size) / Z)
        
        snp['score'] = sc*math.sqrt(rep)        
        if snp['snp']  in self.genophen_score:
            self.genophen_score[snp['snp']] = min(snp['score'] , self.genophen_score[snp['snp']])
        else:
            self.genophen_score[snp['snp']] = snp['score']
    ### END OF CalculateGenophenConfidenceScore
