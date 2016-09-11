#!/usr/bin/python

'''
Created on October, 19 2012
@author: celine
Version: 0.1

Base/concrete Class haplotype_aggregator
- replace variants by haplotype if applicable
- temporary 3.1 will be replace by updated variant_selector to include haplotypes in the search
'''

import logging
import copy
import pymongo
import re

mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genetics"];

### instance from abstract class 
from GeneticEngine import abstract_haplotype_selector

'''
haplotype_aggregator
'''
class haplotype_selector(abstract_haplotype_selector):
    ### __init__ ###
    ### Constructor
    def __init__(self,disease_name,selected_variant):
        logging.debug(' Function: __init__ of haplotype_aggregator class, disease_name %s' % (disease_name))
        '''
        self.disease   = disease_name
        self.version = '0.1'
        self.selected_variant = selected_variant
        self.selected_variant_with_haplotype = {}
        '''
        ### call abstract data_aggregator
        super(haplotype_selector, self).__init__(disease_name,selected_variant)
    ### END OF __init__
    
    ### update_selected_variant_with_haplotype
    ###
    def get_haplotype_genophen_yes(self):
        logging.debug(' Function (get_haplotype_genophen_yes): disease: %s' % (self.disease))
        haplotype_data = {}
        finddoccursor = GENOPHEN30.snpor.aggregate([
                {'$match':{'$and':[
                       {'DiseaseID':self.disease}, 
                       {'dbSNP':{'$ne':'rs*'}},
                       {'Genophen':'KAP'}]}},
                {'$project':{'dbSNP':1,'GeneNextBio':1,'DiseaseID':1 ,'PMID':1,'BroadEthnicity':1,'_id':0 }},
                {'$group':{'_id':{'snp':'$dbSNP','gene':'$GeneNextBio', 'disease':'$DiseaseID','pmid':'$PMID','eth':'$BroadEthnicity','genophen':'$Genophen'}}}])   

        for info in finddoccursor['result']:
            eth = info['_id']['eth']
            snp = info['_id']['snp']
            pmid = info['_id']['pmid']
            gene = info['_id']['gene']
            try: 
                if haplotype_data[eth] :
                    'do nothing'
            except: 
                haplotype_data[eth] = {}
            haplotype_data[eth][snp] = {}
            haplotype_data[eth][snp]['gene'] = gene
            haplotype_data[eth][snp]['pmid'] = pmid
        return haplotype_data
    ### END OF update_selected_variant_with_haplotype
    
    ### add_new_ld_block
    ###
    def add_new_ld_block(self, haplotype, eth):
        logging.debug(' Function (add_new_ld_block)' )
        for hap_name in haplotype:
            logging.info(' Hap name: %s, gene block: %s (add_new_ld_block)' %(hap_name,haplotype[hap_name]['gene']) )
            print(' ADDING BLOCK Hap name: %s, gene block: %s (add_new_ld_block)' %(hap_name,haplotype[hap_name]['gene']) )
            gene = haplotype[hap_name]['gene']
            self.selected_variant_with_haplotype[eth][gene] = {}
            self.selected_variant_with_haplotype[eth][gene]['LD_best_variant'] = hap_name
            self.selected_variant_with_haplotype[eth][gene]['sample_best_OR'] = 'Haplotypic'
            self.selected_variant_with_haplotype[eth][gene]['variant_best_sample'] = haplotype[hap_name]['pmid']
            self.selected_variant_with_haplotype[eth][gene]['max_error'] = None
            self.selected_variant_with_haplotype[eth][gene]['max_pvalue'] = None
            self.selected_variant_with_haplotype[eth][gene]['significant_samples'] = None
            self.selected_variant_with_haplotype[eth][gene]['total_samples'] = None
    ### END OF add_new_ld_block
    
    ### update_selected_variant_with_haplotype
    ###
    def update_selected_variant_with_haplotype(self):
        super(haplotype_selector, self).update_selected_variant_with_haplotype()
        logging.debug(' ************ Function: haplotype_selector.update_selected_variant_with_haplotype -- VERSION: %s, disease: %s **************' % (self.version,self.disease ))
        self.selected_variant_with_haplotype = copy.deepcopy(self.selected_variant)
        haplotype_data = self.get_haplotype_genophen_yes()

        ### loop on eth in haplotypes
        for eth in haplotype_data:
            hap_info = haplotype_data[eth]
            ## loop of genes in selected variant for this ethnicity
            for gene in self.selected_variant[eth]:
                print ' ethnicity: %s, gene_block: %s' % (eth,gene) 
                hap_info_copy = copy.deepcopy(hap_info)
                ## loop on hap in hapinfo
                for hap_name in hap_info:
                    gene_hap = hap_info[hap_name]['gene'].split('_')
                    for g in gene_hap:
                        if g in gene.split('_'):
                            #if len(self.selected_variant_with_haplotype[eth][gene]['LD_best_variant'].split('rs')) > 1:
                            if re.match(re.compile('rs*'),self.selected_variant_with_haplotype[eth][gene]['LD_best_variant']) is not None:
                                logging.info(' ethnicity: %s, gene_block: %s, hap_name: %s, gene_hap: %s, gene:%s (update_selected_variant_with_haplotype) -- replacing variant by haplotype' % (eth, gene, hap_name,gene_hap,g))
                                self.selected_variant_with_haplotype[eth][gene]['LD_best_variant'] = hap_name
                                self.selected_variant_with_haplotype[eth][gene]['sample_best_OR'] = 'Haplotypic'
                                self.selected_variant_with_haplotype[eth][gene]['variant_best_sample'] = hap_info[hap_name]['pmid']
                                del hap_info_copy[hap_name] # remove this haplotype from the list
                                break
                            else:
                                logging.warning(' ethnicity: %s, gene_block: %s, hap_name: %s, gene_hap: %s, gene:%s ( (update_selected_variant_with_haplotype) -- I ALREADY PUT AN HAPLOTYPE IN THIS GENE BLOCK:%s, hap:%s' % (eth, gene, hap_name,gene_hap,g))
                                print (' PROBLEM ethnicity: %s, gene_block: %s, hap_name: %s, gene_hap: %s, gene:%s ( (update_selected_variant_with_haplotype) -- I ALREADY PUT AN HAPLOTYPE IN THIS GENE BLOCK:%s, hap:%s' % (eth, gene, hap_name,gene_hap,g))
                        else:
                            logging.info(' ethnicity: %s, gene_block: %s, hap_name: %s, gene_hap: %s, gene:%s ( (update_selected_variant_with_haplotype) -- do not replace' % (eth, gene, hap_name,gene_hap,g))
                hap_info = copy.deepcopy(hap_info_copy)
                if len(hap_info) == 0: 
                    break
            if len(hap_info) != 0:
                self.add_new_ld_block(hap_info, eth)
    ### END OF update_selected_variant_with_haplotype
