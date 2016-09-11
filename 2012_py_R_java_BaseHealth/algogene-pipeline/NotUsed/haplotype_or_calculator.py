#!/usr/bin/python

'''
Created on October, 22 2012
@author: celine

- recover haplotype data in snpor
- calculate genotypic OR if needed
- calculate adjusted OR using R module
'''

import json
import math
import copy
import pymongo
import logging
import re
import os
import sys

#mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.142.17/admin"
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genetics"]

### instance from abstract class 
from GeneticEngine import abstract_haplotype_or_calculator


class haplotype_or_calculator(abstract_haplotype_or_calculator):
    ''' this class is responsible for calculation adjusted OR for Haps '''
    ### __init__ ###
    ### Constructor
    def __init__(self, disease_name,selected_variant_with_haplotype):
        logging.debug(' Function: __init__ of haplotype_or_calculator class -- disease: %s' % disease_name)
        '''
        self.disease   = disease_name
        self.version = '0.1'
        self.haplotype_adjusted_or = []
        self.haplotype_adjusted_or['snps'] = []
        self.haplotype_adjusted_or['disease'] = disease_name
        self.selected_variant_with_haplotype = selected_variant_with_haplotype
        '''
        self.dest= "" # path
        self.defor = 1.0 #if OR is not significant set it to this
        self.json_genot  = self.dest + str(disease_name) + '_Hap_Geno.json'
        self.json_out    = self.dest + str(disease_name) + '_Hap_Geno_out.json'
        self.haplotype_data = {}
        self.haplotypic_genotype = {}
        self.json_allel  = self.dest + str(disease_name) + '_Hap_Allele.json'
        
        ### call abstract data_aggregator
        super(haplotype_or_calculator, self).__init__(disease_name,selected_variant_with_haplotype)
    ### END OF __init__
    
    ### get_unique_haplotypes ###
    # get haplotypes from snpor
    # called from main
    ### 
    def get_unique_haplotypes(self):
        logging.debug(' Function: get_unique_haplotypes')
        for ethn in self.selected_variant_with_haplotype:
            for gene in self.selected_variant_with_haplotype[ethn]:
                if re.match(re.compile('rs*') , self.selected_variant_with_haplotype[ethn][gene]['LD_best_variant']) == None:
                    if self.selected_variant_with_haplotype[ethn][gene]['LD_best_variant'] not in self.haplotype_data:
                        self.haplotype_data[self.selected_variant_with_haplotype[ethn][gene]['LD_best_variant']] = {}
    ### END OF get_unique_haplotypes   

    ### get_genotype_from_data ###
    # get genotype of haplotypes directly from snpor
    # called by get_haplotype_data_from_db
    ### 
    def get_genotype_from_data(self, data):
        logging.debug(' Function: get_genotype_from_data -- hap: %s, genotype: %s' % (data['dbSNP'], data['AG_Plus']))

        hapname = data['dbSNP']
        eth = data['BroadEthnicity']
        genotype = data['AG_Plus']
        OR = data['ORValue']
        CI = self.get_ci(data)       
        orflag = self.check_ci_allele(CI)
        if orflag == False:
            CIissue = 'No'
        else:
            CIissue = 'Yes'
            OR = 1.0
        
        if hapname not in self.haplotypic_genotype: # check hapnam exists
            self.haplotypic_genotype[hapname] = {}
        if eth not in self.haplotypic_genotype[hapname]: # check ethnicity already exists
            self.haplotypic_genotype[hapname][eth] = {}
        if genotype in self.haplotypic_genotype[hapname]: # check if allele already exist- it should not
            logging.warning(' hapname: %s allele: %s (get_genotype_from_data) -- THIS ALLELE DOES NOT EXIST?' %(hapname,genotype))
        else:
            logging.info(' hapname: %s allele: %s (get_genotype_from_data) -- adding new genotype' %(hapname,genotype))
            self.haplotypic_genotype[hapname][eth][genotype] = {}
            self.haplotypic_genotype[hapname][eth][genotype]["OR"] =  OR
            self.haplotypic_genotype[hapname][eth][genotype]["FCase"] = data['FreqCase']
            self.haplotypic_genotype[hapname][eth][genotype]["FControl"] = data['FreqControl']
            self.haplotypic_genotype[hapname][eth][genotype]["CIissue"] = CIissue
            
            self.haplotypic_genotype[hapname][eth][genotype]["CI"] = CI 
            #self.haplotypic_genotype[hapname][eth][genotype]["CaseNo"] = 2*c['CaseNo'] * c['FreqCase']
            #self.haplotypic_genotype[hapname][eth][genotype]["ControlNo"] = 2*c['ControlNo'] *c['FreqControl']
    ### END OF get_genotype_from_data  
        
    ### get_haplotype_data_from_db ###
    # get haplotypes from snpor
    # called from main
    ###    
    def get_ci(self, data):
        logging.debug(' Function: get_ci')
        ci = data['CI']
        OR= data['ORValue']
        pval = data['Pvalue']
        CI = []
        if ci == None and OR != 1:
            if pval < 0.05:
                logging.info(' (get_ci) -- pval is non-significant')
                CI = [1.0001,1.0002]  
            else:
                logging.info(' (get_ci) -- pval is significant')
                CI = [0.9999,1.0001]  
        elif ci == None and OR == 1:
            logging.info(' (get_ci) -- I am a reference')
            CI = [-1,-1]  
        else:
            logging.info(' (get_ci) -- geting CI')
            for c in ci.split('-'):
                CI.append(float(c))
        return CI
    ### END OF get_ci                         

    ### get_haplotype_data_from_db ###
    # get haplotypes from snpor
    # called from main
    ###    
    def get_haplotype_data_from_db(self):
        logging.debug(' Function: get_haplotype_data_from_db')
        self.get_unique_haplotypes()
        
        for hapname in self.haplotype_data:
            finddoccursor = GENOPHEN30.snpor.aggregate(
               [{'$match':{'$and':
                           [{'DiseaseID':self.disease},
                            {'dbSNP':hapname},
                            {'Genophen': 'KAP'},
                            {'AG_Plus': {'$ne': None}},
                            {'FreqCase': {'$ne': None}},
                            {'FreqControl': {'$ne': None}}]}}])    
            for c in finddoccursor['result']:
                ### record genotypes
                eth = c['BroadEthnicity']
                allele = c['AG_Plus']
                if len(allele.split('/')) == 2:
                    self.get_genotype_from_data(c)
                ### record alleles
                else:
                    if hapname not in self.haplotype_data: # check hapnam exists
                        print 'SHOULD NEVER BE HERE'
                        self.haplotype_data[hapname] = {}
                    if eth not in self.haplotype_data[hapname]: # check ethnicity already exists
                        self.haplotype_data[hapname][eth] = {}
                    if allele in self.haplotype_data[hapname]: # check if allele already exist- it should not
                        logging.warning(' hapname: %s allele: %s (get_haplotype_data_from_db) -- THIS ALLELE DOES NOT EXIST?' %(hapname,allele))
                    else:
                        logging.info(' hapname: %s allele: %s (get_haplotype_data_from_db) -- adding new allele' %(hapname,allele))
                        self.haplotype_data[hapname][eth][allele] = {}
                        self.haplotype_data[hapname][eth][allele]["CI"] = self.get_ci(c)
                        self.haplotype_data[hapname][eth][allele]["or"] = c['ORValue'] 
                        self.haplotype_data[hapname][eth][allele]["FCase"] = c['FreqCase']
                        self.haplotype_data[hapname][eth][allele]["FControl"] = c['FreqControl']
                        ### TO FIX
                        self.haplotype_data[hapname][eth][allele]["CaseNo"] = round(2*c['CaseNo'] * c['FreqCase'], 10)
                        self.haplotype_data[hapname][eth][allele]["ControlNo"] = round(2*c['ControlNo'] *c['FreqControl'], 10)
                        ### TO REMOVE
#                        self.haplotype_data[hapname][eth][allele]["CaseNo"] = round(1*c['CaseNo'] * c['FreqCase'], 10)
#                        self.haplotype_data[hapname][eth][allele]["ControlNo"] = round(1*c['ControlNo'] *c['FreqControl'], 10)
    ### END OF get_haplotype_data_from_db            
            
    ### calculate_or_ci ###
    # caluclate OR and CI for haplotype alleles
    # called from main
    ###    
    def calculate_or_ci(self):
        logging.debug(' Function: calculate_or_ci')
        ''' calculate or from frequencies in case and control '''
        for hap in self.haplotype_data.keys():
            if len(self.haplotype_data[hap]) >0:
                logging.info(' hap: %s (calculate_or_ci) -- need to calculate allelic or' %(hap))
                for ethn in self.haplotype_data[hap].keys():
                    majallele = self.major_allele(hap, ethn) ### find major allele or reference if already have the reference
                    majinfo   = self.haplotype_data[hap][ethn][majallele]
                    for allele in self.haplotype_data[hap][ethn].keys():
                        if "OR" not in self.haplotype_data[hap][ethn][allele]:
                            logging.info(' hap: %s, allele: %s (calculate_or_ci) -- new reference' %(hap, allele))
                            a_info = self.haplotype_data[hap][ethn][allele]
                            a_or   = round(( majinfo["FControl"] / majinfo["FCase"] ) * ( a_info["FCase"] / a_info["FControl"] ),10)
                            ci_2   = math.exp(math.log(a_or) + 1.96 *  math.sqrt(1.0/majinfo["CaseNo"] + 1.0/majinfo["ControlNo"] + 1.0/a_info["CaseNo"] + 1.0/a_info["ControlNo"] ) )
                            ci_1   = math.exp(math.log(a_or) - 1.96 *  math.sqrt(1.0/majinfo["CaseNo"] + 1.0/majinfo["ControlNo"] + 1.0/a_info["CaseNo"] + 1.0/a_info["ControlNo"] ) ) 
                            self.haplotype_data[hap][ethn][allele]["OR"] = a_or
                            self.haplotype_data[hap][ethn][allele]["CI"] = [ci_1, ci_2]
                            if a_or == 1.0:
                                self.haplotype_data[hap][ethn][allele]["CI"] = [-1,-1]
                        else:
                            logging.info(' hap: %s, allele: %s (calculate_or_ci) -- use old OR' %(hap, allele))
            else:
                logging.info(' hap: %s (calculate_or_ci) -- I do not need to calculate allelic or' %(hap))
        ### END OF calculate_or_ci 

    ### calculate_haplotype_genotypic_OR ###
    # caluclate OR and CI for haplotype genotypes
    # called from main
    ###   
    def calculate_haplotype_genotypic_OR(self):
        logging.debug(' Function: calculate_haplotype_genotypic_OR')
        '''allele combinations'''
        
        ### loop on haplotype names
        for hap in self.haplotype_data.keys():
            if len(self.haplotype_data[hap]) >0 :
                logging.info(' hap %s(calculate_haplotype_genotypic_OR) -- need to calculate genotype or' %(hap))
                ## copy strcture of haps object
                self.haplotypic_genotype[hap] = copy.deepcopy(self.haplotype_data[hap])
                ### loop on ethnicities
                for ethn in self.haplotype_data[hap].keys():
                    genotype = {}
                    alleles  = self.haplotype_data[hap][ethn].keys()
    
                    for i in range(len(alleles)):
                        for j in range(i, len(alleles)):
                            logging.info(' hap: %s, ethn: %s, hap1: %s, hap2: %s (calculate_haplotype_genotypic_OR) ' % (hap, ethn,alleles[i],alleles[j]))
                            ### set OR to 1 if confidence interval includes 1
                            if self.check_ci_allele(self.haplotype_data[hap][ethn][alleles[i]]["CI"]) == True:
                                or1     = self.defor
                                orflag1 = True
                            else:
                                or1     = self.haplotype_data[hap][ethn][alleles[i]]["OR"]
                                orflag1 = False
                            if self.check_ci_allele(self.haplotype_data[hap][ethn][alleles[j]]["CI"]) == True:
                                or2     = self.defor
                                orflag2 = True
                            else:
                                or2     = self.haplotype_data[hap][ethn][alleles[j]]["OR"]
                                orflag2 = False
                            ### decide of genotype significance
                            logging.info(' hap: %s, ethn: %s, hap1: %s, hap2: %s (calculate_haplotype_genotypic_OR) -- OR1: %s, OR2: %s' % (hap, ethn,alleles[i],alleles[j],orflag1,orflag2))
                         
                            if orflag2 == True or orflag1 == True:
                                logging.warning(' hap: %s, ethn: %s, hap1: %s, hap2: %s (calculate_haplotype_genotypic_OR) -- ISSUE' % (hap, ethn,alleles[i],alleles[j]))
                                orflag = "Yes"
                                ### DO NOT set genotype OR to 1 as we don't know if genotype is significant or not
    #                            or1 = self.defor
    #                            or2 = self.defor
                            else:
                                logging.warning(' hap: %s, ethn: %s, hap1: %s, hap2: %s (calculate_haplotype_genotypic_OR) -- OK' % (hap, ethn,alleles[i],alleles[j]))
                                orflag = "No"
                            num = 2
                            if i == j:
                                num = 1
                            fcase = round(num *self.haplotype_data[hap][ethn][alleles[i]]["FCase"] * self.haplotype_data[hap][ethn][alleles[j]]["FCase"], 10)
                            fcontr= round(num *self.haplotype_data[hap][ethn][alleles[i]]["FControl"] * self.haplotype_data[hap][ethn][alleles[j]]["FControl"], 10)
                            genotype[alleles[i]+"/"+alleles[j]] = {
                                                                   "FCase"    : fcase, 
                                                                   "FControl" : fcontr,
                                                                   "OR"       : or1 * or2,
                                                                   "CIissue"   : orflag
                                                                   }
                    self.haplotypic_genotype[hap][ethn] = genotype
            else:
                logging.info(' hap: %s (calculate_haplotype_genotypic_OR) -- I already have genotype data' %(hap))
    ### END OF get_haplotype_genotype_or  

    ### major_allele ###
    # decide is major allele
    # called by  calculate_or_ci
    def major_allele(self, hap, ethn):
        logging.debug(' Function: major_allele')
        majallele_freq = 0
        reference_allele = None
        N_allele = None
        for allele,val in self.haplotype_data[hap][ethn].items():
            if val["FControl"] > majallele_freq :# and val['FControl'] >= val["FCase"]: ### TO FIX
                logging.info(' hap: %s (major_allele) -- I found major allele: %s' %(hap,allele))
                majallele_freq = val["FControl"]
                majallele      = allele  
            if re.match(re.compile('N*'),allele):
                N_allele = allele
            if val["or"] == 1 and  val["CI"][0] == -1:
                logging.info(' hap: %s (major_allele) -- I found reference allele: %s' %(hap,allele))
                reference_allele = allele
        ## TO FIX
        if reference_allele != None and N_allele is None:
            logging.info(' hap: %s (major_allele) -- I keep OR from SNPOR' %(hap))
            for allele in self.haplotype_data[hap][ethn]:
                self.haplotype_data[hap][ethn][allele]['OR'] = self.haplotype_data[hap][ethn][allele]['or']
        return majallele
    ### END OF major_allele 
  
    ### check_ci_allele ###
    # if CI includes one -> true
    # called by calculate_haplotype_genotypic_OR
    def check_ci_allele(self, ci):
        logging.debug(' Function: check_ci_allele -- ci %s %s %s' % (ci, ci[0], ci[1]))
        l = ci[0]
        h = ci[1]
        if l == h and h== -1:
            logging.info(' (check_ci_allele) -- I am reference allele')
            return False # its the ref allele
        if l<= 1 <=h:
            logging.info(' (check_ci_allele) -- I am NOT significant')
            return True
        else:
            logging.info(' (check_ci_allele) -- I am significant')
            return False
    ### END OF check_ci_allele 

    ### adjusted_or_post_processing
    #   post processing for Database insertion #
    # after non-linear equation solving is done in R, push things to db
    # Called from main
    def adjusted_or_post_processing(self):
        logging.debug(' Function: adjusted_or_post_processing')
        ### read R output of adjusted OR
        json_data  = open(self.json_out)
        jsn        = json.load(json_data)
        json_data.close()
#         os.system('rm %s ' % self.json_out) ### TO FIX

        ### create genuimap like documents for this haplotype   
        for hap in jsn.keys():
            hap_data = self.get_hap_info(hap)
            #genename = self.snpOrGenename(hap)
            for ethn in jsn[hap].keys():
                for gt in jsn[hap][ethn].keys():
                    newdoc = copy.deepcopy(hap_data)
                    gtinfo  = jsn[hap][ethn][gt]
                    newdoc['gt'] = gt
                    newdoc['gt'] = gt
                    newdoc['freq'] = gtinfo["FControl"]
                    newdoc['oldor'] = gtinfo["OR"]
                    newdoc['or'] = gtinfo["NewOR"]
                    newdoc['lnor'] = math.log(gtinfo["NewOR"])
                    newdoc['ethnicity'] = ethn
                    newdoc['study_in'] = ethn
                    newdoc['disease'] = self.disease
                    self.haplotype_adjusted_or.append(newdoc)
    ### END OF adjusted_or_post_processing 
    
    ### get_hap_info
    #   find chr, pos, gene name from hapID
    def get_hap_info(self,hap):
        logging.debug(' Function: get_hap_info -- hap: %s' % hap)
        hap_data = {}
        finddoccursor = GENOPHEN30.IDHap.rep.aggregate([
            {'$match':{'dbHAP':hap}},
            {'$project':{'dbSNP':1,
                         'GeneNextBio':1,
                         'dbHAP':1,
                         'position':1,
                         'position36':1,
                         'chromosome':1,
                         'relative_position':1,
                         'place':1,
                         '_id':0
                         }},
            {'$group':{'_id':{'hap':'$dbHAP','snp':'$dbSNP','gene':'$GeneNextBio','position':'$position','position36':'$position36','chromosome':'$chromosome', 'relative_position':'$relative_position','place':'$place'}}}])     
        for c in finddoccursor['result']:
            doc = c['_id']
            if len(hap_data) == 0:
                gene = doc['gene'].split('_')
                hap_data['snp'] = hap
                hap_data['gene'] = gene
                hap_data['chromosome'] = doc['chromosome']
            ### find position of haplotype
            if doc['place'] == 1:
                hap_data['position'] = doc['position']
                hap_data['position36'] = doc['position36']
                hap_data['relative_position'] = doc['relative_position']
            ### check chromosome number OK (should be
            if doc['chromosome'] != hap_data['chromosome']:
                logging.debug(' hap: %s, chr: %s (find_hap_info) -- ERROR CHROMOSOME FOR DBSNP: %s (chr: %s)' % (hap,hap_data['chromosome'],doc['snp'],doc['chromosome']))
                sys.exit(' hap: %s, chr: %s (find_hap_info) -- ERROR CHROMOSOME FOR DBSNP: %s (chr: %s) %s' % (hap,hap_data['chromosome'],doc['snp'],doc['chromosome'], doc))
            ### extend gene name if needed
            gene = doc['gene'].split('_')
            for g in gene:
                if g not in hap_data['gene']:
                    hap_data['gene'].append(g)
        ### cat gene names into a sctring
        count = 1
        for g in hap_data['gene']:
            if count == 1:
                gene = g 
            else:
                gene = '%s,%s' % (gene, g)
            count += 1
        hap_data['gene'] = gene
        return hap_data
    ### END OF get_hap_info

    ### calculate_adjusted_OR_with_R ###
    # calculate adjusted using R
    # called by main
    def calculate_adjusted_OR_with_R(self):
        logging.debug(' Function: calculate_adjusted_OR_with_R' )
        cmdline = 'Rscript nonlinear_solver.R %d > Routput' % ( self.disease)
        os.system(cmdline)
        Rfile = open('Routput','r')
        logging.info('\n%s' % Rfile.read()) 
        Rfile.close()
        ## remove input and out for R 
#         os.system('rm Routput %s' % self.json_genot) ### TO FIX
    ### END OF calculate_adjusted_OR_with_R
    
    ### calculate_haplotype_adjusted_or ###
    # calculate adjusted using R
    # called by main
    def calculate_haplotype_adjusted_or(self):
        super(haplotype_or_calculator, self).calculate_haplotype_adjusted_or()
        logging.debug(' ************ Function: haplotype_or_calculator.calculate_haplotype_adjusted_or -- VERSION: %s, disease: %s **************' % (self.version,self.disease ))

        ### get data for haplotypes
        self.get_haplotype_data_from_db()
        
        if len(self.haplotype_data) > 0 :
#            #### TO REMOVE OUTPUT
            f = open(self.json_allel, "w")
            f.write(json.dumps(self.haplotype_data, sort_keys=True, indent=4))
            f.close()
#            ### END TO REMOVE
            
            ## calculate_or_ci if allele - otherwise copy if genotype already
            self.calculate_or_ci()
            
            ## calculate_haplotype_genotypic_OR
            self.calculate_haplotype_genotypic_OR()
            print 
            ### prepare input file for R script that calculate adjusted OR
            f = open(self.json_genot, 'w')
            f.write(json.dumps(self.haplotypic_genotype, sort_keys=True, indent=4))
            f.write('\n')
            f.close()
            ### Run R module here
            self.calculate_adjusted_OR_with_R()
    
            ###final step : push to GenUITable
            self.adjusted_or_post_processing()
    ### END OF calculate_haplotype_adjusted_or
