#!/usr/bin/python

'''
Created on Jul 30, 2012
logs added 8/8/12
@author: celine
'''

import copy
import sys
import json
import logging

### functions and global variable from VariantSelectorSharedFunctions
import VariantSelectorSharedFunctions as SF
STAT = SF.STATISTICS
SharedFunctions = SF.VariantSelectorSharedFunctions()
import AlgoGeneConfig as config

class ScoreVariants(object):
    
    def __init__(self):
        '''
        Constructor
        '''
        self.ethnicity = ''
        self.variant_id = ''
    
    ### scoreVariants ###
    # choose OR for that sample
    # compare to best sample for an ethnicity
    # called by main
    ### 
    def scoreVariants(self,data):
        logging.info('\n================================\nFunction: scoreVariants\n================================')
        
        ### initialization
        variant_score = {} 
        ns = 0
        ### loop on variants -> choose best sample per ethnicity
        for rs in data: 
            logging.info(" variant: %s (scoreVariants)" % rs)
            
            ## variant_score initiation
            variant_score[rs] = self.initializeVariantScore()   
            variant_score[rs]['variant_id'] = rs
            self.variant_id = rs
               
            ### loop on ethnicity -> choose best sample
            for et in data[rs]['broad_ethnicity']: 
                logging.debug(" LOOP: ethnicity, scoreVariants %s, ethnicity: %s" % (rs, et))
            
                ## Init sumary for this ethnicity
                self.ethnicity = et
                STAT['variant_scoring'][self.ethnicity]['total_variant'] += 1
      
                ## check ethnicity ID
                SharedFunctions.checkIdFit(et,data[rs]['broad_ethnicity'][et]['ethnicity_id'])
                
                ## ethnicity_score initiation
                ethnicity_score = self.initializeEthnicityScore(data[rs]['broad_ethnicity'][et])          
                
                ### loop on samples -> score samples
                for sa in ethnicity_score['sample_score']: 
                    logging.info(" ***** variant %s, ethnicity: %s, sample %s (scoreVariants) -- starting sample_scoring *****" % (rs, et, sa))
                    STAT['variant_scoring'][self.ethnicity]['total_sample'] += 1
                    STAT['variant_scoring'][self.ethnicity]['sum_sample_size'] += ethnicity_score['sample_score'][sa]['sample_size']
                    STAT['variant_scoring'][self.ethnicity]['sum_sample_type'] += ethnicity_score['sample_score'][sa]['study_type']
                    STAT['variant_scoring']['max_sample_size'] = max(STAT['variant_scoring']['max_sample_size'],ethnicity_score['sample_score'][sa]['sample_size'])
 
                    ## check sample ID
                    SharedFunctions.checkIdFit(sa,ethnicity_score['sample_score'][sa]['sample_id'])
                    
                    ## score sample
                    sample_score = ethnicity_score['sample_score'][sa]
                    sample_score = self.scoreSample(sample_score,ethnicity_score['hapmap_freq'])
                
                    ## update statistics for this ethnicity
                    ethnicity_score = self.updateEthnicityScore(ethnicity_score,sample_score)
                        
                    ## compare sample to best sample so far
                    ethnicity_score['best_sample'] = SharedFunctions.compareSamples(ethnicity_score['best_sample'], sample_score,ethnicity_score['hapmap_freq'], True)

                ## update meta-analysis (if sample is a subsample of a meta-analysis)
                for sa in ethnicity_score['sample_score']: 
                    ethnicity_score = self.updateMetaAnalysisScore(ethnicity_score, ethnicity_score['sample_score'][sa]) # update meta analysis score
                
                ### Check Best sample has OK ethnicity risk allele
                ethnicity_score = self.checkRiskAlleleForBestSample(ethnicity_score)
                
                ### update variant scores
                variant_score[rs] = self.updateVariantScore(variant_score[rs] , ethnicity_score)

            STAT['variant_scoring']['sum_risk_allele'] += len(variant_score[rs]['variant_ethnicity_risk_allele'])
            ns += 1
            if ns%100 == 0:
                print ns , rs
        logging.info(" ***** I AM DONE SCORING THE self.variant_idS *****")
        return variant_score
    ### END OF scoreVariants

######################### Initialization functions #######################      
    ### initializeVariantScore ###
    # initialize output of varant_scoring function
    # called by scoreVariants
    ###
    def initializeVariantScore(self):        
        logging.debug(' Function: initializeVariantScore')
        variant = {}
        variant['variant_sign_ethnicity'] = 0
        variant['variant_nsign_ethnicity'] = 0
        variant['variant_num_ethnicity'] = 0
        variant['variant_ethnicity_risk_allele'] = {}
        variant['variant_total_risk_allele'] = {}
        variant['ethnicity_score'] = {}

        return variant
    ### END OF initializeVariantScore
    
    ### initializeEthnicityScore ###
    # copy all ethnicity data
    # add all intermediate output/ cooked data for this ethnicity
    # called by scoreVariants
    ###
    def initializeEthnicityScore(self,ethnicity):
        logging.debug(' Function: initializeEthnicityScore: ethnicity %s'% ( ethnicity['ethnicity_id']))
        et_score = copy.deepcopy(ethnicity)
        et_score['mean_study_type'] = 0
        et_score['mean_sample_size'] = 0
        et_score['mean_OR'] = 0
        et_score['mean_error'] = 0
        et_score['mean_pvalue'] = 0
        et_score['ethnicity_sign_sample'] = 0
        et_score['ethnicity_nsign_sample'] = 0
        et_score['ethnicity_num_sample'] = 0   
        et_score['ethnicity_risk_allele'] = {}
        et_score['sample_score'] = (ethnicity['sample'])
        
        ### loop on samples
        for sa in et_score['sample_score']: 
            sa_score = et_score['sample_score'][sa]
            sa_score['sample_risk_allele'] =  {}
            
            ## if meta analysis only
            sa_score['sign_sample'] =  0
            sa_score['nsign_sample'] =  0
            
            ### loop on models
            for mod in sa_score['association']: 
                mod_score = sa_score['association'][mod]
                mod_score['model_type'] = ""
                mod_score['significant'] = False
                mod_score['max_error'] = -1
                mod_score['max_pvalue'] = -1
                mod_score['min_OR'] = -1
                mod_score['risk_allele'] = ""
                mod_score['protective_allele'] = ""
                mod_score['reference_allele'] = ""
                mod_score['is_OK'] = False # true if all ORs have no error
                mod_score['OR_OK'] = False # true if all ORs have no error
                
                ### loop on allele/genotypes
                for ag in mod_score['OR_data']: 
                    ag_score = mod_score['OR_data'][ag]
                    ag_score['is_risk'] = False
                    ag_score['is_prot'] = False
                    ag_score['is_ref'] = False
                    ag_score['is_sign'] = False
                    ag_score['is_OK'] = False # true if no error
                    
            ## first model = best model so far for this sample       
            sa_score['best_model'] =  (sa_score['association'][sa_score['association'].keys()[0]])  
             
        ## first sample = best sample so far for this ethnicity 
        et_score['best_sample'] = (et_score['sample_score'][et_score['sample_score'].keys()[0]]) 
        
        return et_score
    ### END OF initializeEthnicityScore
######################### 
    
    ### scoreSample ###
    # finds the best model of association for this sample
    # Called by scoreVariants
    ### 
    def scoreSample(self, sample, hapfreq):
        logging.debug(' Function: scoreSample: sample %s' % sample['sample_id'])

        ### loop on model of association
        for mo in sample['association']: 
            logging.info(" sample: %s, model: %s (scoreSample) -- starting model_scoring" % (sample['sample_id'],mo))
            STAT['variant_scoring'][self.ethnicity]['total_model'] += 1
            
            ## check model ID
            SharedFunctions.checkIdFit(mo, sample['association'][mo]['model_id'])
            model = sample['association'][mo]

            ## assess quality of model
            model = self.scoreModel(model) 
            
            ## check if is a real snp
            if model['is_OK'] == True and sample['is_snp'] == False:
                logging.warning(" NOT BIALLELIC SNP sample: %s, model: %s (scoreSample) " % (sample['sample_id'],mo))
                model['is_OK'] = False
                
            ## find best model for this sample
            sample['best_model'] = SharedFunctions.compareModels(sample['best_model'], model, hapfreq, True)

        ### update sample_risk_allele 
        risk_allele = sample['best_model']['risk_allele'] 
        if sample['best_model']['is_OK'] == True and risk_allele != "":         
            ## if new allele
            if risk_allele in sample['sample_risk_allele']:
                logging.info(" sample: %s, risk_allele %s (scoreSample) -- updating sample_risk_allele count" % (sample['sample_id'],risk_allele))
                sample['sample_risk_allele'][risk_allele] += 1
                            
            ## if first allele    
            else  :
                logging.info(" sample: %s, risk_allele %s (scoreSample) -- creating sample_risk_allele" % (sample['sample_id'],risk_allele))
                sample['sample_risk_allele'][risk_allele] = 1   
     
        return sample
    ### END OF sample_scoring
    
    ### scoreModel ###
    # score a specific model
    # - define significant, minOR, maxpval and max CI
    # - define risk, protective and reference alleles 
    # Called by sample_scoring
    ###     
    def scoreModel(self,model_data):
        logging.debug(' Function: scoreModel: %s' % model_data['model_id'])
        # Initialization of number of OR
        nOR = 0
        ### Loop on alleles and genotypes
        for all_geno in model_data['OR_data']:
            logging.info(" model: %s, allele/genotype: %s (scoreModel) -- starting scoreOddRatio" % (model_data['model_id'],all_geno))
            nOR += 1
            
            ## check OR ID
            SharedFunctions.checkIdFit(all_geno, model_data['OR_data'][all_geno]['OR_id'])
            
            OR = model_data['OR_data'][all_geno]
            
            ## assess significant, check no error
            OR = self.scoreOddRatio(OR)  
            
            ## update model statistics  
            model_data = self.updateModelScore(model_data, OR)   

        ### report reference allele is risk or protective allele
        model_data = self.updateModelType(model_data,nOR)
        model_data = self.updateModelSignificance(model_data)
        return model_data
    ### END OF model_scoring
    
    ### scoreOddRatio ###
    # score a specific model
    # - define if significant
    # - define risk, protective and reference alleles 
    # Called by model_scoring
    ###     
    def scoreOddRatio(self,OR_data):
        logging.debug(' Function: scoreOddRatio: %s' % OR_data['OR_id'])
        
        ### check type is number
        OR_data['CI'][1]+1
        OR_data['odd_ratio'] +1
        OR_data['pvalue'] +1
        
        ### reference allele of genotype
        if OR_data['CI'][0] == -1 \
            and OR_data['CI'][1] == -1 \
            and OR_data['pvalue'] == -1 \
            and OR_data['odd_ratio'] == 1.0 : 
            logging.info(" allele/genotype: %s (scoreOddRatio) -- I am reference " % OR_data['OR_id'])
            OR_data['is_ref'] = True
            OR_data['is_OK'] = True
            
        ### allele/genotype with reported OR
        else:
            logging.debug(" allele/genotype: %s (scoreOddRatio) -- I have OR data" % OR_data['OR_id'])
            
            ## CI provided -> define risk/protective allele and sign from CI
            if OR_data['CI'][0] < OR_data['CI'][1]:    
                ## OR data found
                if OR_data['odd_ratio'] > 0 :   
                    ## OR is within CI range -> no error, sign or not
                    if OR_data['CI'][1] > OR_data['odd_ratio'] and OR_data['CI'][0] < OR_data['odd_ratio']:   
                        logging.info(" allele/genotype: %s (scoreOddRatio) -- I have full OR data" % OR_data['OR_id'] )              
                        OR_data = self.getSignificanceFromOrAndCi(OR_data)
                        
                    ## OR is outside CI range -> ERROR
                    else:
                        logging.error(" allele/genotype: %s (scoreOddRatio) -- CI & OR DO NOT FIT (OR: %s, CI: %s)" % (OR_data['OR_id'],OR_data['odd_ratio'],OR_data['CI']))
                        SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],'(scoreOddRatio) -- CI & OR DO NOT FIT')
                            
                ## no OR given
                else:     
                    logging.warning(" allele/genotype: %s (scoreOddRatio) -- CANNOT FIND ODD RATIO BUT I FOUND CI (OR: %s, CI: %s)" % (OR_data['OR_id'],OR_data['odd_ratio'],OR_data['CI']))
                    SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(scoreOddRatio) -- CANNOT FIND ODD RATIO BUT I FOUND CI')
                    OR_data = self.getSignificanceFromCiOnly(OR_data)  
                
            ## CI not provided -> define risk/protective allele and sign from pvalue  (TODO CHECK NEXTBIO DATA SHOULD ALWAYS BE SIGNIFICANT?)            
            elif OR_data['CI'][0] == -1 and OR_data['pvalue'] != -1:       
                logging.info(" allele/genotype: %s (scoreOddRatio) -- confidence measured by pvalue only ") 
                OR_data = self.getSignificanceFromPvalueOnly(OR_data) 
            
            ## No OR or CI           
            elif OR_data['odd_ratio'] != -1:       
                logging.warning(" allele/genotype: %s (scoreOddRatio) -- CANNOT TELL IF I AM SIGNIFICANT. MAYBE OK? ") 
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(scoreOddRatio) -- CANNOT TELL IF I AM SIGNIFICANT. MAYBE OK?')
                
            ## ERROR in CI
            else:  
                logging.critical(" VARIANT: %s, ALLELE/GENOTYPE: %s (scoreOddRatio) -- CI ERROR (CI: %s)" % (self.variant_id,OR_data['OR_id'],OR_data['CI']))
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],'(scoreOddRatio) -- VARIANT: %s, ALLELE/GENOTYPE: %s (scoreOddRatio) -- CI ERROR (CI: %s)' % (self.variant_id,OR_data['OR_id'],OR_data['CI']))
#                 sys.exit(" VARIANT: %s, ALLELE/GENOTYPE: %s (scoreOddRatio) -- CI ERROR (CI: %s)" % (self.variant_id,OR_data['OR_id'],OR_data['CI']))
                OR_data['is_OK'] = False
        return OR_data
    ### END OF scoreOddRatio
    
    ### getSignificanceFromOrAndCi ###
    #  define if an OR is risk/protective allele of the model 
    # from Confidence interval only
    # called by scoreOddRatio
    ### 
    def getSignificanceFromOrAndCi(self,OR):
        logging.debug(' Function: getSignificanceFromOrAndCi:  OR %s' % (OR['OR_id']))
    
        ## OR protective
        if  OR['CI'][1] < 1 and OR['odd_ratio'] < 1:
            logging.info(" allele/genotype: %s (OR & CI) -- I am protective"% (OR['OR_id']))
            OR['is_prot'] = True 
            OR['is_sign'] = True  
            OR['is_OK'] = True
            OR['odd_ratio'] = round(1.0/OR['odd_ratio'], 5)
        
        ## OR risk
        elif OR['CI'][0] > 1 and OR['odd_ratio'] > 1:
            logging.info(" allele/genotype: %s (OR & CI) -- I am risky" % (OR['OR_id']))
            OR['is_risk'] = True 
            OR['is_sign'] = True  
            OR['is_OK'] = True
            
        ## non significant but OK?
        else :
            logging.info(" allele/genotype: %s (OR & CI) -- I am not significant"% (OR['OR_id']))
            OR['is_OK'] = True     
            
        ## check that pval is also significant = data validation
        OR = self.checkPvalueCiAgree(OR)    
        
        return OR
    ### END OF getSignificanceFromCiOnly

    ### checkPvalueCiAgree ###
    # called by 
    # - getSignificanceFromCiOnly
    # - getSignificanceFromOrAndCi
    ###  
    def checkPvalueCiAgree(self, OR):
        logging.debug(' Function: checkPvalueCiAgree:  OR %s' % (OR['OR_id']))
        
        ### pvalue not significant, CI significant
        if OR['pvalue'] != -1 and OR['pvalue'] >= 0.05 and OR['is_sign'] == True :
            OR['is_OK'] = False
            logging.warning(" allele/genotype: %s (CI & pvalue) -- CI & PVALUE DO NOT MATCH (+/-)" % (OR['OR_id']))  
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(CI & pvalue) -- CI & PVALUE DO NOT MATCH (+/-)')
        
        ### pvalue significant, CI not significant
        elif OR['pvalue'] != -1 and OR['pvalue'] < 0.05 and OR['is_sign'] == False :
            OR['is_OK'] = False
            logging.warning(" allele/genotype: %s (CI & pvalue) -- CI & PVALUE DO NOT MATCH (-/+)" % (OR['OR_id']) )
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(CI & pvalue) -- CI & PVALUE DO NOT MATCH (-/+)')

        ### pvalue & CI agree significant
        elif OR['pvalue'] != -1 and OR['pvalue'] < 0.05 and OR['is_sign'] == True:
            logging.info(" allele/genotype: %s (CI & pvalue) -- pvalue is significant too" % OR['OR_id'])
        
        ### pvalue & CI agree no significant
        elif OR['pvalue'] != -1 and OR['pvalue'] >= 0.05 and OR['is_sign'] == False:
            logging.info(" allele/genotype: %s (CI & pvalue)  -- pvalue is not significant too" % OR['OR_id'])
         
        ### pvalue not provided   
        elif OR['pvalue'] == -1 : 
            logging.warning(" allele/genotype: %s (CI & pvalue) -- NO PVALUE PROVIDED" % OR['OR_id'])
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(CI & pvalue) -- NO PVALUE PROVIDED')
        
        ### SHOULD NEVER HAPPEN?
        else:
            criticaltext = " VARIANT: %s, ALLELE/GENOTYPE: %s (CI & pvalue) -- SHOULD NOT BE HERE????????? (OR: %s)" % (self.variant_id,OR['OR_id'], json.dumps(OR, sort_keys = True, indent = 4) )
            logging.critical(criticaltext)
            sys.exit(criticaltext)
        return OR
    ### END OF checkPvalueCiAgree
    
    ### getSignificanceFromCiOnly ###
    #  define if an OR is risk/protective allele of the model 
    # from Confidence interval only
    # called by scoreOddRatio
    ### 
    def getSignificanceFromCiOnly(self,OR):
        logging.debug(' Function: getSignificanceFromCiOnly: OR %s' % (OR['OR_id']))
        ## CI protective
        if  OR['CI'][1] < 1 :
            logging.info(" allele/genotype: %s (CI only) -- CI is protective"% (OR['OR_id']))
            OR['is_prot'] = True 
            OR['is_sign'] = True 
            OR['is_OK'] = True
        
        ## OR risk
        elif OR['CI'][0] > 1 :
            logging.info(" allele/genotype: %s (CI only) -- CI is risky"% (OR['OR_id']))
            OR['is_risk'] = True 
            OR['is_sign'] = True 
            OR['is_OK'] = True
            
        ## non significant but OK?
        else :
            logging.info(" allele/genotype: %s (CI only) -- CI is not significant"% (OR['OR_id']))
            OR['is_OK'] = True     
            
        ## check that pval is also significant = data validation
        OR = self.checkPvalueCiAgree(OR)
            
        return OR
    ### END OF getSignificanceFromCiOnly
   
    ### getSignificanceFromPvalueOnly ###
    #  define if an OR is risk/protective allele of the model 
    # from pvalue only
    # called by scoreOddRatio
    ### 
    def getSignificanceFromPvalueOnly(self,OR):  
        logging.debug(' Function: getSignificanceFromPvalueOnly: OR %s' % ( OR['OR_id']))              
        ## pval significant?
        if OR['pvalue'] < config.PVALUE :
            OR['is_sign'] = True
            
            ## protective
            if OR['odd_ratio'] < 1 and OR['odd_ratio'] > 0 : 
                logging.info(" allele/genotype: %s (pvalue only) -- pvalue is protective"% ( OR['OR_id']))
                OR['is_prot'] = True  
                OR['is_OK'] = True
                OR['odd_ratio'] = round(1.0/OR['odd_ratio'],5)
                
            ## risk allele    
            elif  OR['odd_ratio'] > 1: 
                logging.info(" Iallele/genotype: %s (pvalue only) -- pvalue is risky"% ( OR['OR_id']))
                OR['is_risk'] = True 
                OR['is_OK'] = True
            
            ## significant OR= -1 ok
            elif OR['odd_ratio'] <= 0:
                OR['is_OK'] = True
                OR['is_sign'] = True
                logging.warning(" allele/genotype: %s (pvalue only) -- PVALUE SIGNIFICANT BUT I DO NOT HAVE OR" % ( OR['OR_id']))
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(pvalue only) -- PVALUE SIGNIFICANT BUT I DO NOT HAVE OR')
            ## Not significant OR=1 is OK
            else :
                OR['is_OK'] = False
                OR['is_sign'] = False
                logging.warning(" allele/genotype: %s (pvalue only) -- PVALUE SIGNIFICANT BUT I AM A REFERENCE" % ( OR['OR_id']))
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(pvalue only) -- PVALUE SIGNIFICANT BUT I AM A REFERENCE')
        
        ## Not significant is OK
        else :
            logging.info("  allele/genotype: %s (pvalue only) -- pvalue is not significant")
            OR['is_OK'] = True
            
        return OR
    ### END OF checkPvalueCiAgree    
    
    ### updateModelScore ###
    #  define if an OR is risk/protective allele of the model 
    # update model scores given worth values
    # called by model_scoring
    ###     
    def updateModelScore(self,model, OR):
        logging.debug(' Function: updateModelScore: model %s, OR %s' % (model['model_id'], OR['OR_id']))
        
        ### Define allele from allele or genotype
        test_allele_genotype = len(OR['OR_id']) == 1 or (len(OR['OR_id']) == 2 and OR['OR_id'][0] == OR['OR_id'][1])
        
        ### OR has no error
        if OR['is_OK'] == True :   
            ### OR for reference allele 
            if OR['is_ref'] == True:
                ## record ref allele if allele or homozygous genotype
                if test_allele_genotype:
                    logging.info(" model: %s, allele: %s (updateModelScore) -- I am the reference allele" % (model['model_id'], OR['OR_id'][0]))
                    model['reference_allele'] = OR['OR_id'][0]
                
                ## problem allele definition
                else :
                    logging.info(" model: %s, genotype: %s (updateModelScore) -- I am a reference heterozygous genotype" % (model['model_id'], OR['OR_id']))
                  
            ### OR sifnificant
            elif  OR['is_sign'] == True :
                model['significant'] = True
    
                ## allele or homozygous genotype risk 
                if test_allele_genotype == True and  OR['is_risk'] == True:
                        logging.info(" model: %s, allele: %s (updateModelScore) -- I am the risk allele" % (model['model_id'], OR['OR_id'][0]))
                        model['risk_allele'] = OR['OR_id'][0]
                
                ## allele or homozygous genotype  protective allele
                elif test_allele_genotype == True and OR['is_prot'] == True:
                    logging.info(" model: %s, allele: %s (updateModelScore) -- I am the protective allele" % (model['model_id'], OR['OR_id'][0]))  
                    model['protective_allele'] = OR['OR_id'][0]
                
                ## heterozygous
                else :
                    logging.info(" model: %s, genotype: %s (updateModelScore) -- I am an heterozygous genotype with OR data" % (model['model_id'], OR['OR_id']))
  
                ## check type is number
                OR['CI'][1]+1
                OR['odd_ratio'] +1
                OR['pvalue'] +1
                
                ### define worth (max) CI error for model comparison  
                error = round(OR['CI'][1] - OR['CI'][0],5) # CI different
                if  OR['CI'][1] != -1 and (model['max_error'] == -1 or error > model['max_error']): #max error if available
                    logging.info(" model: %s, allele/genotype: %s (updateModelScore) -- I have max_error" % (model['model_id'], OR['OR_id']))           
                    model['max_error'] = error
                else:
                    logging.info(" model: %s, allele/genotype: %s (updateModelScore) -- I do not have max_error" % (model['model_id'], OR['OR_id'])) 
                
                ### define worth (max) pvalue for model comparison 
                if  OR['pvalue'] != -1  and (model['max_pvalue'] == -1 or OR['pvalue'] > model['max_pvalue']): # max pvalue if available
                    logging.info(" model: %s, allele/genotype: %s (updateModelScore) -- I have max_pvalue" % (model['model_id'], OR['OR_id']))           
                    model['max_pvalue'] = OR['pvalue']
                else:
                    logging.info(" model: %s, allele/genotype: %s (updateModelScore) -- I do not have max_pvalue " % (model['model_id'], OR['OR_id']))  
                
                ### define worth (min) OR for model comparison 
                if OR['odd_ratio'] != -1 and OR['odd_ratio'] > 1.0 and (model['min_OR'] == -1 or OR['odd_ratio']  < model['min_OR']):
                    logging.info(" model: %s, allele/genotype: %s (updateModelScore) -- I have min_OR " % (model['model_id'], OR['OR_id']))
                    model['min_OR'] = OR['odd_ratio'] 
                    model['OR_OK'] = True
                    STAT['variant_scoring'][self.ethnicity]['total_OK_OR_model'] += 1
                
                ### No OR data
                else:
                    logging.info(" model: %s, allele/genotype: %s (updateModelScore) -- I do not have min_OR" % (model['model_id'], OR['OR_id']))
                    STAT['variant_scoring'][self.ethnicity]['total_NOTOK_OR_model'] += 1
            
            ### not significant and not reference
            else :  
                logging.warning(" model: %s, allele/genotype: %s (updateModelScore) -- I AM NOT SIGNIFICANT"  % (model['model_id'], OR['OR_id']))
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(updateModelScore) -- I AM NOT SIGNIFICANT')

        ### OR data error
        else :  
            warningtext = '(updateModelScore) -- OR DATA ERROR'
            logging.error(" model: %s, allele/genotype: %s %s (OR: %s)"  % (model['model_id'], OR['OR_id'], warningtext, json.dumps(OR, sort_keys = True, indent = 4)))    
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],warningtext)
                    
        return model       
    ### END OF updateModelScore
    
    ### updateModelType ###
    # - finalizes if reference if risk or protective allele for model
    # called by model_scoring
    ###   
    def updateModelType(self,model,nOR):
        logging.debug(' Function: updateModelType: model: %s, number of OR: %s' % (model['model_id'], nOR))

        ### allelic OR
        if nOR == 2:
            logging.info(' model: %s (updateModelType) -- I am allelic' % (model['model_id']))
            model['model_type'] = 0
        
        ### genotypic OR
        elif nOR == 3:
            logging.info(' model: %s (updateModelType) -- I am genotypic' % (model['model_id']))
            model['model_type'] = 1
            
        ### haplotypic OR
        elif nOR > 3:
            logging.info(' model: %s (updateModelType) -- I am haplotypic' % (model['model_id']))
            model['model_type'] = 2
        else:
            criticaltext = " VARIANT: %s, MODEL: %s (updateModelType) --I DO NOT KNOW WHAT TYPE OF OR I AM (MODEL: %s)" % (self.variant_id,model['model_id'],json.dumps(model, sort_keys=True, indent=4))
            logging.critical(criticaltext)
            sys.exit(criticaltext)
        
        return model
    ### END OF updateModelSignificance
    
    ### updateModelSignificance ###
    # - finalizes if reference if risk or protective allele for model
    # called by model_scoring
    ###   
    def updateModelSignificance(self,model):
        logging.debug(' Function: updateModelSignificance: model %s' % (model['model_id']))

        ### significant association-> proceed
        if model['significant'] == True:     
            logging.info(" model: %s (updateModelSignificance) -- I am significant" % model['model_id'])  
            STAT['variant_scoring'][self.ethnicity]['total_sign_model']  += 1  
            model = self.defineAllelesSignificance(model)
                    
        ### not significant 
        else:
            logging.info(" model: %s (updateModelSignificance) -- not significant" % model['model_id'])  
            STAT['variant_scoring'][self.ethnicity]['total_nsign_model']  += 1
            model = self.defineAllelesSignificance(model)
        return model
    ### END OF updateModelSignificance
    
    ### defineAllelesSignificance ###
    # - define if reference is risk or protective
    # called by updateModelSignificance
    ###  
    def defineAllelesSignificance(self,model):
        logging.debug(' Function: defineAllelesSignificance: model %s' % (model['model_id']))
        ## ref=risk allele
        if model['risk_allele'] == "" and model['reference_allele'] != "" and model['protective_allele'] !="":
            logging.info(" model: %s, allele: %s (defineAllelesSignificance) -- I am the reference and risk allele" % (model['model_id'],model['reference_allele']))               
            model['risk_allele'] = model['reference_allele']
            model['is_OK'] = True
            STAT['variant_scoring'][self.ethnicity]['total_OK_model'] += 1
    
        ## ref=protective allele
        elif model['protective_allele'] == "" and model['reference_allele'] != "" and model['risk_allele'] !="":
            logging.info(" model: %s, allele: %s (defineAllelesSignificance) -- I am the reference and protective" % (model['model_id'],model['reference_allele']))               
            model['protective_allele'] = model['reference_allele']
            model['is_OK'] = True
            STAT['variant_scoring'][self.ethnicity]['total_OK_model'] += 1
             
        ## error: couldn't find reference allele
        elif model['reference_allele'] == "" :
            errortext = '(defineAllelesSignificance) -- COULD NOT FIND REFERENCE ALLELE'
            logging.error(" VARIANT: %s, MODEL: %s %s (MODEL: %s)" % (self.variant_id,model['model_id'],errortext,json.dumps(model,sort_keys=True, indent=4)))
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],'%s' % errortext)
            STAT['variant_scoring'][self.ethnicity]['total_NOTOK_model'] += 1
        
        ## model is ok is non significant
        else:
            model['is_OK'] = True
            for OR in model['OR_data']:
                if model['OR_data'][OR]['is_OK'] == False:
                    model['is_OK'] = False
                    break
            ## error: SHOULD NEVER BE HERE
            if  model['is_OK'] == False and model['OR_OK'] == False:
                errortext = '(defineAllelesSignificance) -- COULD NOT FIND ANY ALLELE'
                logging.error(" VARIANT: %s, MODEL: %s %s (MODEL: %s)" % (self.variant_id,model['model_id'],errortext,json.dumps(model,sort_keys=True, indent=4)))
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],errortext)
                STAT['variant_scoring'][self.ethnicity]['total_NOTOK_model'] += 1
                #sys.exit(" VARIANT: %s, MODEL: %s %s (MODEL: %s)" % (self.variant_id,model['model_id'],errortext,json.dumps(model,sort_keys=True, indent=4)))
            else :
                if model['significant'] == True:
                    errortext = '(defineAllelesSignificance) -- ADD TO CHANGE SIFGNIFICANCE'
                    logging.error(" VARIANT: %s, MODEL: %s %s (MODEL: %s)" % (self.variant_id,model['model_id'],errortext,json.dumps(model,sort_keys=True, indent=4)))
                    SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],errortext)
#                     model['significant'] = False
                logging.info(" model: %s (defineAllelesSignificance) -- I am not significant but OK" % (model['model_id']))  
                STAT['variant_scoring'][self.ethnicity]['total_OK_model'] += 1
                
        return model
    ### END OF defineAllelesSignificance

    
    ### updateEthnicityScore ###
    # update the means and significance for that ethnicity given this sample
    # called by scoreVariants
    ### 
    def updateEthnicityScore(self, ethnicity, sample):
        logging.debug(' Function: updateEthnicityScore: ethnicity %s, sample %s' % (ethnicity['ethnicity_id'], sample['sample_id']))
                
        ## update meta-analysis (if sample is a subsample of a meta-analysis)
#         ethnicity = self.updateMetaAnalysisScore(ethnicity, sample) # update meta analysis score
        
        ### check number types
        sample['study_type'] +1
        sample['sample_size'] +1 
        
        ### no best model error -> proceed to calculate means and other statistics
        if sample['best_model']['is_OK'] == True :
            #if sample['best_model']['min_OR'] != -1 or sample['best_model']['max_error'] != -1 or sample['best_model']['max_pvalue'] != -1:
            if sample['study_type'] != -1 :
                ethnicity['mean_study_type'] += sample['study_type']
            if sample['sample_size'] != -1 :
                ethnicity['mean_sample_size'] += sample['sample_size']
            if sample['best_model']['min_OR'] != -1 :
                ethnicity['mean_OR'] += sample['best_model']['min_OR']
            if  sample['best_model']['max_error'] != -1:
                ethnicity['mean_error'] += sample['best_model']['max_error']
            if  sample['best_model']['max_pvalue'] != -1:
                ethnicity['mean_pvalue'] += sample['best_model']['max_pvalue']
            
            ## count significant samples for this ethnicity
            if sample['best_model']['significant'] == True:
                ethnicity['ethnicity_sign_sample'] += 1
                ethnicity['ethnicity_num_sample'] += 1
                STAT['variant_scoring'][self.ethnicity]['total_sign_sample'] += 1
                STAT['variant_scoring'][self.ethnicity]['total_num_sample'] += 1
            else: 
                ethnicity['ethnicity_nsign_sample'] += 1
                ethnicity['ethnicity_num_sample'] += 1
                STAT['variant_scoring'][self.ethnicity]['total_nsign_sample'] += 1
                STAT['variant_scoring'][self.ethnicity]['total_num_sample'] += 1
            
            ## loop on risk alleles found in this sample -> update risk allele count in this ethnicity
            for risk_allele in sample['sample_risk_allele']:
                if risk_allele in  ethnicity['ethnicity_risk_allele']:
                    ethnicity['ethnicity_risk_allele'][risk_allele] += 1 ##! add up only 1 count per sample otherwise over count with meta analysis
                else  :
                    ethnicity['ethnicity_risk_allele'][risk_allele] = 1
                    
        return ethnicity
    ### END OF updateEthnicityScore
    
    ### updateMetaAnalysisScore ###
    # update score of meta analysis
    ### 
    def updateMetaAnalysisScore(self, ethnicity, sample):        
        logging.debug(' Function: updateMetaAnalysisScore: sample %s' % (sample['sample_id']))
        smeta = SharedFunctions.convertToNone(sample['meta_sample'])
        ### no error and sample part of a meta analysis
        if sample['best_model']['is_OK'] == True and smeta is not None:        
            logging.info(" sample: %s (updateMetaAnalysisScore) -- part of meta samples" % (sample['sample_id']))
            
            ### loop on meta analysis samples
            for meta_sa in sample['meta_sample']:
                sample_name  = '%s_%s'% (sample['sample_id'].split('_')[0], meta_sa)
                if sample_name in ethnicity['sample_score']:
                    meta= ethnicity['sample_score'][sample_name]
                    ### Assess significance of sample
                    if sample['best_model']['significant'] == True:
                        meta['sign_sample'] += 1
                    else:
                        meta['nsign_sample'] += 1
                    ### update meta sample risk allele
                    risk_allele = sample['best_model']['risk_allele'] 
                    if risk_allele != "":   
                        ## known risk allele
                        if risk_allele in  meta['sample_risk_allele']:
                            logging.info(" sample: %s, risk allele: %s (updateMetaAnalysisScore) -- updating sample_risk_allele count"  % (sample['sample_id'], risk_allele))
                            meta['sample_risk_allele'][risk_allele] += 1
                        ## new risk allele
                        else  :
                            logging.info(" sample: %s, risk allele: %s (updateMetaAnalysisScore) -- creating sample_risk_allele"  % (sample['sample_id'], risk_allele))
                            meta['sample_risk_allele'][risk_allele] = 1   
                else:
                    logging.warning(" NO META SAMPLE sample: %s (updateMetaAnalysisScore)  - meta name: %s"  % (sample['sample_id'], sample_name))    
        else :
            logging.info(" sample: %s (updateMetaAnalysisScore) -- Not part of meta sample"  % (sample['sample_id']))
        return ethnicity   
    ### END OF updateMetaAnalysisScore
        
    ### checkRiskAlleleForBestSample ###
    # find max risk allele for ethnicity
    # check that best sample_risk_allele fits
    # if not find best sample again
    # called by scoreVariants
    ###
    def checkRiskAlleleForBestSample(self,ethnicity):
        logging.debug(' Function: checkRiskAlleleForBestSample: ethnicity %s' % (ethnicity['ethnicity_id']))
        
        ### Find best Risk alleles of ethnicity and its best sample
        ethnicity_risk_allele = self.defineRiskAllele(ethnicity['ethnicity_risk_allele'] )
        sample_risk_allele = self.defineRiskAllele(ethnicity['best_sample']['sample_risk_allele'] )
    
        ### undetermines risk allele
        if ethnicity_risk_allele == None or  sample_risk_allele == None: 
            logging.warning(" ethnicity: %s (risk allele: %s), best_sample: %s (risk allele: %s) (checkRiskAlleleForBestSample) -- CANNOT KEEP SAMPLE: BAD RISK ALLELE" % (ethnicity['ethnicity_id'],ethnicity_risk_allele,ethnicity['best_sample']['sample_id'], sample_risk_allele))
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(checkRiskAlleleForBestSample) -- CANNOT KEEP SAMPLE: BAD RISK ALLELE')
            ethnicity['best_sample']['best_model']['risk_allele_OK'] = False
         
        ### coud determine risk allele
        else: 
            ### test that best sample risk allele = best risk allele for ethnicity
            ok_risk_allele = self.checkRiskAlleleSampleFitEthnicity(sample_risk_allele, ethnicity_risk_allele)
            
            ### if best risk alleles don't match - find next best sample
            if ok_risk_allele != True and len(ethnicity['sample_score']) > 1:
                logging.warning(" ethnicity: %s (risk allele: %s), best_sample: %s (risk allele: %s) (checkRiskAlleleForBestSample) -- UNMATCHNING RISK ALLELES" % (ethnicity['ethnicity_id'],ethnicity_risk_allele,ethnicity['best_sample']['sample_id'], sample_risk_allele))
                SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(checkRiskAlleleForBestSample) -- UNMATCHNING RISK ALLELES')
                
                ## delete current best sample
                del ethnicity['sample_score'][str(ethnicity['best_sample']['sample_id'])]
                
                ## reinit best sample
                ethnicity['best_sample'] = (ethnicity['sample_score'][ethnicity['sample_score'].keys()[0]]) 
                
                ## Find new best sample
                for sa in ethnicity['sample_score']: 
                    logging.info("ethnicity: %s, sample: %s (checkRiskAlleleForBestSample) -- call compareSamples again" % (ethnicity['ethnicity_id'],sa ))
                    sample_score = ethnicity['sample_score'][sa]  
                    ethnicity['best_sample'] = SharedFunctions.compareSamples(ethnicity['best_sample'], sample_score,ethnicity['hapmap_freq'],True)
                
                ## rechek that best sample has good risk allele
                ethnicity = self.checkRiskAlleleForBestSample(ethnicity)
            else :
                logging.info(" ethnicity: %s (risk allele: %s), best_sample: %s (risk allele: %s) (checkRiskAlleleForBestSample) -- risk allele OK" % (ethnicity['ethnicity_id'],ethnicity_risk_allele,ethnicity['best_sample']['sample_id'], sample_risk_allele))
                ethnicity['best_sample']['best_model']['risk_allele_OK'] = True
                     
        return ethnicity
    ### END OF checkRiskAlleleForBestSample
    
    ### defineRiskAllele ###
    # find max risk allele for that array of risk alleles
    # called by checkRiskAlleleForBestSample
    ###
    def defineRiskAllele(self,RA_list):
        logging.debug(' Function: defineRiskAllele: risk_allele list %s' % (RA_list))
      
        ### initiate 
        max_count = -1.1
        n_allele = 0.0
        max_RA = []
        
        ### loop on risk alleles found so far
        for risk_allele in RA_list:
            n_allele += 1
            count = RA_list[risk_allele] 
            if count > max_count:
                max_RA = []
                max_count = count
                max_RA.append(risk_allele)
            elif count == max_count:
                max_RA.append(risk_allele) 
        
        ### check if can define risk allele
        if ( n_allele > 0 and  round( max_count/n_allele , 5) <= 0.5) or (n_allele == 0):
            logging.warning(" RISK ALLELE: %s (defineRiskAllele) -- CANNOT DEFINE RISK ALLELE" % (RA_list))
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['WARNING'],'(defineRiskAllele) -- CANNOT DEFINE RISK ALLELE')
            return None
        else: 
            logging.info(" risk allele: %s (defineRiskAllele) -- found risk allele" % (max_RA))
            return max_RA
    ### END OF defineRiskAllele
    
    ### checkRiskAlleleSampleFitEthnicity ###
    # returns whether a risk allele 1 is present in RA2
    # called by checkRiskAlleleForBestSample
    ###
    def checkRiskAlleleSampleFitEthnicity(self,RA1,RA2):
        logging.debug(' Function: checkRiskAlleleSampleFitEthnicity: risk_allele1 %s, risk_allele2 %s' % (RA1,RA2))
        
        for allele1 in RA1:
            for allele2 in RA2:
                if allele1 == allele2:
                    return True 
        return False
    ### END OF checkRiskAlleleSampleFitEthnicity
    
    ### updateVariantScore ###
    # update the sign ethnicity and risk allele for that variant given this sample
    # called by scoreVariants
    ### 
    def updateVariantScore(self,variant, ethnicity):
        logging.debug(' Function: updateVariantScore: variant %s, ethnicity %s' % (variant['variant_id'], ethnicity['ethnicity_id']))
          
        ### calculate mean (divide by nsign+sign)
        ethnicity = self.calculate_means(ethnicity)
        
        ### no best model error -> proceed to calculate means and other statistics
        sample = ethnicity['best_sample']
        if sample['best_model']['is_OK'] == True:
            ## count significant samples for this ethnicity
            if sample['best_model']['significant'] == True:
                variant['variant_sign_ethnicity'] += 1
                variant['variant_num_ethnicity'] += 1
                STAT['variant_scoring']['total_sign_ethnicity'] += 1
                STAT['variant_scoring']['total_num_ethnicity'] += 1
            else: 
                variant['variant_nsign_ethnicity'] += 1
                variant['variant_num_ethnicity'] += 1
                STAT['variant_scoring']['total_nsign_ethnicity'] += 1  
                STAT['variant_scoring']['total_num_ethnicity'] += 1
          
        ### loop on risk alleles found in this sample -> update risk allele count in this ethnicity
        for risk_allele in ethnicity['ethnicity_risk_allele']:
            STAT['variant_scoring'][self.ethnicity]['total_risk_allele'] += 1
            if risk_allele in  variant['variant_ethnicity_risk_allele']:
                variant['variant_ethnicity_risk_allele'][risk_allele] += 1
                variant['variant_total_risk_allele'][risk_allele] += ethnicity['ethnicity_risk_allele'][risk_allele]
            else :
                variant['variant_ethnicity_risk_allele'][risk_allele] = 1
                variant['variant_total_risk_allele'][risk_allele] = ethnicity['ethnicity_risk_allele'][risk_allele]
           
        ### prepare input for next module
        variant['ethnicity_score'][ethnicity['ethnicity_id']] = self.formatVariantScoreForSelection(ethnicity)
   
        return variant
    ### END OF updateVariantScore
    
    ### calculate_means ###
    # divide all means by sign+nsign
    # called by updateVariantScore
    ### 
    def calculate_means(self,ethnicity):
        logging.debug(' Function: calculate_means: ethnicity %s' % (ethnicity['ethnicity_id']))

        mean_list = {}
        total = 0
        ### find total and mean keys
        for key in ethnicity:
            if len(key.split('_')) > 1 :
                if key.split('_')[1] == 'sign' or key.split('_')[1] == 'nsign':
                    total += ethnicity[key]
                if key.split('_')[0] == 'mean':
                    mean_list[key]={} 
                    
        ### calculate means
        if total > 0:
            for key in mean_list:
                ethnicity[key] = (ethnicity[key] / total)
        else:
            errortext = '(calculate_means) -- COULD NOT CALCULATE MEAN SCORES'
            logging.error('VARIANT: %s, ETHNICITY: %s %s (ETHNICITY: %s)' % (self.variant_id,ethnicity['ethnicity_id'],errortext, json.dumps(ethnicity, sort_keys = True, indent =4)))
            SharedFunctions.fillWarningSummary(STAT['variant_scoring']['ERROR'],errortext)              
            #sys.exit('VARIANT: %s, ETHNICITY: %s %s (ETHNICITY: %s)' % (self.variant_id,ethnicity['ethnicity_id'], errortext, json.dumps(ethnicity, sort_keys = True, indent =4)))
            
        return ethnicity
    ### END OF calculate_means
    
    ### formatVariantScoreForSelection ###
    # prepare variant score for output and variant selection
    # clean/remove now useless key/values
    # called by updateVariantScore
    ### 
    def formatVariantScoreForSelection(self,eth):
        logging.debug(' Function: formatVariantScoreForSelection: ethniticy  %s' % (eth))
        
        variant_et = copy.deepcopy(eth)
        SF.WORKING[self.variant_id+eth['ethnicity_id']] = copy.deepcopy(eth) ### copy for test output
 
        del variant_et['hapmap_freq']
        del variant_et['sample_score'] 
        del variant_et['sample'] 
        del variant_et['best_sample']['association']
        del variant_et['best_sample']['best_model']['OR_data']
        
        ### clean best sample
        if 'dif' in variant_et['best_sample']:
            del variant_et['best_sample']['dif']
        if 'test' in variant_et['best_sample']:
            del variant_et['best_sample']['test']
        
        ### clean best model
        if 'dif' in variant_et['best_sample']['best_model']:
            del variant_et['best_sample']['best_model']['dif']
        if 'test' in variant_et['best_sample']['best_model']:
            del variant_et['best_sample']['best_model']['test']

        return variant_et