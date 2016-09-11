#!/usr/bin/python

'''
Created on Jul 30, 2012
Logs added 8/8/12
@author: celine
'''

import copy
import sys
import logging
import json

### Working dictionary
global WORKING
WORKING = {}
 
### Statistics
global STATISTICS 
STATISTICS = {}
STATISTICS['shared_functions'] = {}
STATISTICS['shared_functions']['WARNING'] = {}
STATISTICS['shared_functions']['ERROR'] = {}

class VariantSelectorSharedFunctions(object):
    
    def __init__(self):
        ''' '''

    ### checkIdFit ###
    # Check that the _id value fit the key or parent
    # called by:
    # - variant_scoring
    # - sample_scoring
    # - LD_variant_selection
    ### 
    def checkIdFit(self,id_key,id_value):
        textbase = 'key: %s, id_value: %s' % (id_key,id_value)
        logging.debug(' Function: checkIdFit: %s' % textbase )
        
        if id_key == id_value or id_key == str(id_value):
            logging.info(' %s (checkIdFit) -- key and id value fit' % (textbase))
        else:
            warningtext = '(checkIdFit) -- KEY AND ID VALUE DO NOT FIT'
            logging.critical("%s %s " %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
            sys.exit("%s %s " %(textbase,warningtext))
    ### END OF checkIdFit

    ### compareModels ###
    # - find best model  so far
    # called by:
    # - sample_scoring
    # - test_best_model
    ### 
    def compareModels(self,best, model, hapfreq,same_variant):
        textbase = 'best: %s, model: %s' % (best['model_id'], model['model_id'])
        logging.debug(' Function: compareModels: %s, same variant: %s' % (textbase,same_variant))
        
        if (same_variant == True and best['model_id'] != model['model_id'] or same_variant == False): 
            ####  check if have frequencies 
            hfreq = self.convertToBoolean(hapfreq)
            bfreq = self.convertToBoolean(best['control_freq'])
            mfreq = self.convertToBoolean(model['control_freq'])
            
            ### proceed if model and best are sign and no error
            if (hfreq == True or (hfreq == False and mfreq == True and bfreq == True)) and\
                    model['significant'] == True and \
                    model['is_OK'] == True and \
                    model['OR_OK'] == True and \
                    best['significant'] == True and \
                    model['is_OK'] == True and \
                    best['OR_OK'] == True:  
                return self.perform_model_comparison(best, model,same_variant)
                    
            ### best model not chosen
            elif  model['significant'] == True and \
                   model['is_OK'] == True and \
                   model['OR_OK'] == True and \
                   (best['significant'] == False or\
                    best['OR_OK'] == False or\
                    best['is_OK'] == False or\
                    (hfreq == False and mfreq == True and bfreq == False)):
                self.checkWhyBestModel(best,model, same_variant,hfreq ,mfreq,bfreq)
                return model 
            
            ### comparison not needed
            else :
                self.checkWhyModelComparisonImpossible(best,model, same_variant,hfreq ,mfreq,bfreq)
                return best
        else:
            logging.info(" %s (compareModels - same_variant %s) -- comparison not needed: best = model" % (textbase,same_variant))
            return best
    ### END OF compareModels
    
    ### checkWhyModelComparisonImpossible ###
    # warning on why can't compare models
    # Called by compareModels
    ###
    def checkWhyModelComparisonImpossible(self,b,m,sv,hfreq ,mfreq,bfreq ):
        textbase = ' best: %s, model: %s' % (b['model_id'], m['model_id'])
        logging.debug(' Function: checkWhyModelComparisonImpossible: %s, same variant: %s' % (textbase, sv))
        
        if m['significant'] == False:
            warningtext = '(compareModels - same_variant: %s) -- COMPARISON NOT NEEDED: MODEL IS NOT SIGNIFICANT' % sv
            logging.warning("%s %s " %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        elif m['OR_OK'] == False:
            warningtext = '(compareModels - same_variant: %s) -- COMPARISON NOT NEEDED: MODEL HAS NO OR DATA' % sv
            logging.warning("%s %s " %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        elif m['OR_OK'] == False:
            warningtext = '(compareModels - same_variant: %s) -- COMPARISON NOT NEEDED: MODEL IS NOT OK' % sv
            logging.warning("%s %" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        elif hfreq == False and mfreq == False:
            warningtext = '(compareModels - same_variant: %s) -- COMPARISON NOT NEEDED: MODEL HAS NO FREQUENCIES' % sv
            logging.warning("%s %s" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        else:
            criticaltext = " %s (compareModels - same_variant: %s) -- CANNOT PERFORM COMPARISON: WHY NOT (BEST/MODEL: {\"%s\":%s,\"%s\":%s})" % (textbase,sv, b['model_id'],json.dumps(b, sort_keys= True, indent = 4),m['model_id'],json.dumps(m, sort_keys= True, indent = 4) )
            logging.critical( "%s " %(criticaltext)) 
            sys.exit(" %s" %(criticaltext))              
    ### END OF checkWhyBestModel 
    
    ### checkWhyBestModel ###
    # warning on why new model is best
    # Called by compareModels
    ###
    def checkWhyBestModel(self,b,m,sv,hfreq ,mfreq,bfreq ):
        textbase = ' best: %s, model: %s' % (b['model_id'], m['model_id'])
        logging.debug(' Function: checkWhyBestModel: %s, same variant: %s' % (textbase,sv))
        
        if b['significant'] == False:
            warningtext = "(compareModels - same_variant: %s) -- BEST MODEL IS NOT SIGNIFICANT" % sv
            logging.warning("%s %s" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        elif b['OR_OK'] == False:
            warningtext = "(compareModels - same_variant: %s) -- BEST MODEL HAS NO OR DATA" % sv
            logging.warning("%s %s" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        elif b['OR_OK'] == False:
            warningtext = "(compareModels - same_variant: %s) -- BEST MODEL IS NOT OK" % sv
            logging.warning("%s %s" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        elif hfreq == False and mfreq == True and bfreq == False:
            warningtext = "(compareModels - same_variant: %s) -- BEST MODEL HAS NO FREQUENCIES" % sv
            logging.warning("%s %s" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        else:
            criticaltext = " %s (compareModels - same_variant: %s) -- MODEL BETTER THAN BEST: WHY? (BEST/MODEL: {\"%s\":%s,\"%s\":%s})" % (textbase,sv, b['model_id'],json.dumps(b, sort_keys= True, indent = 4),m['model_id'],json.dumps(m, sort_keys= True, indent = 4))
            logging.critical("%s" %(criticaltext))
            sys.exit("%s" %(criticaltext))
        
        ### check if model types fit
        if sv == True and b['model_type'] > m['model_type']:
            warningtext = '(compareModels - same_variant: %s) -- BEST MODEL TYPE < NEW MODEL TYPE' % sv
            logging.warning("%s %s (best type:%s, model_type: %s)" % (textbase,warningtext,b['model_type'], m['model_type']))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
    ### END OF checkWhyBestModel 
    
    ### perform_model_comparison ###
    # compare 2 different models 
    # Called by compareModels
    ###
    def perform_model_comparison(self,b,m,sv):
        textbase =' best: %s, model: %s' % (b['model_id'], m['model_id'])
        logging.debug(' Function: perform_model_comparison: %s, same variant: %s' % (textbase , sv))
        
        ### Initiate test and differences tables
        b['test'] = []
        b['dif'] = []
        m['test'] = []
        m['dif'] = []
        
        ### if genotype/allelic comparison and same variant -> take genotype data
        if (sv == True and b['model_type'] == m['model_type']) or sv == False :
            logging.info(" %s (compareModels - same variant: %s) -- perform comparison" % (textbase, sv))
            
            ## test which of best and model has best error, pvalue, OR
            b, m = self.test('max_error', b,m, "min", 1.5)
            b, m = self.test('max_pvalue', b,m, "min", 1)
            b, m = self.test('min_OR', b,m, "max", .5) 
            
            ## compare model type between 2 variants
            if b['model_type'] != m['model_type'] and sv == False:
                b, m = self.test('model_type', b,m, "max", 1) 
            
            ## find best model
            temp = self.findBest(b, m,'model_id')    
            ## variant_scoring keep best if undefined best model
            if temp == None and sv == True:
                warningtext = "(compareModels - same variant: %s) -- I KEEP BEST DESPITE INDECISION" % (sv)
                logging.warning("%s %s" %(textbase,warningtext))
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
                return b
            ## LD_selection do not choose best if cant define it
            if temp == None and sv == False:
                warningtext = "(compareModels - same variant: %s) -- INDECISION" % (sv)
                logging.warning("%s %s" %(textbase,warningtext))
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
                return None
            ##
            else:
                logging.info(" %s, temp: %s (compareModels - same variant: %s) -- clear best model" % (textbase,temp['model_id'], sv))
                return temp
        
        ### if genotype/allelic comparison and same variant -> take genotype data
        elif sv == True and b['model_type'] != m['model_type']:            
            warningtext = "(compareModels - same variant:%s) -- CHOOSING GENOTYPIC OR HAPLOTYPIC MODEL" % (sv)
            logging.warning("%s %s (best type: %s, model_type: %s)" %(textbase,warningtext,b['model_type'], m['model_type']))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
            
            ## find best model -> genotype or haplotype
            b, m = self.test('model_type', b,m, "max", 1) 
            temp = self.findBest(b, m,'model_id')    
            if temp == None and sv == True: # variant_scoring keep best
                return b
            if temp == None and sv == False:# LD selection NONE
                return None
            else: # return best
                return temp
            
        ### Can't figure why here SHOULD NEVER HAPPEN
        else:
            criticaltext = " %s (compareModels - same variant: %s) -- DONT KNOW HOW I GOT HERE (BEST/MODEL: {\"%s\":%s,\"%s\":%s})" % (textbase,b['model_id'], m['model_id'],b['model_id'],json.dumps(b, sort_keys= True, indent = 4),m['model_id'],json.dumps(m, sort_keys= True, indent = 4))
            logging.critical("%s %s" %(criticaltext))
            sys.exit("%s %s" %(criticaltext))
    ### END OF perform_model_comparison   
    
    ### convertToBoolean ###
    # convert a value into a boolean if in correct format
    # - true/false (json style boolean) or 
    # - "true"/"false" string (any cases are OK) or
    # - "yes"/"no" strings (any cases are OK) or
    # - 0/1
    # called by:
    # - compareModels
    # - compareSamples
    ###
    def convertToBoolean(self, value):
        logging.debug(' Function: convert_2_bolean: value: %s' % (value))
        
        ### true or yes
        if str(value).lower() == 'true' or str(value).lower() == 'yes' :
            return True
        ### false or no
        elif str(value).lower() == 'false' or str(value).lower() == 'no' :
            return False
        ### 0 or 1
        else :
            try:
                if  int(value) == 0:
                    return False  
                elif int(value) == 1:
                    return True 
                ## critical error if not 0 or 1 number
                else:
                    criticaltext = " VALUE: %s (convertToBoolean) -- COULD NOT CONVERT INTO BOOLEAN OR INTEGER"  % value
                    logging.critical(criticaltext)
                    sys.exit(criticaltext)
            
            ## critical error if not an interger or boolean like value
            except:
                criticaltext = " VALUE: %s (convertToBoolean) -- COULD NOT CONVERT INTO BOOLEAN"  % value
                logging.critical(criticaltext)
                sys.exit(criticaltext)
    ### END OF convert_2_boolean    
    
    ### correct_pvalue_test ###
    # weight down pvalue product and weight with max of pvalues
    # called by test
    ### 
    def testCorrectedPvalue(self,score, best, x, prod):
        logging.debug(' Function: correct_weight_for_pvalue -- score: %s, best: %s, x: %s, weight: %s' % (score, best, x, prod))
  
        if len(score.split('_')) > 1 :
            if score.split('_')[1] == 'pvalue':
                if best > x:
                    logging.info(' score: %s (correct_weight_for_pvalue) -- newweight: %s' % (score, round(prod * best, 1000) ))
                    return round(prod * best, 1000) 
                else:
                    logging.info(' score: %s (correct_weight_for_pvalue) -- newweight: %s' % (score, round(prod * x, 1000) ))
                    return round(prod * x, 1000) 
            else: 
                return prod            
        else:        
            return prod
    ### END OF testCorrectedPvalue 
    
    ### test ###
    # finds the min, or max (define by "comp") between 
    # the best value and a new value for a specific statictics "score".
    # weights define the importance of 'score' for future decisions
    # Calculate the difference between best and x values
    # record if best or x "win" the comparison test
    # called by:
    # - compareModels
    # - compareSamples
    # - compare_variants
    ### 
    def test(self,score, best, x, comp, weight):
        textbase = 'score: %s' % score
        logging.debug(' Function: test -- %s, best: %s, x: %s, weight: %s' % (textbase, best[score], x[score], weight))

        ## find statistic value to compare
        b = self.convertToNone(best[score])
        c = self.convertToNone(x[score])
        
        ### best has not data: x wins
        if  (b is None ) and c is not None :
            logging.info(" %s (test) -- x wins by default" % (textbase))
            best['test'].append(0)
            x['test'].append(weight)
            best['dif'].append(0)
            x['dif'].append(c)   
        
        ### x has not data best: wins 
        elif  (c is None ) and b is not None :
            logging.info(" %s (test) -- best wins by default" % (textbase))
            best['test'].append(weight)
            x['test'].append(0)
            best['dif'].append(b)
            x['dif'].append(0)   
        
        ### x and best have data -> compare
        elif b is not None and c is not None:
            
            dif = round((b - c) * weight, 1000)
            if c !=0 and b !=0:
                if dif>0:
                    prod = (round(b,1000)/round(c,1000)) * round(weight,1000)
                else :
                    prod = (round(c,1000)/round(b,1000)) * round(weight,1000)
            else:
                prod = abs(dif)
            logging.info(' %s %s' %(dif, prod))
            prod = self.testCorrectedPvalue(score,b,c,prod)
            ## b>c
            if dif > 0 : 
                ## x has min score
                if comp == 'min': #c
                    
                    logging.info(" %s (test) -- x is minimum (dif>0)" % (textbase))
                    best['test'].append(0)
                    x['test'].append(weight)
                    best['dif'].append(0)
                    x['dif'].append(prod)
                
                ## best has max score
                elif comp == 'max': #b
                    logging.info(" %s (test) -- best is maximum (dif>0)" % (textbase))
                    best['test'].append(weight)
                    x['test'].append(0)
                    best['dif'].append(prod)
                    x['dif'].append(0)
                    
                ## comp error SHOULD NEVER HAPPEN    
                else:
                    textbase = '%s, comp: %s' % (textbase, comp)
                    errortext = "(test) -- I DONT UNDERSTANT COMPARISON TAG (dif>0)" 
                    logging.error(" %s %s" %(textbase,errortext))
                    self.fillWarningSummary(STATISTICS['shared_functions']['ERROR'], errortext)
                
            ## b<c
            elif dif <0 :
                ## best has min score
                if comp == 'min': #b
                    logging.info(" %s (test) -- best is minimum (dif<0)" % (textbase))
                    best['test'].append(weight)
                    x['test'].append(0)
                    best['dif'].append(prod)
                    x['dif'].append(0)
                
                ## x has max score
                elif comp == 'max': #c
                    logging.info(" %s (test) -- x is maximum (dif<0)" % (textbase))
                    best['test'].append(0)
                    x['test'].append(weight)
                    best['dif'].append(0)
                    x['dif'].append(prod)
                    
                ## comp error SHOULD NEVER HAPPEN    
                else:
                    textbase = '%s, comp: %s' % (textbase, comp)
                    errortext = "(test) -- I DONT UNDERSTANT COMPARISON TAG (dif<0)" 
                    logging.error(" %s %s" %(textbase,errortext))
                    self.fillWarningSummary(STATISTICS['shared_functions']['ERROR'], errortext)
            ## no difference  
            else:
                logging.info(" %s (test) -- no difference" % (textbase))
                
        ## no value to compare  
        else:
            warningtext = "(test) -- CAN NOT COMPARE VALUES: NO DATA FOUND"
            logging.warning(" %s %s" %(textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
        return best, x
    ### END OF test
    
    ### findBest ###
    # find the best out of pre-tested objects
    # call itself if unconclusive between test and dif
    # Called by
    # - findBest
    # - compareModels
    # - compareSamples
    # - compare_variants
    ### 
    def findBest(self,b, c, name_id):
        logging.debug(' Function: findBest: best_test %s, x_test %s, best_dif %s, x_dif %s' % ((b['test']),(c['test']) ,  (b['dif']), (c['dif'])))
        
        ### find object and id compared
        b_key, b_name = self.findNameId(b, name_id)
        c_key, c_name = self.findNameId(c, name_id)
        textbase = " object: %s, best %s, new %s" % (b_key, b_name, c_name)
        if b_key == c_key:
            ### best test and dif < c -> c is best
            if sum(b['test']) < sum(c['test']) and  sum(b['dif']) < sum(c['dif']):
                logging.info(" %s (findBest) -- new is best" % (textbase))
                b = c
                
            ### best test and best >= c -> best remain best
            elif sum(b['test']) > sum(c['test']) and  sum(b['dif']) > sum(c['dif']):
                logging.info(" %s (findBest) -- best is still best" % (textbase))
            elif sum(b['test']) == sum(c['test']) and  sum(b['dif']) > sum(c['dif']):
                warningtext = "(findBest) --  best is still best from products"
                logging.warning(" %s %s (test: %s, %s, dif %s, %s))" %(textbase,warningtext,b['test'],c['test'] , b['dif'], c['dif']))
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
     
                #logging.info(" %s (findBest) -- best is still best from products" % (textbase))
            elif sum(b['test']) == sum(c['test']) and  sum(b['dif']) < sum(c['dif']):
                warningtext = "(findBest) -- new is best from products"
                logging.warning(" %s %s (test: %s, %s, dif %s, %s))" %(textbase,warningtext,b['test'],c['test'] , b['dif'], c['dif']))
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
                #logging.info(" %s (findBest) --  new is best from products" % (textbase))
                b = c
            elif  len(b['test']) == 0:
                warningtext = "(findBest) -- CANNOT DECIDE WHO IS BEST"
                logging.warning(" %s %s (test: %s, %s, dif %s, %s))" %(textbase,warningtext,b['test'],c['test'] , b['dif'], c['dif']))
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
                return None
           
            ### test and dif disagree -> remove last (less important) statistics and try again
            elif len(b['test']) > 1:                
                warningtext = "(findBest) -- NEED TO CALL FIND_BEST AGAIN"
                logging.warning(" %s %s (test: %s, %s, dif %s, %s))" %(textbase,warningtext,b['test'],c['test'] , b['dif'], c['dif']))
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'], warningtext)
                b['test'].pop()
                c['test'].pop()
                b['dif'].pop()
                c['dif'].pop()
                b = self.findBest(b,c,name_id)
                
            ### test and dif still disagree ERROR- SHOULD NEVER BE HERE
            else:
                c['is_OK'] = False
                criticaltext = " %s (findBest) -- COULD NOT FIND BEST" %textbase
                logging.critical("%s" %(criticaltext))
                sys.exit("%s" %(criticaltext))
        ### KEYS DONT MATCH- SHOULD NEVER BE HERE
        else:
            criticaltext = " OBJECT1: %s (VALUE: %s), OBJECT2: %s (VALUE: %s) (findBest) -- I AM COMPARING UNICORN WITH LEPRECHAUN!! " % (b_key,b_name, c_key, c_name)
            logging.critical("%s" %(criticaltext))
            sys.exit("%s" %(criticaltext))
        return b
    ### END OF findBest
         
    ### compareSamples ###
    # loop on properties to find key and id
    # Called by findBest
    ### 
    def findNameId(self,x, name_id):
        logging.debug(' Function: findNameId: x: %s, name_id: %s ' % (x,name_id))
        k = ""
        val = -1
        for key in x:
            if name_id in key:
                logging.info(' (findNameId) -- got an ID %s'% key)
                k = key
                val = x[key]
                break
        return k, val
    ### END OF findNameId    
                                       
    ### compareSamples ###
    # compare a sample to best so far
    # Called by:
    # - variant_scoring
    # - find_best_variant_from_sample_score
    ### 
    def compareSamples(self,best, sample, hapfreq,same_variant):
        logging.debug(' Function: compareSamples: best: %s, sample %s, same_variant: %s' % (best['sample_id'], sample['sample_id'],same_variant))
        
        ### compare only different samples
        if (best['sample_id'] != sample['sample_id'] and same_variant == True) or same_variant == False:  

            ####  check if have frequencies 
            hfreq = self.convertToBoolean(hapfreq)
            bfreq = self.convertToBoolean(best['best_model']['control_freq'])
            sfreq = self.convertToBoolean(sample['best_model']['control_freq'])
            
            ### Meta-analysis of each other is same variant?
            if same_variant == True:
                bmeta = self.isSampleMetaAnalysis(best['sample_id'], sample['meta_sample'])
                smeta = self.isSampleMetaAnalysis(sample['sample_id'], best['meta_sample'])
            else:
                bmeta = False
                smeta = False
            
            ### proceed if sample.best model and bestsample.bestmodel are sign and no error and freq
            if ((bmeta == False and smeta == False) or (bmeta == True and smeta == True)) and\
                (hfreq == True or (hfreq == False and sfreq == True and bfreq == True)) and\
                sample['best_model']['significant'] == True and\
                sample['best_model']['OR_OK'] == True and\
                sample['best_model']['is_OK'] == True and\
                best['best_model']['significant'] == True and\
                best['best_model']['OR_OK'] == True and \
                best['best_model']['is_OK'] == True:
                logging.info(" best: %s, sample: %s (compareSamples - same_variant: %s) -- perform comparison" % (best['sample_id'], sample['sample_id'],same_variant))
                
                ## compare sample_score
                best, sample = self.findBestScoreOf2Samples(best,sample)
                
                ## compare best models only if needed
                best, sample = self.findBestModelOf2Samples(best, sample,hfreq,same_variant)
                
                ## find best sample
                best = self.findBest(best, sample,'sample_id' )    
                
            ### best sample was not scored yet (should never happen) or no frequencies
            elif  sample['best_model']['significant'] == True and \
                   sample['best_model']['is_OK'] == True and \
                   sample['best_model']['OR_OK'] == True and\
                   bmeta == False and \
                   (best['best_model']['significant'] == False or \
                    best['best_model']['OR_OK'] == False or\
                    best['best_model']['is_OK'] == False or\
                    (hfreq == False and sfreq == True and bfreq == False) or\
                    smeta == True ):
                self.checkWhyBestSample(best,sample,hfreq ,sfreq,bfreq,bmeta , smeta ,same_variant )
                best = sample 
            else:
                self.checkWhySampleComparisonImpossible(best,sample,hfreq ,sfreq,bfreq,bmeta , smeta ,same_variant )
        return best
    ### END OF compareSamples
     
    ### checkWhySampleComparisonImpossible ###
    # warning on why can't compare samples
    # Called by compareSamples
    ###
    def checkWhySampleComparisonImpossible(self,b,s,hfreq ,sfreq,bfreq,bmeta,smeta,sv ):
        textbase = 'best :%s, sample: %s' % (b['sample_id'],s['sample_id'])
        logging.debug(' Function: checkWhySampleComparisonImpossible: %s' % textbase)
        
        if s['best_model']['significant'] == False:
            warningtext= " (checkWhySampleComparisonImpossible - sample_variant: %s) -- COMPARISON NOT NEEDED: SAMPLE BEST_MODEL IS NOT SIGNIFICANT" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
        
        elif s['best_model']['OR_OK'] == False:
            warningtext= " (checkWhySampleComparisonImpossible - sample_variant: %s) -- COMPARISON NOT NEEDED: SAMPLE BEST_MODEL HAS NO OR DATA" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
        
        elif s['best_model']['is_OK'] == False:
            warningtext= " (checkWhySampleComparisonImpossible - sample_variant: %s) -- SAMPLE BEST_MODEL NOT GOOD" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        elif (hfreq == False and bfreq == False) or (hfreq == False and sfreq == False):
            warningtext= " (checkWhySampleComparisonImpossible - sample_variant: %s) -- COMPARISON NOT NEEDED: ONE SAMPLE HAS NO FREQUENCIES AND NO HAPMAP" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        elif bmeta == True and smeta == False:
            warningtext= " (checkWhySampleComparisonImpossible - sample_variant: %s) -- COMPARISON NOT NEEDED: BEST IS THE META-ANALYSIS OF SAMPLE" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
        else:
            criticaltext = " %s (compareSamples) -- CANNOT PERFORM COMPARISON: WHY? -- BEST/SAMPLE: {\"%s\":%s,\"%s\":%s}" % (textbase,b['sample_id'],json.dumps(b, sort_keys= True, indent = 4),s['sample_id'],json.dumps(s, sort_keys= True, indent = 4))
            logging.critical("%s" %(criticaltext))
            sys.exit("%s" %(criticaltext))
    ### END OF checkWhySampleComparisonImpossible 
    
    ### checkWhyBestModel ###
    # warning on why new sample is best
    # Called by compareSamples
    ###
    def checkWhyBestSample(self,b,s,hfreq ,sfreq,bfreq,bmeta,smeta,sv ): 
        textbase = 'best :%s, sample: %s' % (b['sample_id'],s['sample_id'])
        logging.debug(' Function: checkWhyBestSample: %s' % textbase)
        
        if b['best_model']['significant'] == False:
            warningtext= " (checkWhyBestSample - sample_variant: %s) -- BEST_SAMPLE BEST_MODEL IS NOT SIGNIFICANT" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        elif b['best_model']['OR_OK'] == False:
            warningtext= " (checkWhyBestSample - sample_variant: %s) -- BEST_SAMPLE BEST_MODEL HAS NO OR DATA" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        elif b['best_model']['is_OK'] == False:
            warningtext= " (checkWhyBestSample - sample_variant: %s) -- BEST_SAMPLE BEST_MODEL IS NOT GOOD" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        elif hfreq == False and sfreq == True and bfreq == False:
            warningtext= " (checkWhyBestSample - sample_variant: %s) -- BEST_SAMPLE BEST_MODEL HAS NO FREQUENCIES" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        elif smeta == True :
            warningtext= " (checkWhyBestSample - sample_variant: %s) -- SAMPLE IS THE META-ANALYSIS OF BEST_SAMPLE" %sv
            logging.warning(" %s %s" % (textbase,warningtext))
            self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],warningtext)
            
        else:
            criticaltext = " %s (compareSamples) -- SAMPLE BETTER THAN BEST: WHY? -- BEST/SAMPLE: {\"%s\":%s,\"%s\":%s}" % (textbase,b['sample_id'],s['sample_id'],b['sample_id'],json.dumps(b, sort_keys= True, indent = 4),s['sample_id'],json.dumps(s, sort_keys= True, indent = 4))                       
            logging.critical("%s" %(criticaltext))
            sys.exit("%s" %(criticaltext))
    ### END OF checkWhyBestModel 
    
    ### isSampleMetaAnalysis ###
    # is this sample_id a meta-analysis
    # Called by compare_sample
    ###  
    def isSampleMetaAnalysis(self,sid, meta):
        logging.debug(' Function: isSampleMetaAnalysis: sample_id :%s, meta_list: %s' % (sid, meta))
        is_meta = False

        if self.convertToNone(meta) is not None:
            for m in meta:
                if str(sid) == str(m):
                    logging.info(' Meta sample: found sample_id: %s in meta_sample list: %s' % (sid, meta))
                    is_meta = True
        return is_meta
    ### END OF isSampleMetaAnalysis 

    ### convertToInteger ###
    # convert a value into integer
    # called by :
    # - isSampleMetaAnalysis
    # - update_meta
    ###
    def convertToInteger(self, value):
        logging.debug(' Function: convertToNone: value %s' % (value))
        
        ### none or null -> none
        if str(value).lower() == 'none' or str(value).lower() == 'null' or str(value).lower() == 'na':
            return None
        
        ## -1 instead of None
        else:
            try:
                if  int(value) == -1:
                    return None  
               
                ##  error if not -1
                else:
                    logging.info(" value: %s (convertToNone) -- found data (number)"  % value)  
                    return value
            ##  error if not an interger or None like value
            except:
                logging.warning(" value: %s (convertToNone) -- found data (string)"  % value)
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],'(convertToNone) -- found data (string)')
                return value
    ### END OF convert_2_integer

    ### convert_2_None ###
    # convert a value into None
    # called by :
    # - is_meta
    # - update_meta
    ###
    def convertToNone(self, value):
        logging.debug(' Function: convertToNone: value %s' % (value))
        
        ### none or null -> none
        if str(value).lower() == 'none' or str(value).lower() == 'null' or str(value).lower() == 'na':
            return None
        
        ## -1 instead of None
        else:
            try:
                if  int(value) == -1:
                    return None  
               
                ##  error if not -1
                else:
                    logging.debug(" value: %s (convertToNone) -- found data (number)"  % value)  
                    return value
            ##  error if not an interger or None like value
            except:
                logging.warning(" value: %s (convertToNone) -- found data (string)"  % value)
                self.fillWarningSummary(STATISTICS['shared_functions']['WARNING'],'(convertToNone) -- found data (string)')
                return value
    ### END OF convert_2_None
      
    ### test_sample_score ###
    # find best sample from sample specific scores
    # called by compare_samples
    ### 
    def findBestScoreOf2Samples(self, b, s):
        logging.debug(' Function: findBestScoreOf2Samples: best: %s, sample %s' % (b['sample_id'],s['sample_id']))

        ## Initialize test and dif tables
        b['test'] = []
        b['dif'] = []
        s['test'] = []
        s['dif'] = []
        
        ## test sample specific values
        b,s = self.test('study_type', b,s, "max", 1)
        b,s = self.test('sample_size', b,s, "max", 2)
        b,s = self.test('sign_sample', b,s, "max", 1) 
        b,s = self.test('nsign_sample', b,s, "min", 1)
        
        return b, s    
    ### END OF findBestScoreOf2Samples  
                               
    ### test_best_model ###
    # - compare best model of each sample
    # - update test and dif tables for these samples
    # Called by compare_sample
    ### 
    def findBestModelOf2Samples(self,sample1, sample2,hapfreq,same_variant):
        logging.debug(' Function: findBestModelOf2Samples: best: %s, sample %s, hapfreq %s' % (sample1['sample_id'], sample2['sample_id'],hapfreq))
        
        ### Initialization
        temp_best_model = copy.deepcopy(sample1['best_model'])
        temp_sample_model = copy.deepcopy(sample2['best_model'])
        temp_best_model['sample_id'] = sample1['sample_id']
        temp_sample_model['sample_id'] = sample2['sample_id']
######### TODO CHANGE SO THAT SAMPLE_ID HAVE RS NAME IN IT?
        temp_best_model = self.compareModels(temp_best_model,temp_sample_model, hapfreq,False)
        
        if temp_best_model == None :
            temp_best_model = sample1

        ### best model wins
        if temp_best_model['sample_id'] == sample1['sample_id']:
            logging.info(" best: %s, sample: %s (findBestModelOf2Samples) -- best sample has best_model" % (sample1['sample_id'], sample2['sample_id']))
            sample1['test'].append(1)
            sample2['test'].append(0)
            sample1['dif'].append(1)
            sample2['dif'].append(-1)   
            
        ### sample model wins
        elif temp_best_model['sample_id'] == sample2['sample_id']:
            logging.info(" best: %s, sample: %s (findBestModelOf2Samples) -- new sample has best_model" % (sample1['sample_id'], sample2['sample_id']))
            sample1['test'].append(0)
            sample2['test'].append(1)
            sample1['dif'].append(-1)
            sample2['dif'].append(1) 
        
        ### neither samples are best SHOULD NEVER HAPPEN
        else:
            criticaltext = " BEST: %s, SAMPLE: %s (findBestModelOf2Samples) -- I CANNOT FIND THE SAMPLE WITH BEST_MODEL (MODEL: %s, BEST/SAMPLE: \"%s\":%s,\"%s\":%s)" % (sample1['sample_id'], sample2['sample_id'], temp_best_model['sample_id'],sample1['sample_id'],json.dumps(sample1, sort_keys= True, indent = 4),sample2['sample_id'],json.dumps(sample2, sort_keys= True, indent = 4))
            logging.critical(criticaltext)
            sys.exit(criticaltext)
            
        return sample1, sample2
    ### END OF test_best_model
    
    ### fillWarningSummary ###
    def fillWarningSummary(self,STAT,comment):
        #logging.debug(' Function: fillWarningSummary: comment %s' % (comment))
        try:
            STAT[comment] += 1
        except:
            STAT[comment] = 1
    ### END OF fill_warning_summary
    

