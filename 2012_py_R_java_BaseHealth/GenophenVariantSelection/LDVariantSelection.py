#!/usr/bin/python

'''
Created on Jul 30, 2012
logs added 8/9/12
@author: celine
'''

import copy
import json
import logging

### functions and global variable from SharedFunctions
import SharedFunctions as SF
STAT = SF.STATISTICS
SF_class = SF.SharedFunctions()

### internal global variable
PER_ETHNICITY = ''

class LDVariantSelection(object):

    def __init__(self):
        '''
        Constructor
        '''
 
    ### initialize_LD_summaries ###
    def initialize_LD_summaries(self, STAT, num):
        logging.debug(' Function: initialize_LD_summaries')
        if num == 0:
            STAT['WARNING'] = {} #IN
            STAT['ERROR'] = {} #IN
            STAT['ethnicity_with_variant'] =[] # IN
            STAT['ethnicity_no_variant'] =[] # IN
            STAT['total_ethnicity'] = 0 # IN
        
        STAT['total_LD_block'] = 0 #IN
        STAT['LD_block_with_variant'] = 0#IN
        STAT['LD_block_no_variant'] = 0#IN
        STAT['OK_variant'] = 0#IN
        STAT['not_OK_variant'] = 0#IN
        STAT['total_variant'] = 0#IN
        STAT['total_sign_variant'] = 0#IN
        STAT['total_nsign_variant'] = 0#IN
            
    ### END OF initialize_LD_summaries

    ### LD_variant_selection ###
    # selects the best variants within an LD block
    # called by main
    ### 
    def LD_variant_selection(self,variant_score,variant_info, LD_blocks):
        logging.info('\n================================\nFunction: LD_variant_selection\n================================')
        
        ## Init output
        LD_selection = {} 
        STAT['LD_selection'] = {}
        self.initialize_LD_summaries(STAT['LD_selection'], 0)
        
        ### loop on ethnicities -> choose best variants per LD bocks
        for et in LD_blocks: 
            ## Initi ethnicity
            LD_selection[et] = {} 
            global PER_ETHNICITY 
            PER_ETHNICITY = et
            STAT['LD_selection']['total_ethnicity'] +=1 
            
            STAT['LD_selection'][PER_ETHNICITY] = {}
            self.initialize_LD_summaries(STAT['LD_selection'][PER_ETHNICITY], 1)
            
            ### loop on LD blocks
            for LD in LD_blocks[et]: 
                logging.info(" *****: ethnicity: %s, LD block: %s (LD_variant_selection)-- start variant selection *****" % (et,LD))
                STAT['LD_selection'][PER_ETHNICITY]['total_LD_block'] += 1
                ## check ethnicity ID
                #from SharedFunctions import SharedFunctions
                SF_class.check_id(LD,LD_blocks[et][LD]['LD_block_id'])
                
                ### Init LD selection output
                selection = self.LD_selection_initiation(LD_blocks[et][LD],variant_score,variant_info, et)
                
                ### loop on variant data in this block - choose best variant
                for variant in selection['variant_score']: 
                    logging.info(" ***** ethnicity: %s, LD block: %s, variant: %s (LD_variant_selection) -- start compare_variants function" % (et,LD, variant))
                    selection['LD_best_variant'] =  self.compare_variants(selection['LD_best_variant'],selection['variant_score'][variant] )
                    
                
                ###  update LD_selection 
                LD_selection[et] = self.update_LD_selection(LD,LD_selection[et], selection)
        logging.info(" ***** I AM DONE SELECTING THE VARIANTS FOR ALL THE LD BLOCKS *****")
        return LD_selection
    ### END OF LD_variant_selection
    
    ### LD_selection_initiation ###
    # copy variant scores within LD
    # add all intermediate output/ cooked data for this variant
    # called by LD_variant_selection
    ###
    def LD_selection_initiation(self,LD_block, v_score, v_info, ethnicity):
        logging.debug(' Function: LD_selection_initiation: block %s' % (LD_block['LD_block_id']))
        
        ### Init LD block selection
        LD_select = copy.deepcopy(LD_block)
    
        LD_select['LD_best_variant'] = {}
        LD_select['variant_score'] = {}
        LD_select['num_ok_variant'] = 0
        
        ### loop on variants in this block
        for variant in LD_select["variant_list"]: 
            logging.debug("variant, LD_selection_initiation: variant: %s" % (variant))
            STAT['LD_selection'][PER_ETHNICITY]['total_variant'] += 1

            ## Check if variant is an option
            OKvariant = self.check_LD_variant(v_score, v_info, ethnicity, variant)
            if OKvariant == True:      
                logging.info(" LD block: %s, variant (initiation): %s -- fit for comparison" % (LD_block['LD_block_id'],variant))       
                STAT['LD_selection'][PER_ETHNICITY]['OK_variant'] += 1
                   
                ## record variant_score
                LD_select['num_ok_variant'] += 1
                LD_select['variant_score'][variant] = copy.deepcopy(v_score[variant]) 
                
                ## record hapmap freq and nexbio scores
                LD_select['variant_score'][variant]['nextbio_score'] =  v_info[variant]['nextbio_score']     
                LD_select['variant_score'][variant]['nextbio_pvalue'] =  v_info[variant]['nextbio_pvalue']  
                LD_select['variant_score'][variant]['hapmap_freq'] = v_info[variant]['broad_ethnicity'][ethnicity]
    
                ##  extract only this ethnicity_score
                LD_select['variant_score'][variant]['ethnicity_score'] = copy.deepcopy(v_score[variant]['ethnicity_score'][ethnicity])
                LD_select['variant_score'][variant]['ethnicity_score']['variant_id'] = variant
            else:
                logging.warning(" LD block: %s, variant: %s (initiation) -- UNFIT VARIANT FOR COMPARISON" % (LD_block['LD_block_id'],variant))    
                SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(initiation) -- UNFIT VARIANT FOR COMPARISON')
                STAT['LD_selection'][PER_ETHNICITY]['not_OK_variant'] += 1
                    
        ### initialise best variant                    
        if len(LD_select['variant_score']) > 0:
            logging.info(" LD block: %s, variant: %s (initiation) -- Initial best variant" % (LD_block['LD_block_id'],LD_select['variant_score'][LD_select['variant_score'].keys()[0]]['variant_id']))
            LD_select['LD_best_variant'] =  LD_select['variant_score'][LD_select['variant_score'].keys()[0]] 
            STAT['LD_selection'][PER_ETHNICITY]['LD_block_with_variant'] += 1

        else: 
            logging.warning(" LD block: %s (initiation) -- CANNOT FIND ANY VARIANT" % LD_block['LD_block_id'])
            SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(initiation) -- UNFIT VARIANT FOR COMPARISON')
            STAT['LD_selection'][PER_ETHNICITY]['LD_block_no_variant'] += 1
        ### delete list of variant for this block
        del LD_select['variant_list']
        
        return LD_select
    ### END OF LD_selection_initiation
    
    ### LD_variant_selection ###
    # check that the OR is sign, ok and there is frequencies
    # called by main
    ### 
    def check_LD_variant(self,vscore, vinfo, ethnicity, vid):
        logging.debug(' Function: check_LD_variant: ethnicity: %s, variant: %s' % (ethnicity, vid))

        ### find variant id in score and information
        OKv = False
        if vid in vscore and vid in vinfo:
            ## find ethnicity in score and LD info
            if ethnicity in vscore[vid]['ethnicity_score'] and ethnicity in vinfo[vid]['broad_ethnicity']:
                ## Initialization
                hapfreq =  vinfo[vid]['broad_ethnicity'][ethnicity]
                model = vscore[vid]['ethnicity_score'][ethnicity]['best_sample']['best_model']
                mfreq = model['control_freq']

                ## sign, OK and frequencies found
                if (hapfreq == True or (hapfreq == False and mfreq == True)) and\
                    model['significant'] == True and model['is_OK'] == True and model['OR_OK'] == True and model['risk_allele_OK'] == True:
                    logging.info(" ethnicity: %s, variant: %s (check_LD_variant) -- variant OK for comparison" % (ethnicity,vid))
                    OKv = True 
                    STAT['LD_selection'][PER_ETHNICITY]['total_sign_variant'] += 1
                ## can't use this variant
                else:
                    STAT['LD_selection'][PER_ETHNICITY]['total_sign_variant'] += 1
                    if hapfreq == False and mfreq == False:
                        logging.warning(" ethnicity: %s, variant: %s (check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (no frequencies)" % (ethnicity,vid))
                        SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (no frequencies)')
                    elif model['significant'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (not significant)" % (ethnicity,vid))
                        SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (not significant)')
                        STAT['LD_selection'][PER_ETHNICITY]['total_nsign_variant'] += 1
                        STAT['LD_selection'][PER_ETHNICITY]['total_sign_variant'] -= 1
                    elif  model['is_OK'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (model error)" % (ethnicity,vid))
                        SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (model error)')
                    elif  model['risk_allele_OK'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (risk allele error)" % (ethnicity,vid))
                        SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (risk allele error)')
                    elif model['OR_OK'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (no OR data found)" % (ethnicity,vid))
                        SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- VARIANT UNFIT FOR COMPARISON (no OR data found)')
                    else :
                        errortext = '(check_LD_variant) -- VARIANT UNFIT FOR COMPARISON. WHY?'
                        logging.error(" ethnicity: %s, variant: %s %s (SCORE/INFO: {\"data\":%s,\"info\":%s})" % (ethnicity,vid,errortext,json.dumps(vscore[vid],sort_keys= True, indent = 4),json.dumps(vinfo[vid], sort_keys= True, indent = 4)))
                        SF_class.fill_warning_summary(STAT['LD_selection']['ERROR'],errortext)

            ### can't find ethnicity 
            else:
                logging.warning("  ethnicity: %s, variant: %s (check_LD_variant) -- COULD NOT FIND DATA FOR VARIANT IN THIS ETHNICITY" % (ethnicity, vid))
                SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- COULD NOT FIND DATA FOR VARIANT IN THIS ETHNICITY')
        
        ### can't find variant id in variant_score or LD info data
        else:
            logging.warning(" variant: %s (check_LD_variant) -- COULD NOT FIND SCORE OR INFORMATION FOR THIS VARIANT" % (vid))
            SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(check_LD_variant) -- COULD NOT FIND SCORE OR INFORMATION FOR THIS VARIANT')
        
        return OKv
    ### END OF check_LD_variant    
    
    ### compare_variants ###
    # choose the best out of 2 variants
    # both variant are significants and have frequencies (checked at initialization)
    # Uses
    # - best variant scores
    # - best ethnicity scores
    # - best sample scores
    # called by LD_variant_selection
    ### 
    def compare_variants(self,best,variant):
        logging.debug(' Function: compare_variants: best: %s, variant: %s' % (best['variant_id'], variant['variant_id']))
    
        best_id = []
        ### compare only different variants
        if best['variant_id'] != variant['variant_id']:
            logging.info(" best: %s, variant: %s (compare_variants) -- compare variants" % (best['variant_id'], variant['variant_id']))
            ###  find best variant based on variant-specific scores
            best_id.append(self.find_best_variant_from_variant_score(best,variant))
    
            ### compare ethnicity_score
            best_id.append(self.find_best_variant_from_ethnicity_score(best, variant))
    
            ### compare samples
            best_id.append(self.find_best_variant_from_sample_score(best, variant))
    
            ### find best variant based on 3 comparisons
            best = self.choose_best_variant(best_id, best, variant)
        else:
            logging.info(" best: %s, variant: %s (compare_variants) -- no comparison needed, best = variant" % (best['variant_id'], variant['variant_id']))
        return best
    ### END OF compare_variants
    
    ### find_best_variant_score ###
    # find best variant from variant specific scores
    # called by compare_variants
    ### 
    def find_best_variant_from_variant_score(self, b, v):
        logging.debug(' Function: find_best_variant_from_variant_score: best: %s, variant: %s' % (b['variant_id'], v['variant_id']))
        
        ### Initialize test and dif tables
        b['test'] = []
        b['dif'] = []
        v['test'] = []
        v['dif'] = []
    
        ### test sample specific values
        b,v = SF_class.test('nextbio_score', b,v, "max", 1)
        b,v = SF_class.test('nextbio_pvalue', b,v, "max", 1 )
        b,v = SF_class.test('variant_nsign_ethnicity', b,v, "min", 1)
        b,v = SF_class.test('variant_sign_ethnicity', b,v, "max", 1)
        
        ###  compare risk alleles if >1 choose the other
        b,v = self.test_risk_allele('variant_risk_allele', b,v)
    
        ###  find best variant best of variant-specific scores
        logging.info(' best: %s, variant: %s (find_best_variant_from_variant_score)' % (b['variant_id'], v['variant_id']))
        temp = SF_class.find_best(b,v)
        if temp == None:
            return None
        else:
            return temp['variant_id']
        
    ### END OF find_best_variant_from_variant_score
    
    ### test_risk_allele ###
    # test which variant as most consistant risk alleles
    # all result to test and dif
    # called by find_best_variant_score
    ### 
    def test_risk_allele(self, score, b, v):
        logging.info(" score: %s, best_RA %s, new_RA: %s (test_risk_allele)" % (score,b[score],v[score]))
        
        ### define level
        level = score.split('_risk_allele')[0]
        
        ### compare number of risk alleles
        b['num_allele'] = len(b[score]) 
        v['num_allele'] = len(v[score])
        b , v = SF_class.test('num_allele', b,v, "min", 1)
        
        ### compare counts of risk alleles
        if level == 'ethnicity':
            bra = b['best_sample']['best_model']['risk_allele']
            vra = v['best_sample']['best_model']['risk_allele']
        elif level == 'variant':
            bra = b['ethnicity_score']['best_sample']['best_model']['risk_allele']
            vra = v['ethnicity_score']['best_sample']['best_model']['risk_allele']
            
        b['count_allele'] = b[score][bra]
        v['count_allele'] = v[score][vra]
        b , v = SF_class.test('count_allele', b,v, "max", 1)
        
        return b, v
    ### END OF test_risk_allele
    
    ### find_best_variant_from_ethnicity_score ###
    # get the id of variant with best ethnicity score
    # called by compare_variants
    ### 
    def find_best_variant_from_ethnicity_score(self,best, variant):
        logging.debug(' Function: find_best_variant_from_ethnicity_score: best: %s, variant: %s' % (best['variant_id'], variant['variant_id']))
        
        ### Initialization
        bethscore = best['ethnicity_score']
        vethscore = variant['ethnicity_score']
        
        bethscore['test'] = []
        bethscore['dif'] = []
        vethscore['test'] = []
        vethscore['dif'] = []
        
        ### test on ethnicity_score
        bethscore,vethscore = SF_class.test('ethnicity_sign_sample', bethscore,vethscore, "max", 1)
        bethscore,vethscore = SF_class.test('ethnicity_nsign_sample', bethscore,vethscore, "min", 1)
        bethscore,vethscore = SF_class.test('mean_study_type', bethscore,vethscore, "max", 1)
        bethscore,vethscore = SF_class.test('mean_sample_size', bethscore,vethscore, "max", 1)
        bethscore,vethscore = SF_class.test('mean_error', bethscore,vethscore, "min", 1)
        bethscore,vethscore = SF_class.test('mean_pvalue', bethscore,vethscore, "min", 1)
        bethscore,vethscore = SF_class.test('mean_OR', bethscore,vethscore, "max", 1)
    
        ###  compare risk alleles if >1 choose the other
        bethscore,vethscore = self.test_risk_allele('ethnicity_risk_allele', bethscore,vethscore)
    
        ### find best variant based on ethnicity_score
        logging.info(' best: %s, variant: %s (find_best_variant_from_ethnicity_score)' % (best['variant_id'], variant['variant_id']))
        temp = SF_class.find_best(bethscore,vethscore) 
        if temp == None:
            return None
        else:
            return temp['variant_id']
    
    ### END OF find_best_variant_from_ethnicity_score
    
    ### find_best_variant_from_sample_score ###
    # get the variant with best model OR
    # called by compare_variants
    ### 
    def find_best_variant_from_sample_score(self,best, variant):
        logging.debug(' Function: find_best_variant_from_sample_score: best%s, variant: %s' % (best['variant_id'], variant['variant_id']))
 
        ### initialization
        bsample = best['ethnicity_score']['best_sample']
        vsample = variant['ethnicity_score']['best_sample']
        bsample['variant_id'] = best['variant_id']
        vsample['variant_id'] = variant['variant_id']
        
        ### compare samples
        logging.info(' best: %s, variant: %s (find_best_variant_from_sample_score)' % (best['variant_id'], variant['variant_id']))
        temp =  SF_class.compare_samples(bsample, vsample, True, False)
        if temp == None:
            return None
        else:
            return temp['variant_id']
    ### END OF find_best_variant_from_sample_score
      
    ### choose_best_variant ###
    # find variant that was best the most out of 3 comparisons
    # called by compare_variants
    ### 
    def choose_best_variant(self,b_id, b, v):
        logging.debug(' Function: choose_best_variant: list: %s best: %s, variant: %s' % (b_id,b['variant_id'], v['variant_id']))

        ## Initialization
        b['test'] = []
        b['dif'] = []
        v['test'] = []
        v['dif'] = []
        b['best'] = 0
        v['best'] = 0
        
        ### count best_variant
        for rs in b_id:
            if rs == b['variant_id']:
                logging.info(" best: %s, variant: %s (choose_best_variant) -- best variant still best" % (b['variant_id'], v['variant_id']))
                b['best'] += 1
            elif rs == v['variant_id']:
                logging.info(" best: %s, variant: %s (choose_best_variant) -- new variant is best" % (b['variant_id'], v['variant_id']))
                v['best'] += 1
            
            elif rs == None:
                logging.warning(" best: %s, variant: %s (choose_best_variant) -- VARIANT IS NONE" % (b['variant_id'], v['variant_id']))
                SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(choose_best_variant) -- VARIANT IS NONE')
        
            else:
                logging.error(" UNKNOWN VARIANT: %s, best: %s, variant: %s (choose_best_variant) -- CAN NOT FIND AMOUNGST BEST AND VARIANT" % (rs,b['variant_id'], v['variant_id']))
                SF_class.fill_warning_summary(STAT['LD_selection']['ERROR'],'(choose_best_variant) -- CAN NOT FIND AMOUNGST BEST AND VARIANT')
        
        ### find best variant     
        b,v = SF_class.test('best',b,v,'max',1)  
        if len(b['test']) >1: 
            logging.info(" best: %s, variant: %s (choose_best_variant) -- use find best" % (b['variant_id'], v['variant_id']))
            b = SF_class.find_best(b,v)
        else:
            logging.info(" best: %s, variant: %s (choose_best_variant) -- use first non null variant" % (b['variant_id'], v['variant_id']))
            for rs in b_id:
                if rs != None:
                    if rs == b['variant_id']:
                        return b
                    elif rs == v['variant_id']:
                        return v
                    

        return b
    ### END OF find_best_variant_score
    
    ### update_LD_selection ###
    # save best variants for this LD block
    # called by LD_variant_selection
    ###
    def update_LD_selection(self,LD,selection_et, selection):
        logging.debug(' Function: update_LD_selection: LD block  %s' % (LD))

        if len(selection['LD_best_variant']) > 0:
            ## update list of ethnicities
            if PER_ETHNICITY not in STAT['LD_selection']['ethnicity_with_variant']:
                STAT['LD_selection']['ethnicity_with_variant'].append(PER_ETHNICITY)
            STAT['LD_selection']['LD_block_with_variant'] += 1
            logging.info(" LD block: %s, variant: %s (update) -- best variant" % (LD,selection['LD_best_variant']['variant_id']))
            selection_et[LD]= {}
            selection_et[LD]['LD_best_variant'] = selection['LD_best_variant']['variant_id']
            selection_et[LD]['variant_best_sample'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['sample_id']
            selection_et[LD]['sample_best_OR'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['model_id']
            selection_et[LD]['max_error'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['max_error']
            selection_et[LD]['max_pvalue'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['max_pvalue'] 
            ### output replication factor
            total_samples = selection['LD_best_variant']['ethnicity_score']['ethnicity_sign_sample'] + selection['LD_best_variant']['ethnicity_score']['ethnicity_nsign_sample']
            selection_et[LD]['significant_samples'] = selection['LD_best_variant']['ethnicity_score']['ethnicity_sign_sample']
            selection_et[LD]['total_samples'] = total_samples
        else :
            logging.warning(" LD block: %s (update) -- COULD NOT FIND ANY VARIANT" % (LD))
            SF_class.fill_warning_summary(STAT['LD_selection']['WARNING'],'(update) -- COULD NOT FIND ANY VARIANT')
            if PER_ETHNICITY not in STAT['LD_selection']['ethnicity_no_variant']:
                STAT['LD_selection']['ethnicity_no_variant'].append(PER_ETHNICITY)
            STAT['LD_selection']['LD_block_no_variant'] += 1
        return selection_et
    ### END OF update_LD_selection