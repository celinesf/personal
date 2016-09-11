#!/usr/bin/python

'''
Created on Jul 30, 2012
logs added 8/9/12
@author: celine
'''

import copy
import json
import logging
import AlgoGeneConfig as config
from AlgoGeneUtil import AlgoGeneUtil
import CalculateGenophenConfidenceScore as CalculateGenophenConfidenceScore
CONF_SCORE = CalculateGenophenConfidenceScore.CalculateGenophenConfidenceScore()

### functions and global variable from VariantSelectorSharedFunctions
import VariantSelectorSharedFunctions as SF
STAT = SF.STATISTICS
SharedFunctions = SF.VariantSelectorSharedFunctions()


class SelectVariantPerLdBlock(object):

    def __init__(self):
        ### internal global variable
        self.ethnicity = ''
        self.util = AlgoGeneUtil()

 
    ### initializeLdBlockSummaries ###
    # called by selectVariantPerLdBlock
    def initializeLdBlockSummaries(self, STAT, num):
        logging.debug(' Function: initializeLdBlockSummaries')
        if num == 0:
            STAT['WARNING'] = {} 
            STAT['ERROR'] = {} 
            STAT['ethnicity_with_variant'] =[] 
            STAT['ethnicity_no_variant'] =[] 
            STAT['total_ethnicity'] = 0 
        
        STAT['total_LD_block'] = 0 
        STAT['LD_block_with_variant'] = 0
        STAT['LD_block_no_variant'] = 0
        STAT['OK_variant'] = 0
        STAT['not_OK_variant'] = 0
        STAT['total_variant'] = 0
        STAT['total_sign_variant'] = 0
        STAT['total_nsign_variant'] = 0
    ### END OF initializeLdBlockSummaries

    ### selectVariantPerLdBlock ###
    # selects the best variants within an LD block
    # called by main
    ### 
    def selectVariantPerLdBlock(self,variant_score,variant_info, LD_blocks):
        logging.info('\n================================\nFunction: selectVariantPerLdBlock\n================================')
        
        ## Init output
        LD_selection = {} 
        disease_genophen_score = 0
        STAT['LD_selection'] = {}
        self.initializeLdBlockSummaries(STAT['LD_selection'], 0)
      
        self.treshold = {'Caucasian': max(config.POP_LIMIT,min(STAT['variant_scoring']['max_sample_size'] / 2, config.POP_LIMIT*20)), 'Asian': max(config.POP_LIMIT,min(STAT['variant_scoring']['max_sample_size'] / 2, config.POP_LIMIT*10))} 
        print self.treshold 
        logging.info(" self.treshold: %s (selectVariantPerLdBlock)" % self.treshold)
        
        ### loop on ethnicities -> choose best variants per LD bocks
        for et in LD_blocks: 
            if et != 'mixed':
                ## Initi ethnicity
                LD_selection[et] = {} 
                
                self.ethnicity = et
                STAT['LD_selection']['total_ethnicity'] +=1 
                
                STAT['LD_selection'][self.ethnicity] = {}
                self.initializeLdBlockSummaries(STAT['LD_selection'][self.ethnicity], 1)
                
                ### loop on LD blocks
                for LD in LD_blocks[et]: 
                    logging.info(" *****: ethnicity: %s, LD block: %s (selectVariantPerLdBlock)-- start variant selection *****" % (et,LD))
                    STAT['LD_selection'][self.ethnicity]['total_LD_block'] += 1
                    ## check ethnicity ID
                    #from VariantSelectorSharedFunctions import VariantSelectorSharedFunctions
                    SharedFunctions.checkIdFit(LD,LD_blocks[et][LD]['LD_block_id'])
                    
                    ### Init LD selection output
                    selection = self.initializeLdBlockSelection(LD_blocks[et][LD],variant_score,variant_info, et)
                    
                    ### loop on variant data in this block - choose best variant
                    for variant in selection['variant_score']: 
                        logging.info(" ***** ethnicity: %s, LD block: %s, variant: %s (selectVariantPerLdBlock) -- start compareVariants function" % (et,LD, variant))
                        selection['LD_best_variant'] =  self.compareVariants(selection['LD_best_variant'],selection['variant_score'][variant] )
          
                    ### check LD block had enough evidence for this disease
                    selection = self.checkLdBlockScore(selection)
                    ###  update LD_selection 
                    LD_selection[et] = self.updateLdBlockSelection(LD,LD_selection[et], selection)
        LD_selection = self.updateEthnicitySelection(LD_selection)
        num = 0
        for snp in CONF_SCORE.genophen_score:
            num += 1
            disease_genophen_score += CONF_SCORE.genophen_score[snp]
        if num >0:
            disease_genophen_score /= num * 4.5/5
        self.util.warnMe('critical', 'disease_genophen_score: %s' % disease_genophen_score)
        ####################### TO ADD - update disease info with genophen genetic score

        logging.info(" ***** I AM DONE SELECTING THE VARIANTS FOR ALL THE LD BLOCKS *****")
        return LD_selection, disease_genophen_score
    ### END OF selectVariantPerLdBlock

    ### updateEthnicitySelection ###
    # check>1 SNP for the ethnicity
    ### 
    def updateEthnicitySelection(self, data):
        logging.debug(' Function: updateEthnicitySelection: block %s' % (data))
        new_data = copy.deepcopy(data)
        for et in data:
            if len(data[et]) <2:
                self.util.warnMe('warning', ' %s DONT HAVE ENOUGH SNP (updateEthnicitySelection) ' % et)
                new_data.pop(et)
            else:
                self.util.warnMe('warning',' %s I have enough SNP: %s (updateEthnicitySelection) ' % (et, len(data[et])))
        return new_data
    
    ### checkLdBlockScore ###
    # check that the OR is sign, ok and there is frequencies
    # called by selectVariantPerLdBlock
    ### 
    def checkLdBlockScore(self, block):
        logging.debug(' Function: checkLdBlockScore: block %s' % (block['LD_block_id']))
        block['is_ok'] = False
        if len(block['LD_best_variant']) > 0:
            logging.debug( 'var %s, var %s/%s=%s, ethn %s/%s=%s, sample %s/%s=%s' % (block['LD_block_id'], block['sign_variant'] , block['nsign_variant'],block['num_variant'] , block['sign_ethnicity'], block['nsign_ethnicity'] , block['num_ethnicity'] , block['sign_sample'], block['nsign_sample'] ,block['num_sample'] )) 
            
            sign = block['sign_variant']+  block['sign_ethnicity'] + block['sign_sample']
            nsign = block['nsign_variant'] + block['nsign_ethnicity'] + block['nsign_sample']
            total = block['num_variant'] + block['num_ethnicity'] + block['num_sample']
            dif = [block['sign_variant']-block['nsign_variant'], block['sign_ethnicity']-block['nsign_ethnicity'],block['sign_sample']-block['nsign_sample']]
            text = " - ss: %s, sign: %s nsign: %s, tot: %s dif: %s - var %s/%s=%s, ethn %s/%s=%s, sample %s/%s=%s" % (block['sample_size'],sign, nsign, total, dif,  block['sign_variant'] , block['nsign_variant'],block['num_variant'] , block['sign_ethnicity'], block['nsign_ethnicity'] , block['num_ethnicity'] , block['sign_sample'], block['nsign_sample'] ,block['num_sample'])

            if  dif[0] > 0 and dif[1] > 0 and dif[2] > 0 and sign > config.SIGN and total > config.SIGN :
                logging.info(" block %s (checkLdBlockScore) -- good block %s" % (block['LD_block_id'], text))
                block['is_ok'] = True
            elif dif[0] > 0 and sum(dif) > 0 and dif[2] > 0 and sign >= config.SIGN and total >= config.SIGN:# dicount non significant ethnicity if still a good snp
                if self.ethnicity not in self.treshold:
                    limit = 1000
                else:
                    limit = self.treshold[self.ethnicity]
                if block['sample_size'] > limit and block['nsign_ethnicity'] == 0 : # rescue by sample size
                    logging.warning(" LD block: %s (checkLdBlockScore) -- RESCUED BLOCK sample size %s" % (block['LD_block_id'], text))
                    SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- RESCUED BLOCK sample size')
                    block['is_ok'] = True
                elif block['num_ethnicity'] > 1:
                    sign = block['sign_variant']+ block['sign_sample']
                    nsign = block['nsign_variant'] + block['nsign_sample']
                    total = block['num_variant'] + block['num_sample']
                    dif = sign - nsign
                    text = " - ss: %s, sign: %s nsign: %s, tot: %s dif: %s - var %s/%s=%s, ethn %s/%s=%s, sample %s/%s=%s" % (block['sample_size'],sign, nsign, total, dif,  block['sign_variant'] , block['nsign_variant'],block['num_variant'] , block['sign_ethnicity'], block['nsign_ethnicity'] , block['num_ethnicity'] , block['sign_sample'], block['nsign_sample'] ,block['num_sample'])

                    if dif > 0 and sign > (config.SIGN -1) and total > (config.SIGN -1) :
                        logging.warning(" LD block: %s (checkLdBlockScore) -- RESCUED BLOCK REMOVE ETH%s" % (block['LD_block_id'], text))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- RESCUED BLOCK REMOVE ETH ')
                        block['is_ok'] = True
                    elif dif > 0 and sign >= (config.SIGN -1) and total >= (config.SIGN -1): # rescue by sample size
                        if block['sample_size'] > limit  :
                            logging.warning(" LD block: %s (checkLdBlockScore) -- RESCUED BLOCK sample size ETH%s" % (block['LD_block_id'], text))
                            SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- RESCUED BLOCK sample size ETH')
                            block['is_ok'] = True
                        else:
                            logging.critical(' block %s (checkLdBlockScore) -- BAD BLOCK ETH%s' % (block['LD_block_id'], text))
                            SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- BAD BLOCK ETH')
                    else:
                        logging.critical(' block %s (checkLdBlockScore) -- BAD BLOCK MORE ETH%s' % (block['LD_block_id'], text))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- BAD BLOCK MORE ETH')
                else:
                    logging.critical(' block %s (checkLdBlockScore) -- BAD BLOCK%s' % (block['LD_block_id'], text))
                    SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- BAD BLOCK')
            else :
                logging.critical(' block %s (checkLdBlockScore) -- NOT SIGN REP%s' % (block['LD_block_id'], text))
                SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- NOT SIGN REP ')
        else:
            logging.warning(" LD block: %s (checkLdBlockScore) -- COULD NOT FIND ANY VARIANT" % (block['LD_block_id']))
            SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkLdBlockScore) -- COULD NOT FIND ANY VARIANT')
            if self.ethnicity not in STAT['LD_selection']['ethnicity_no_variant']:
                STAT['LD_selection']['ethnicity_no_variant'].append(self.ethnicity)
            STAT['LD_selection']['LD_block_no_variant'] += 1
        return block
    ### END OF checkLdBlockScore
    
    ### initializeLdBlockSelection ###
    # copy variant scores within LD
    # add all intermediate output/ cooked data for this variant
    # called by selectVariantPerLdBlock
    ###
    def initializeLdBlockSelection(self,LD_block, v_score, v_info, ethnicity):
        logging.debug(' Function: initializeLdBlockSelection: block %s' % (LD_block['LD_block_id']))
        
        ### Init LD block selection
        LD_select = copy.deepcopy(LD_block)
        LD_select['LD_best_variant'] = {}
        LD_select['variant_score'] = {}
        LD_select['num_ok_variant'] = 0
        LD_select['sign_variant'] = 0
        LD_select['nsign_variant'] = 0
        LD_select['num_variant'] = 0
        LD_select['sign_ethnicity'] = 0
        LD_select['nsign_ethnicity'] = 0
        LD_select['num_ethnicity'] = 0
        LD_select['sign_sample'] = 0
        LD_select['nsign_sample'] = 0
        LD_select['num_sample'] = 0
        LD_select['sample_size'] = 0
        
        ### loop on variants in this block
        for variant in LD_select["variant_list"]: 
            logging.debug("variant, initializeLdBlockSelection: variant: %s" % (variant))
            STAT['LD_selection'][self.ethnicity]['total_variant'] += 1

            ## Check if variant is an option
            LD_select,OKvariant = self.checkVariantOkForLdBlock(LD_select,v_score, v_info, ethnicity, variant)
            if OKvariant == True:      
                logging.info(" LD block: %s, variant (initiation): %s -- fit for comparison" % (LD_block['LD_block_id'],variant))       
                STAT['LD_selection'][self.ethnicity]['OK_variant'] += 1
                   
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
                SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(initiation) -- UNFIT VARIANT FOR COMPARISON')
                STAT['LD_selection'][self.ethnicity]['not_OK_variant'] += 1
                    
        ### initialise best variant                    
        if len(LD_select['variant_score']) > 0:
            logging.info(" LD block: %s, variant: %s (initiation) -- Initial best variant" % (LD_block['LD_block_id'],LD_select['variant_score'][LD_select['variant_score'].keys()[0]]['variant_id']))
            LD_select['LD_best_variant'] =  LD_select['variant_score'][LD_select['variant_score'].keys()[0]] 
            STAT['LD_selection'][self.ethnicity]['LD_block_with_variant'] += 1

        else: 
            logging.warning(" LD block: %s (initiation) -- CANNOT FIND ANY VARIANT" % LD_block['LD_block_id'])
            SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(initiation) -- CANNOT FIND ANY VARIANT')
            STAT['LD_selection'][self.ethnicity]['LD_block_no_variant'] += 1
        ### delete list of variant for this block
        del LD_select['variant_list']
        
        return LD_select
    ### END OF initializeLdBlockSelection
    
    ### checkVariantOkForLdBlock ###
    # check that the OR is sign, ok and there is frequencies
    # called by main
    ### 
    def checkVariantOkForLdBlock(self,block, vscore, vinfo, ethnicity, vid):
        logging.debug(' Function: checkVariantOkForLdBlock: ethnicity: %s, variant: %s' % (ethnicity, vid))

        ### find variant id in score and information
        OKv = False
        if vid in vscore and vid in vinfo:
            ## block score
            block['num_variant'] += 1
            block['sign_ethnicity'] += vscore[vid]['variant_sign_ethnicity']
            block['nsign_ethnicity'] += vscore[vid]['variant_nsign_ethnicity']
            block['num_ethnicity'] += vscore[vid]['variant_num_ethnicity']
            block['sign_sample'] += vscore[vid]['ethnicity_score'][ethnicity]['ethnicity_sign_sample']
            block['nsign_sample'] += vscore[vid]['ethnicity_score'][ethnicity]['ethnicity_nsign_sample']
            block['num_sample'] += vscore[vid]['ethnicity_score'][ethnicity]['ethnicity_num_sample']
            block['sample_size'] = max(block['sample_size'] , vscore[vid]['ethnicity_score'][ethnicity]['best_sample']['sample_size'])

            ## find ethnicity in score and LD info
            if ethnicity in vscore[vid]['ethnicity_score'] and ethnicity in vinfo[vid]['broad_ethnicity']:
                ## Initialization
                hapfreq =  vinfo[vid]['broad_ethnicity'][ethnicity]
                model = vscore[vid]['ethnicity_score'][ethnicity]['best_sample']['best_model']
                mfreq = model['control_freq']

                ## sign, OK and frequencies found
                if (hapfreq == True or (hapfreq == False and mfreq == True)) and\
                    model['significant'] == True and model['is_OK'] == True and model['OR_OK'] == True and model['risk_allele_OK'] == True:
                    logging.info(" ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- variant OK for comparison" % (ethnicity,vid))
                    OKv = True 
                    STAT['LD_selection'][self.ethnicity]['total_sign_variant'] += 1

                    ## getLdBlockScore for variant good for comparison
                    block['sign_variant'] += 1
                    
                ## can't use this variant
                else:
                    block = self.getLdBlockScore(block,vscore[vid],ethnicity )
                    STAT['LD_selection'][self.ethnicity]['total_sign_variant'] += 1
                    if hapfreq == False and mfreq == False:
                        logging.warning(" ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (no frequencies)" % (ethnicity,vid))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (no frequencies)')
                    elif model['significant'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (not significant)" % (ethnicity,vid))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (not significant)')
                        STAT['LD_selection'][self.ethnicity]['total_nsign_variant'] += 1
                        STAT['LD_selection'][self.ethnicity]['total_sign_variant'] -= 1
                    elif  model['is_OK'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (model error)" % (ethnicity,vid))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (model error)')
                    elif  model['risk_allele_OK'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (risk allele error)" % (ethnicity,vid))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (risk allele error)')
                    elif model['OR_OK'] == False:
                        logging.warning(" ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (no OR data found)" % (ethnicity,vid))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON (no OR data found)')
                    else :
                        errortext = '(checkVariantOkForLdBlock) -- VARIANT UNFIT FOR COMPARISON. WHY?'
                        logging.error(" ethnicity: %s, variant: %s %s (SCORE/INFO: {\"data\":%s,\"info\":%s})" % (ethnicity,vid,errortext,json.dumps(vscore[vid],sort_keys= True, indent = 4),json.dumps(vinfo[vid], sort_keys= True, indent = 4)))
                        SharedFunctions.fillWarningSummary(STAT['LD_selection']['ERROR'],errortext)

            ### can't find ethnicity 
            else:
                logging.warning("  ethnicity: %s, variant: %s (checkVariantOkForLdBlock) -- COULD NOT FIND DATA FOR VARIANT IN THIS ETHNICITY" % (ethnicity, vid))
                SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- COULD NOT FIND DATA FOR VARIANT IN THIS ETHNICITY')
        
        ### can't find variant id in variant_score or LD info data
        else:
            logging.warning(" variant: %s (checkVariantOkForLdBlock) -- COULD NOT FIND SCORE OR INFORMATION FOR THIS VARIANT" % (vid))
            SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(checkVariantOkForLdBlock) -- COULD NOT FIND SCORE OR INFORMATION FOR THIS VARIANT')
        
        return block , OKv
    ### END OF checkVariantOkForLdBlock    
   
    ### getLdBlockScore ###
    # calculate block score when variant is not good for comparison
    # called by checkVariantOkForLdBlock
    ### 
    def getLdBlockScore(self,block,variant, ethnicity):
        logging.debug(' Function: getLdBlockScore: variant: %s, ethhnicity: %s' % (variant['variant_id'], ethnicity))
        model = variant['ethnicity_score'][ethnicity]['best_sample']['best_model']
        if  model['is_OK'] == True : #and model['OR_OK'] == True:
            if model['significant'] == True and model['risk_allele_OK'] == True :
                block['sign_variant'] += 1
            else:
                block['nsign_variant'] += 1
        else:
            block['nsign_variant'] += 1 
        return block
    ### END OF getLdBlockScore    
    
    ### compareVariants ###
    # choose the best out of 2 variants
    # both variant are significants and have frequencies (checked at initialization)
    # Uses
    # - best variant scores
    # - best ethnicity scores
    # - best sample scores
    # called by selectVariantPerLdBlock
    ### 
    def compareVariants(self,best,variant):
        logging.debug(' Function: compareVariants: best: %s, variant: %s' % (best['variant_id'], variant['variant_id']))
    
        best_id = []
        ### compare only different variants
        if best['variant_id'] != variant['variant_id']:
            logging.info(" best: %s, variant: %s (compareVariants) -- compare variants" % (best['variant_id'], variant['variant_id']))
            ###  find best variant based on variant-specific scores
            best_id.append(self.findBestVariantFromVariantScore(best,variant))
    
            ### compare ethnicity_score
            best_id.append(self.findBestVariantFromEthnicityScore(best, variant))
    
            ### compare samples
            best_id.append(self.findBestVariantFromSampleScore(best, variant))
    
            ### find best variant based on 3 comparisons
            best = self.chooseBestVariantForLdBlock(best_id, best, variant)
        else:
            logging.info(" best: %s, variant: %s (compareVariants) -- no comparison needed, best = variant" % (best['variant_id'], variant['variant_id']))
        return best
    ### END OF compareVariants
    
    ### find_best_variant_score ###
    # find best variant from variant specific scores
    # called by compareVariants
    ### 
    def findBestVariantFromVariantScore(self, b, v):
        logging.debug(' Function: findBestVariantFromVariantScore: best: %s, variant: %s' % (b['variant_id'], v['variant_id']))
        
        ### Initialize test and dif tables
        b['test'] = []
        b['dif'] = []
        v['test'] = []
        v['dif'] = []

        ### test sample specific values
        b,v = SharedFunctions.test('nextbio_score', b,v, "max", 1)
        b,v = SharedFunctions.test('nextbio_pvalue', b,v, "max", 1 )
        b,v = SharedFunctions.test('variant_sign_ethnicity', b,v, "max", 2)
        b,v = SharedFunctions.test('variant_nsign_ethnicity', b,v, "min", .5)
        b,v = SharedFunctions.test('variant_num_ethnicity', b,v, "max", 5)
        
        ###  compare risk alleles if >1 choose the other
        b,v = self.testRiskAllele('variant_ethnicity_risk_allele', b,v)
        b,v = self.testRiskAllele('variant_total_risk_allele', b,v)
    
        ###  find best variant best of variant-specific scores
        logging.info(' best: %s, variant: %s (findBestVariantFromVariantScore)' % (b['variant_id'], v['variant_id']))
        temp = SharedFunctions.findBest(b,v,'variant_id')
        if temp == None:
            return None
        else:
            return temp['variant_id']     
    ### END OF findBestVariantFromVariantScore
    
    ### testRiskAllele ###
    # test which variant as most consistant risk alleles
    # all result to test and dif
    # called by find_best_variant_score
    ### 
    def testRiskAllele(self, score, b, v):
        logging.info(" score: %s, best_RA %s, new_RA: %s (testRiskAllele)" % (score,b[score],v[score]))
        
        ### define level
        level = score.split('_risk_allele')[0]

        ### compare number of risk alleles
        b['num_allele'] = len(b[score]) 
        v['num_allele'] = len(v[score])
        b , v = SharedFunctions.test('num_allele', b,v, "min", 1)
        
        ### compare counts of risk alleles
        if level == 'ethnicity':
            bra = b['best_sample']['best_model']['risk_allele']
            vra = v['best_sample']['best_model']['risk_allele']
        elif level == 'variant_ethnicity' or level == 'variant_total' :
            bra = b['ethnicity_score']['best_sample']['best_model']['risk_allele']
            vra = v['ethnicity_score']['best_sample']['best_model']['risk_allele']
          
        b['count_allele'] = b[score][bra]
        v['count_allele'] = v[score][vra]
        b , v = SharedFunctions.test('count_allele', b,v, "max", 1)
        
        return b, v
    ### END OF testRiskAllele
    
    ### findBestVariantFromEthnicityScore ###
    # get the id of variant with best ethnicity score
    # called by compareVariants
    ### 
    def findBestVariantFromEthnicityScore(self,best, variant):
        logging.debug(' Function: findBestVariantFromEthnicityScore: best: %s, variant: %s' % (best['variant_id'], variant['variant_id']))
        
        ### Initialization
        bethscore = best['ethnicity_score']
        vethscore = variant['ethnicity_score']
        
        bethscore['test'] = []
        bethscore['dif'] = []
        vethscore['test'] = []
        vethscore['dif'] = []
        
        ### test on ethnicity_score
        bethscore,vethscore = SharedFunctions.test('ethnicity_sign_sample', bethscore,vethscore, "max", 1)
        bethscore,vethscore = SharedFunctions.test('ethnicity_nsign_sample', bethscore,vethscore, "min", 1)
        bethscore,vethscore = SharedFunctions.test('ethnicity_num_sample', bethscore,vethscore, "max", 5)
        
        bethscore,vethscore = SharedFunctions.test('mean_study_type', bethscore,vethscore, "max", 1)
        bethscore,vethscore = SharedFunctions.test('mean_sample_size', bethscore,vethscore, "max", 2)
        bethscore,vethscore = SharedFunctions.test('mean_error', bethscore,vethscore, "min", 1.5)
        bethscore,vethscore = SharedFunctions.test('mean_pvalue', bethscore,vethscore, "min", 1)
        bethscore,vethscore = SharedFunctions.test('mean_OR', bethscore,vethscore, "max", .5)
    
        ###  compare risk alleles if >1 choose the other
        bethscore,vethscore = self.testRiskAllele('ethnicity_risk_allele', bethscore,vethscore)
    
        ### find best variant based on ethnicity_score
        logging.info(' best: %s, variant: %s (findBestVariantFromEthnicityScore)' % (best['variant_id'], variant['variant_id']))
        temp = SharedFunctions.findBest(bethscore,vethscore,'variant_id') 
        if temp == None:
            return None
        else:
            return temp['variant_id']
    
    ### END OF findBestVariantFromEthnicityScore
    
    ### findBestVariantFromSampleScore ###
    # get the variant with best model OR
    # called by compareVariants
    ### 
    def findBestVariantFromSampleScore(self,best, variant):
        logging.debug(' Function: findBestVariantFromSampleScore: best%s, variant: %s' % (best['variant_id'], variant['variant_id']))
 
        ### initialization
        bsample = best['ethnicity_score']['best_sample']
        vsample = variant['ethnicity_score']['best_sample']
        bsample['variant_id'] = best['variant_id']
        vsample['variant_id'] = variant['variant_id']
        
        ### compare samples
        logging.info(' best: %s, variant: %s (findBestVariantFromSampleScore)' % (best['variant_id'], variant['variant_id']))
        temp =  SharedFunctions.compareSamples(bsample, vsample, True, False)
        if temp == None:
            return None
        else:
            return temp['variant_id']
    ### END OF findBestVariantFromSampleScore
      
    ### chooseBestVariantForLdBlock ###
    # find variant that was best the most out of 3 comparisons
    # called by compareVariants
    ### 
    def chooseBestVariantForLdBlock(self,b_id, b, v):
        logging.debug(' Function: chooseBestVariantForLdBlock: list: %s best: %s, variant: %s' % (b_id,b['variant_id'], v['variant_id']))

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
                logging.info(" best: %s, variant: %s (chooseBestVariantForLdBlock) -- best variant still best" % (b['variant_id'], v['variant_id']))
                b['best'] += 1
            elif rs == v['variant_id']:
                logging.info(" best: %s, variant: %s (chooseBestVariantForLdBlock) -- new variant is best" % (b['variant_id'], v['variant_id']))
                v['best'] += 1
            
            elif rs == None:
                logging.warning(" best: %s, variant: %s (chooseBestVariantForLdBlock) -- VARIANT IS NONE" % (b['variant_id'], v['variant_id']))
                SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(chooseBestVariantForLdBlock) -- VARIANT IS NONE')
        
            else:
                logging.error(" UNKNOWN VARIANT: %s, best: %s, variant: %s (chooseBestVariantForLdBlock) -- CAN NOT FIND AMOUNGST BEST AND VARIANT" % (rs,b['variant_id'], v['variant_id']))
                SharedFunctions.fillWarningSummary(STAT['LD_selection']['ERROR'],'(chooseBestVariantForLdBlock) -- CAN NOT FIND AMOUNGST BEST AND VARIANT')
        
        ### find best variant    
        b,v = SharedFunctions.test('best',b,v,'max',1)  

        if len(b['test']) > 0: 
            logging.info(" best: %s, variant: %s (chooseBestVariantForLdBlock) -- use find best" % (b['variant_id'], v['variant_id']))
            b = SharedFunctions.findBest(b,v,'variant_id')
        else:
            logging.info(" best: %s, variant: %s (chooseBestVariantForLdBlock) -- use first non null variant" % (b['variant_id'], v['variant_id']))
            for rs in b_id:
                if rs != None:
                    if rs == b['variant_id']:
                        return b
                    elif rs == v['variant_id']:
                        return v
                    

        return b
    ### END OF find_best_variant_score
    
    ### updateLdBlockSelection ###
    # save best variants for this LD block
    # called by selectVariantPerLdBlock
    ###
    def updateLdBlockSelection(self,LD,selection_et, selection):
        logging.debug(' Function: updateLdBlockSelection: LD block  %s' % (LD))

        #if len(selection['LD_best_variant']) > 0 :
        if selection['is_ok'] == True:
            ## update list of ethnicities
            if self.ethnicity not in STAT['LD_selection']['ethnicity_with_variant']:
                STAT['LD_selection']['ethnicity_with_variant'].append(self.ethnicity)
            STAT['LD_selection']['LD_block_with_variant'] += 1
            logging.info(" LD block: %s, variant: %s (update) -- best variant" % (LD,selection['LD_best_variant']['variant_id']))
            selection_et[LD]= {}
            selection_et[LD]['LD_best_variant'] = selection['LD_best_variant']['variant_id']
            selection_et[LD]['variant_best_sample'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['sample_id'].split('_')[1]
            selection_et[LD]['sample_best_OR'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['model_id'].split('_')[2]
            selection_et[LD]['max_error'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['max_error']
            selection_et[LD]['max_pvalue'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['max_pvalue'] 
            selection_et[LD]['min_OR'] = selection['LD_best_variant']['ethnicity_score']['best_sample']['best_model']['min_OR'] 
            ### output replication factors
            selection_et[LD]['significant_variants'] = selection['sign_variant']
            selection_et[LD]['total_variants'] = selection['num_variant']            
            total_samples = selection['LD_best_variant']['ethnicity_score']['ethnicity_sign_sample'] + selection['LD_best_variant']['ethnicity_score']['ethnicity_nsign_sample']
            selection_et[LD]['significant_samples'] = selection['LD_best_variant']['ethnicity_score']['ethnicity_sign_sample']
            selection_et[LD]['total_samples'] = total_samples
            CONF_SCORE.calculateGenophenConfidenceScore(selection)            
        else :
            logging.warning(" LD block: %s (update) -- NOT GOOD ENOUGH" % (LD))
            SharedFunctions.fillWarningSummary(STAT['LD_selection']['WARNING'],'(update) -- NOT GOOD ENOUGH')
            if self.ethnicity not in STAT['LD_selection']['ethnicity_no_variant']:
                STAT['LD_selection']['ethnicity_no_variant'].append(self.ethnicity)
            STAT['LD_selection']['LD_block_no_variant'] += 1
        return selection_et
    ### END OF updateLdBlockSelection