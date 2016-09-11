#!/usr/bin/env python
"""
   Adjusted odds ratio calculations
   connects to genetic engine 
   09/15/12 - 0.0.1 : 
"""


"""
   Algorithm & Flow
   1. add required pieces (confidence interval, etc.)
   2. check confidence, mask to 1 if includes 1
   3. form the equations
   4. solve the non-linear equations
   
   and log properly
"""
from GeneticEngine import AbstractVariantOrAdjuster
from SolveEquations import SolveEquations
from RelativePosition import RelativePosition
import pandas as pd  #@UnresolvedImport
pd.set_printoptions(max_colwidth=200)
import logging
import copy, math
import AlgoGeneConfig as config
from AlgoGeneUtil import AlgoGeneUtil


class VariantOrAdjuster(AbstractVariantOrAdjuster):
    
    def __init__(self,disease_name, lifetime_prevalence, selected_variant):
        logging.debug(' Function: __init__ of VariantOrAdjuster class for disease %s', disease_name)
        '''
        self.disease = disease_name
        self.selected_variant     = selected_variant
        self.variant_adjusted_or = {}
        '''
        self.util = AlgoGeneUtil()
        self.pmid = []
                
        super(VariantOrAdjuster, self).__init__(disease_name,lifetime_prevalence, selected_variant)
        self.selected_variant_df  = None   # selected_variant frame with information from qc results

    def readDataFrame(self):
        ''' read the json into a dataframe '''
        logging.debug(' Function: readDataFrame')
        logging.info(' reading QCed selected_variant ...')
        tmp = []
        for ethn, val in self.selected_variant.items():
            for gene, val2 in val.items():
                if val2['sample_best_OR'].lower() == "haplotypic":
                    print "its a HAP! %s " % val
                else:
                    tmp.append({'ethn':ethn, 'gene':gene, 'var':val2['LD_best_variant'],
                            'sample': val2['sample_best_OR'], 'bioset':val2['variant_best_sample']})    
        
        self.qc_data_df = pd.DataFrame(tmp)
        logging.info(' reading QCed selected_variant Done. Dataframe size is : %s' % str(self.qc_data_df.shape) )
    
    ''' finds information for the variant '''
    def getVariantInformation(self, var, bioset, ethn):
        logging.debug(' Function: getVariantInformation - var: %s, bioset: %s, ethn: %s' %(var, bioset,ethn))
        query = {'disease_id':self.disease,'dbsnp':var, 'bioset_id' : bioset , 'broad_ethnicity':ethn}
        crsr  = config.NEXTBIO_DB.nextbio.cooked.find( query)
        
        # if alleles
        if crsr.count() == 2:
            for dic in crsr:
                self.snp_info = copy.deepcopy(dic)
                if dic['or_value'] != 1:
                    risk_or     = float(dic['or_value'])
                    risk_allele = dic["allele"]
                    if dic["frequency_pop"]:
                        f2 = float(dic["frequency_pop"])
                    elif dic["frequency_control"]:
                        f2 = float(dic["frequency_control"])
                    else:
                        f2 = None
                else:
                    ref_or     = dic['or_value']
                    ref_allele = dic["allele"]
                    if dic["frequency_pop"]:
                        f1 = float(dic["frequency_pop"])
                    elif dic["frequency_control"]:
                        f1 = float(dic["frequency_control"])
                    else:
                        f1 = None
            logging.info("before going to hardy-weinberg : ")
            logging.info("ref_allele: " + ref_allele + " risk_allele: " + risk_allele + " ref or: "+ str(ref_or) + " risk or: " + str(risk_or) + " Model : " + dic["model_type"] )
            var_info = self.calculateHardyWeinbergEquilibrium(ref_allele, ref_or, risk_allele, risk_or,  dic["model_type"] , f1, f2)
                        
        ### if genotypes take OR directly
        elif crsr.count() == 3 or crsr.count() == 5:
            genotype={}
            for dic in crsr:
                if len(dic["allele"]) == 2 :
                    if dic["or_value"] == 1 and not genotype.has_key('ref_gt'):
                        genotype['ref_or'] = dic["or_value"]
                        genotype['ref_gt'] = dic["allele"]
                        if dic["frequency_pop"]:
                            f1 = float(dic["frequency_pop"])
                        elif dic["frequency_control"]:
                            f1 = float(dic["frequency_control"])
                        else:
                            f1 = None
                        genotype['ref_freq'] = f1
                    
                    elif not genotype.has_key('gt2_gt'):
                        genotype['gt2_or'] = float(dic["or_value"])
                        genotype['gt2_gt'] = dic["allele"]
                        if dic["frequency_pop"]:
                            f2 = float(dic["frequency_pop"])
                        elif dic["frequency_control"]:
                            f2 = float(dic["frequency_control"])
                        else:
                            f2 = None
                        genotype['gt2_freq']=f2
                        
                    elif not genotype.has_key('gt3_gt'):
                        genotype['gt3_or'] = float(dic["or_value"])
                        genotype['gt3_gt'] = dic["allele"]
                        if dic["frequency_pop"]:
                            f3 = float(dic["frequency_pop"])
                        elif dic["frequency_control"]:
                            f3 = float(dic["frequency_control"])
                        else:
                            f3 = None
                        genotype['gt3_freq'] = f3
            var_info = genotype
        else:
            message = "encountered a puzzling situation, dictionary recieved has %s  elements" % crsr.count()
            raise Exception(message)
    
        # if no frequency found, go to hapmap
        if not any([var_info['gt3_freq'],var_info['gt2_freq'], var_info["ref_freq"] ] ) :
            hap_freqs = self.getHapMapFrequencies(var_info["ref_gt"], var_info["gt2_freq"], var_info["gt3_freq"], var, ethn ,dic["narrow_ethnicity"] )  
            var_info["ref_freq"] = hap_freqs[0]
            var_info["gt2_freq"] = hap_freqs[1]
            var_info["gt3_freq"] = hap_freqs[2]
        return var_info
    
    
    def getHapMapFrequencies(self, gt1, gt2, gt3, snp, ethn, nethn):
        ''' fetch frequencies from hapmap '''
        logging.debug(' Function: getHapMapFrequencies -- variant: %s, ethnicity: %s, narrow: %s' % (snp, ethn,nethn))
        maja  = self.snp_info ["hapmap_major_allele"]
        mina  = self.snp_info ["hapmap_minor_allele"]
        ###### TODO ADD CHECK WITH OTHER MINOR/MAJOR IF HERE

        ### get hapmap frequencies
        if self.snp_info["hapmap_genotype_freq"] is not None:
            logging.info(' (getHapMapFrequencies) -- I have HapMap Freq' )
            words = self.snp_info ["hapmap_genotype_freq"].split('|')
            freqs = {}
            for index in range(0,len(words)):
                if '/' in words[index]:
                    gt = words[index].replace('/','')
                    freqs[gt] = float(words[index+1])
            freqs_order = [0,0,0]
            for gt in freqs:
                if gt1 == gt:
                    freqs_order[0] = float(freqs[gt])
                elif maja+maja == gt or gt == mina+mina:
                    freqs_order[2] = float(freqs[gt])
                else:
                    freqs_order[1] = float(freqs[gt])
            return freqs_order
        elif self.snp_info["hapmap_allele_freq"] is not None:
            self.util.warnMe('warning',' calculating HapMap gt_freq from allele_freq(getHapMapFrequencies) - variant: %s, ethnicity: %s, narrow: %s' % (snp, ethn,nethn))
            freqs_order = [0,0,0]
            if float(self.snp_info["hapmap_allele_freq"]) <0.5:
                max_freq = 1-float(self.snp_info["hapmap_allele_freq"])
                min_freq = float(self.snp_info["hapmap_allele_freq"])
            else:
                max_freq = float(self.snp_info["hapmap_allele_freq"])
                min_freq = 1-float(self.snp_info["hapmap_allele_freq"])
            freqs_order[0] = min_freq*min_freq
            freqs_order[2] = max_freq*max_freq
            freqs_order[1] = 2*min_freq*max_freq
            return freqs_order
        else:
            self.util.warnMe('warning',' I DO NOT have HapMap Freq (getHapMapFrequencies) -- variant: %s, ethnicity: %s, narrow: %s' % (snp, ethn,nethn))
            return None
        

        
    def calculateHardyWeinbergEquilibrium(self, ref_a1, or1, a2, or2, model, f1, f2):
        ''' to calculate genotype frequency from allele frequencies
        ref_a1 : reference allele '''
        logging.debug(' Function: calculateHardyWeinbergEquilibrium - bioset: %s, snp: %s' % (self.bioset, self.snp))
        genotype = {}
        genotype['ref_gt'] = ref_a1 + ref_a1
        genotype['gt2_gt'] = a2 + a2
        genotype['gt3_gt'] = ref_a1 + a2
        ### calculate genotype frequencies if possible
        if f1 is not None or f2 is not None:
            if f1 is not None and f2 is  None:
                f2 = 1 - f1
            elif f1 is None and f2 is not None:
                f1 = 1-f2
            genotype["gt2_freq"] = f2*f2
            genotype["gt3_freq"] = 2*f1*f2
            genotype["ref_freq"] = f1*f1
        else:
            genotype["gt2_freq"] = None
            genotype["gt3_freq"] = None
            genotype["ref_freq"] = None
        ###### end of Celine added

        if model.lower() == 'additive':
            genotype['ref_or'] = or1**2
            genotype['gt2_or'] = or2**2
            genotype['gt3_or'] = 2 * or1 * or2           
        elif model.lower() == 'dominant':
            genotype['ref_or'] = or1
            genotype['gt2_or'] = or2**2
            genotype['gt3_or'] = or2           
        elif model.lower() == 'recessive':
            genotype['ref_or'] = or1
            genotype['gt2_or'] = or2
            genotype['gt3_or'] = or1   
        else:
            message = 'Model cant be recognized : %s - bioset: %s, snp: %s' % ( model,self.bioset, self.snp)
            raise(message)
        return genotype
            
    
    def calculateAdjustedOr(self):
        ''' main module tha glues it all together '''
        logging.debug(' Function: calculateAdjustedOr')
        # read the qc results in our frame
        self.readDataFrame()
        var_info = []
        ns = 0
        for idx, row in self.qc_data_df.iterrows():
            try:
                s=None; tt={};
                self.snp = row['var']
                self.bioset = row['bioset']
                s  = self.getVariantInformation(row['var'],row['bioset'], row['ethn'] )
                tt = self.qc_data_df.ix[idx].to_dict()
                tt.update(s)          
                tt['pmid'] = self.snp_info['pmid']
                tt['chromosome'] = self.snp_info['chromosome']
                tt['position'] = self.snp_info['position_37']
                tt['position36'] = self.snp_info['position_36']
                var_info.append(tt)
                ns+=1
                if ns%50 == 0:
                    print ns, self.snp, self.bioset
            except Exception, e:
                print str(e)
        self.var_info = pd.DataFrame(var_info)
        logging.info("passing everything to non-linear solver")
        # pass it to equation solver:
        solutions = SolveEquations(self.var_info, self.lifetime_prevalence)
        self.reformating(solutions.variants_data_frame.T.to_dict())

    ### addGenotype
    def addGenotype(self, genotype, doc, data):
        logging.debug(' Function: addGenotype - bioset: %s, snp: %s' % (self.bioset, self.snp))
        newdoc_gt = copy.deepcopy(doc)
        newdoc_gt['gt'] = data['%s_gt' % genotype]
        newdoc_gt['oldor'] = data['%s_or' % genotype]
        newdoc_gt['expor'] = data['%s_adj' % genotype]
        newdoc_gt['freq'] = data['%s_freq' % genotype]
        newdoc_gt['or'] = math.log(data['%s_adj' % genotype])
        return newdoc_gt
    ### END addGenotype
        
    ### format the adjusted Or data for input for the next function
    def reformating(self,solutions ):
        logging.debug(' Function: reformating')
        for data in solutions:
            if solutions[data]["ref_adj"] > 0 and solutions[data]["gt2_adj"] >0 and solutions[data]["gt3_adj"] >0:
                self.snp = solutions[data]['var']
                self.bioset = solutions[data]['bioset']
                logging.info(' (reformating) - bioset: %s, snp: %s' % (self.bioset, self.snp))
                newdoc = {}
                newdoc['disease'] = self.disease
                newdoc['snp'] = solutions[data]['var']
                newdoc['study_in'] = solutions[data]['ethn']
                newdoc['ethnicity'] = solutions[data]['ethn']
                newdoc['gene'] = solutions[data]['gene']
                newdoc['pmid'] = solutions[data]['pmid']
                if newdoc['pmid'] not in self.pmid:
                    self.pmid.append(newdoc['pmid'])
                newdoc['chromosome'] = solutions[data]['chromosome']
                newdoc['position'] = int(solutions[data]['position'])
                newdoc['position36'] =  int(solutions[data]['position36'])
                relative_position = RelativePosition()
                newdoc['relative_position'] = relative_position.get_relative_position(newdoc['chromosome'], newdoc['position'] ) 
                newdoc['bioset'] = solutions[data]['bioset']
                gt = 'ref'
                newdoc_gt = self.addGenotype(gt, newdoc,solutions[data])
                self.variant_adjusted_or.append(newdoc_gt)
                
                gt = 'gt2'
                newdoc_gt = self.addGenotype(gt, newdoc,solutions[data])
                self.variant_adjusted_or.append(newdoc_gt)
                
                gt = 'gt3'
                newdoc_gt = self.addGenotype(gt, newdoc,solutions[data])
                self.variant_adjusted_or.append(newdoc_gt)  
            else:
                self.snp = solutions[data]['var']
                self.bioset = solutions[data]['bioset']
                self.util.warnMe('warning' ,' NOT ADDED - ref: %s, gt2: %s, gt3: %s (reformating) - bioset: %s, snp: %s' % (solutions[data]["ref_adj"] , solutions[data]["gt2_adj"] ,  solutions[data]["gt3_adj"] ,self.bioset, self.snp))    
    ### END reformating
                
        
