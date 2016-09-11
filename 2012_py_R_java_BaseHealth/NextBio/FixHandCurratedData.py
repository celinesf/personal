#!/usr/bin/env python

"""
   Fix error of what allele info are for in hand currated data
   06/17/13: version 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) Get map for SNP table snp_table_header
    2) Recover data from excel file
    3) fix data
"""

import logging,copy, json,  xlwt
import NextBioConfig as config
from NextBioUtils import NextBioUtils
from SnpUtils import SnpUtils
from ExcelUtils import ExcelUtils

class FixHandCurratedData():
    def __init__(self,disease_name):
        self.disease_name = disease_name
        self.column_key ={}
        self.util = NextBioUtils() 
        self.snp_util = SnpUtils()
        self.excel_util = ExcelUtils()
        self.snp_issue = {} ### record snps with specific issue to recall for other biosets
        self.cooked_data=[] ## record fixed daa
        self.tofix ={} ### record bioset reference and bioset they fix
        self.frequencies = None


    ''' isMinorAllele1 '''
    def isMinorAllele1(self, snp_data):
        logging.debug(' Function:  isMinorAllele1 - snp:%s, bioset: %s' % (snp_data['snp'], self.bioset_id))
        is_minor = []
        for allele in range(1,3):
            is_minor.append(self.snp_util.isSameAllele(snp_data['minor_allele'], snp_data['allele_%s'%allele]))
        if sum(is_minor) == 0:
            is_minor = []
            snp_data['minor_allele'] = self.snp_util.reverseAlleles(snp_data['minor_allele'])
            if 'major_allele' in snp_data:
                snp_data['major_allele'] = self.snp_util.reverseAlleles(snp_data['major_allele'])
            for allele in range(1,3):
                is_minor.append(self.snp_util.isSameAllele(snp_data['minor_allele'], snp_data['allele_%s'%allele]))
            
        if sum(is_minor) == 0 or sum(is_minor) >1 :
            self.util.warnMe('warning', ' FWD AND REV CAN"T FIT MINOR ALLELE (checkIfDataAllele1OrMinor)- sum(is_minor):%s data:%s' % (sum(is_minor),snp_data))
        else:
            ### Allele 1 = minor -> ok no change
            if is_minor[0] == 1:
                logging.info(' allele 1 is minor (checkIfDataAllele1OrMinor) - snp:%s' % snp_data['snp'])
            else:
                logging.warning( ' ALLELE_2 IS MINOR (checkIfDataAllele1OrMinor)- is_minor:%s data:%s' % (is_minor,snp_data))
                if snp_data['snp'] not in self.snp_issue:
                    if 'major_allele' in snp_data:
                        self.snp_issue[snp_data['snp']]={'minor_allele':snp_data['minor_allele'],'major_allele':snp_data['major_allele'],'issue':['major_allele', 'minor_allele','_odds_ratio', '95%ci','_frequency_case','_frequency_control']}
                    else: 
                        self.snp_issue[snp_data['snp']]={'minor_allele':snp_data['minor_allele'],'major_allele':None,'issue':['major_allele', 'minor_allele','_odds_ratio', '95%ci','_frequency_case','_frequency_control']}

    ''' checkOrFrequencies '''
    def checkOrFrequencies(self, data):
        logging.debug(' Function:  checkOrFrequencies - snp:%s, bioset: %s' % (data['snp'],self.bioset_id))
        for a in range(1,3):            
            if 'allele_%s_frequency_case' % a in data :#and 'allele_%s_frequency_control' % a in data 'allele_%s_frequency_odds_ratio' % a in data:
                freq_case = (data['allele_%s_frequency_case' % a])
                freq_control = (data['allele_%s_frequency_control' % a])
                OR = data['allele_%s_odds_ratio' % a]
                pvalue = float(data['p-value'])
                if freq_case is not None and freq_control is not None and OR is not None:
                    freq_case = float(freq_case)
                    freq_control = float(freq_control)
                    OR = float(OR)
                    if OR == 1.0 and freq_case ==  freq_control:
                        logging.info(' OR=1 fit frequencies case:%s = control:%s (checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif OR == 1.0 and round(freq_case,2) == round(freq_control,2):
                        logging.warning(' ROUND OR=1 fit frequencies case:%s = control:%s (checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif OR == 1.0 and round(freq_case,2) == round(freq_control+0.005,2):
                        logging.warning(' ROUND+0.5 OR=1 fit frequencies case:%s = control:%s (checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif OR == 1.0 and round(freq_case,2) == round(freq_control-0.005,2):
                        logging.warning(' ROUND-0.5 OR=1 fit frequencies case:%s = control:%s (checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif freq_case >  freq_control and OR > 1.0:
                        logging.info(' OR>1 fit frequencies case:%s >control:%s ((checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif round(freq_case,2) >  round(freq_control,2) and OR > 1.0:
                        logging.warning(' ROUND OR>1 fit frequencies case:%s >control:%s ((checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif round(freq_case,2) >  round(freq_control,2) and OR <= 1.0 and pvalue < 0.05:
#                         print OR, freq_case, freq_control,round(freq_case,2) , round(freq_control,2), round(freq_control+0.005,2),round(freq_control-0.005,2)
                        self.snp_issue[data['snp']] = data
                        self.snp_issue[data['snp']]['issue'] = ['_odds_ratio','95%ci']
                        self.util.warnMe('critical',' OR = %s DONT fit frequencies case:%s >control:%s (checkOrFrequencies) - snp:%s, data: %s' % (float(OR),freq_case, freq_control ,data['snp'],data))
                    elif  freq_case <  freq_control and OR < 1.0:
                        logging.info(' OR<1 fit frequencies case:%s < control:%s (checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))
                    elif round(freq_case,2) <  round(freq_control,2) and OR < 1.0:
                        logging.warning(' ROUND OR<1 fit frequencies case:%s < control:%s (checkOrFrequencies) - snp:%s' % (freq_case, freq_control ,data['snp']))         
                    elif round(freq_case,2) <  round(freq_control,2) and OR >= 1.0  and pvalue < 0.05:
#                         print OR, freq_case, freq_control,round(freq_case,2) , round(freq_control,2), round(freq_control+0.005,2),round(freq_control-0.005,2)
                        self.snp_issue[data['snp']] = data
                        self.snp_issue[data['snp']]['issue'] = ['_odds_ratio','95%ci']
                        self.util.warnMe('critical',' OR = %s DONT fit frequencies case:%s < control:%s (checkOrFrequencies) - snp:%s, data: %s' % (float(OR),freq_case , freq_control ,data['snp'],data))
                    elif pvalue < 0.05:
                        self.util.warnMe('critical',' WHERE AM I?(checkOrFrequencies) - OR = %s, case:%s, control:%s, snp:%s, data: %s' % (float(OR),freq_case , freq_control ,data['snp'],data))
                    else:
                        logging.warning(' NOT SIGNIFICANT (checkOrFrequencies) - OR = %s, case:%s, control:%s, snp:%s, data: %s' % (float(OR),freq_case , freq_control ,data['snp'],data))

    ''' fixDataFromAllele2Issue '''
    def fixDataFromAllele2Issue(self, data):
        logging.debug(' Function:  fixDataFromAllele2Issue - snp:%s,bioset: %s, ' % (data['snp'], self.bioset_id))
        snp_data = copy.deepcopy(data)
        for key in data:
            flag = False
            for k in self.snp_issue[snp_data['snp']]['issue']:
                if k in key:
                    flag = True
#                     if 'fixed'.upper() not in  self.column_key[self.bioset_id]:                        
#                         self.column_key[self.bioset_id].append('fixed'.upper())
            if 'allele_1_' in key.lower() and flag == True:
                if key != 'allele_1_frequency' :
                    new_key = key.replace('_1_', '_2_')
                    snp_data[new_key] = snp_data[key]
                    snp_data[key] = None
#                     if 'fixed' not in snp_data:
#                         snp_data['fixed'] = []
#                     snp_data['fixed'].append(new_key)
                    ### add newkey to list of column
                    if new_key.upper() not in self.column_key[self.bioset_id]:                        
                        index = self.column_key[self.bioset_id].index(key.upper())
                        self.column_key[self.bioset_id].insert(index+1, new_key.upper())
                    logging.warning(' fixed keyA: %s to new_key: %s (fixDataFromAllele2Issue) - snp: %s,bioset: %s, ' % (key, new_key,data['snp'], self.bioset_id))
            if key == 'major_allele' and flag == True:
#                 if 'fixed' not in snp_data:
#                     snp_data['fixed'] = []
#                 snp_data['fixed'].append(key)
                snp_data[key] = self.snp_issue[snp_data['snp']][key]
                logging.warning(' fixed keyB: %s  (fixDataFromAllele2Issue) - snp: %s,bioset: %s, ' % (key, data['snp'], self.bioset_id))
            if key == 'minor_allele' and flag == True:
#                 if 'fixed' not in snp_data:
#                     snp_data['fixed'] = []
#                 snp_data['fixed'].append(key)
                snp_data[key] = self.snp_issue[snp_data['snp']][key]
                logging.warning(' fixed keyC: %s  (fixDataFromAllele2Issue) - snp: %s,bioset: %s, ' % (key, data['snp'], self.bioset_id))     
        self.util.warnMe('warning', ' fixed keysD: %s  (fixDataFromAllele2Issue) - snp: %s,bioset: %s, ' % (self.snp_issue[snp_data['snp']]['issue'], data['snp'], self.bioset_id))
        return copy.deepcopy(snp_data) 
   
    ''' checkGenotypeFrequencies '''
    def checkGenotypeFrequencies(self, data):
        logging.debug(' Function:  checkGenotypeFrequencies - snp:%s, bioset: %s' % (data['snp'],self.bioset_id))

        frequencies = {'case':{'total_alleles':0, 'genotype_count':[],'allele_count':[0,0], 'genotype_freq':[0,0,0],'allele_freq':[0,0]},
                    'control':{'total_alleles':0, 'genotype_count':[],'allele_count':[0,0], 'genotype_freq':[0,0,0],'allele_freq':[0,0]} }
        ### loop on genotype
        for a1 in range(1,3):
            for a2 in range(a1,3):
                for type in frequencies:
                    frequencies[type]['genotype_count'].append(float(data['g%s%s_count_%s' % (a1, a2,type)]) )
                    frequencies[type]['total_alleles'] +=   float(data['g%s%s_count_%s' % (a1, a2,type)]) *2      
                    if a1 == a2:
                        frequencies[type]['allele_count'][a1-1] += float(data['g%s%s_count_%s' % (a1, a2,type)]) *2 
                    else:
                        frequencies[type]['allele_count'][a1-1] += float(data['g%s%s_count_%s' % (a1, a2,type)]) 
                        frequencies[type]['allele_count'][a2-1] += float(data['g%s%s_count_%s' % (a1, a2,type)]) 
        for type in frequencies:
            for a1 in range(0,2):
                for a2 in range(a1,2):
                    index = a1 + a2 
                    frequencies[type]['genotype_freq'][index] = frequencies[type]['genotype_count'][index]/frequencies[type]['total_alleles']
                frequencies[type]['allele_freq'][a1] = frequencies[type]['allele_count'][a1]/frequencies[type]['total_alleles']
                if 'allele_%s_frequency_%s' % (a1+1, type) in data :
#                     self.util.warnMe('warning', ' allele_%s_frequency_%s: %s != allele_freq[%s]: %s (checkGenotypeFrequencies) - snp:%s, bioset: %s' % (a1+1, type,data['allele_%s_frequency_%s' % (a1+1, type) ], a1, frequencies[type]['allele_freq'][a1], data['snp'],self.bioset_id))
                    if round(float(data['allele_%s_frequency_%s'% (a1+1, type)]),2) != round(frequencies[type]['allele_freq'][a1],2):
                        self.util.warnMe('error', ' allele_%s_frequency_%s: %s != allele_freq[%s]: %s (checkGenotypeFrequencies) - snp:%s, bioset: %s' % (a1+1, type,data['allele_%s_frequency_%s' % (a1+1, type) ], a1, frequencies[type]['allele_freq'][a1], data['snp'],self.bioset_id))
                    else:
                        round(float(data['allele_%s_frequency_%s'% (a1+1, type)]),2) != round(frequencies[type]['allele_freq'][a1],2)
                        data['allele_%s_frequency_%s'% (a1+1, type)] = str(round(frequencies[type]['allele_freq'][a1],3))
        self.frequencies = copy.deepcopy(frequencies)
    
    ''' fdataMinorNotAllele1 
    1) check allele_1 is minor_allele
    2) if not all data are for allele_2
    3) confirm with control/case frequencies and OR value
    '''
    def checkIfDataAllele1OrMinor(self, data):
        logging.debug(' Function:  checkIfDataAllele1OrMinor - %s' %self.bioset_id)
        for snp in data['snps']:
            snp_data = copy.deepcopy(data['snps'][snp])
            ### find if allele_1 is minor
            if 'G11_COUNT_CASE'.lower() in snp_data :
                self.frequencies = None
                if snp_data['G11_COUNT_CASE'.lower()] is not None:
                    self.checkGenotypeFrequencies(snp_data)  
            if 'minor_allele' in snp_data:
                if snp_data['minor_allele'] is not None:
                    if snp_data['snp'] not in self.snp_issue:
                        self.isMinorAllele1(snp_data)
            ### frequencies fit OR
            if 'allele_1_frequency_case' in snp_data:
                self.checkOrFrequencies(snp_data)
            ### Add Allele_2 data is not
            if snp_data['snp'] in self.snp_issue:
                snp_data = self.fixDataFromAllele2Issue(snp_data)
                snp_data['fixed'] = json.dumps(self.snp_issue[snp_data['snp']]['issue'])
                if 'fixed'.upper() not in self.column_key[self.bioset_id]:                        
                    self.column_key[self.bioset_id].append('fixed'.upper())
                if json.dumps(self.snp_issue[snp_data['snp']]['issue']) not in self.bioset_info[self.bioset_id]['fixed']:
                    self.bioset_info[self.bioset_id]['fixed'].append(json.dumps(self.snp_issue[snp_data['snp']]['issue']))
            
            data['snps'][snp] = copy.deepcopy(snp_data)

    ''' addGenotypeFrequencies 
    '''
    def addGenotypeFrequencies(self, ref_bioset, bioset_tofix):
        logging.debug(' Function:  addGenotypeFrequencies- %s' %self.bioset_id)
        print 'ref_bioset',ref_bioset
        print  'bioset_tofix',bioset_tofix 
        for snp in bioset_tofix['snps']:
            if bioset_tofix['snps'][snp]['genotype_frequency'] is None and snp in ref_bioset['snps']:
                bioset_tofix['snps'][snp]['genotype_frequency'] = ref_bioset['snps'][snp]['genotype_frequency']
  
  
    ''' fixBasedOnReferenceBioset '''
    def fixBasedOnReferenceBioset(self):
        logging.debug(' Function:  fixBasedOnReferenceBioset- %s' %self.bioset_id)
        self.fixBasedOnReferenceBioset()
        ### fix using reference biosets
        for bioset in self.nextbio_data:
            if bioset in self.tofix:
                for bioset_tofix in self.tofix[bioset]:
                    if 'missing genotype -> bioset_id ' == self.tofix[bioset][bioset_tofix]:
                        self.addGenotypeFrequencies(self.nextbio_data[bioset], self.nextbio_data[bioset_tofix])
                        self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
  
    ''' fixLastSnpEthnicity '''
    def fixLastSnpEthnicity(self, data):
        logging.debug(' Function:  fixLastSnpEthnicity- %s' %self.bioset_id)
        for snp in data['snps']:
            snp_data = copy.deepcopy(data['snps'][snp])
            if data['hapmap_pop'] != data['snps'][snp]['snp_population']:
                snp_data['snp_population'] = data['hapmap_pop'] 
                self.util.warnMe('warning',' CHANGED  snp_population (fixLastSnpEthnicity) - snp_population:%s, hapmap_pop:%s' % ( data['snps'][snp]['snp_population'], data['hapmap_pop'] ))
            data['snps'][snp] = copy.deepcopy(snp_data)
 

    ''' dataForHapMapMinorAllele 
    '''
    def dataForHapMapMinorAllele(self, data):
        logging.debug(' Function:  dataForHapMapMinorAllele - %s' %self.bioset_id)
        for snp in data['snps']:
            snp_data = copy.deepcopy(data['snps'][snp])
            ### find minor allele in hapmap
            if snp_data['allele_1_frequency'] >  snp_data['allele_2_frequency']:
                self.snp_issue[snp] = snp_data
                self.snp_issue[snp]['issue'] = ['_odds_ratio','95%ci']
                snp_data = self.fixDataFromAllele2Issue(snp_data) 
            data['snps'][snp] = copy.deepcopy(snp_data)


    ''' fixDataToFixedEffect 
    '''
    def fixDataToFixedEffect(self, data):
        logging.debug(' Function:  fixDataToFixedEffect - %s' %self.bioset_id)
        for snp in data['snps']:
            snp_data = copy.deepcopy(data['snps'][snp])
            snp_data['p-value'] = data['snps'][snp]['p-value_(fixed_effect)']
            snp_data['allele_1_odds_ratio'] = data['snps'][snp]['odds_ratio_(fixed_effect)']
            snp_data['allele_1_or95%ci'] = data['snps'][snp]['or_95%ci_(fixed_effect)']
            snp_data['p-value(original)'] = data['snps'][snp]['p-value']
            snp_data['allele_1_odds_ratio(original)'] = data['snps'][snp]['allele_1_odds_ratio']
            
            ### add newkey to list of column
            for new_key in snp_data:
                if new_key.upper() not in self.column_key[self.bioset_id]:  
                    self.column_key[self.bioset_id].append(new_key.upper())
            data['snps'][snp] = copy.deepcopy(snp_data)

    ''' fixPvalue 
    '''
    def fixPvalue(self, data):
        logging.debug(' Function:  fixPvalue - %s' %self.bioset_id)
        for snp in data['snps']:
            snp_data = copy.deepcopy(data['snps'][snp])
            snp_data['p-value'] = data['snps'][snp]['p-value_(cmh)']
            data['snps'][snp] = copy.deepcopy(snp_data)
    
    def fixDataForCurration(self,data):   
        logging.debug(' Function:  fixDataForCurration' )
        self.nextbio_data = {data['bioset_id'] : data}  
        self.bioset_info = {data['bioset_id'] :data['info'] }
        self.column_key = {data['bioset_id'] : [element.upper() for element in data['columns']] }      
        self.fixData() 
        self.cooked_data[0]['columns'] =  self.column_key[data['bioset_id']]
                             
    ''' fixData '''
    def fixData(self):
        logging.debug(' Function:  fixData' )

        self.snp_issue = {} ### record snps with specific issue to recall for other biosets
        self.cooked_data=[] ## record fixed daa
        self.tofix ={} ### record bioset reference and bioset they fix
        flag = False
        for bioset in self.nextbio_data:
            self.bioset_id = bioset
            if 'snp_issue' in self.nextbio_data[bioset]:
                self.snp_issue = self.nextbio_data[bioset]['snp_issue']
            if 'snp_issue' in self.bioset_info[bioset]:
                self.snp_issue = self.bioset_info[bioset]['snp_issue']
            if self.bioset_info[bioset]['issue'].lower() == "minor_allele != allele_1 -> allele_2_<data>":
                print 'A', bioset
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
                flag = True
            elif 'allele_2_<data> -> bioset_id ' in self.bioset_info[bioset]['issue'].lower() :
                print 'B', bioset
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
            elif 'missing genotype -> bioset_id ' in self.bioset_info[bioset]['issue'].lower() :
                print 'C', bioset
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                ref_bioset = self.bioset_info[bioset]['issue'].lower().strip().split('missing genotype -> bioset_id ')[1]
                if ref_bioset not in self.tofix:
                    print 'D', bioset
                    self.tofix[ref_bioset] = {bioset:'missing genotype -> bioset_id '}
                else:
                    print 'F', bioset
                    self.tofix[ref_bioset][bioset] ='missing genotype -> bioset_id '
            elif 'global' in  self.bioset_info[bioset]['issue'].lower() :
                print 'I', bioset 
                self.fixLastSnpEthnicity(self.nextbio_data[bioset])
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
            elif '-> or for allele_2' in self.bioset_info[bioset]['issue'].lower() :
                snp = self.bioset_info[bioset]['issue'].lower().split(' ')[0]
                print 'K snp allele 2 data',bioset, snp, self.bioset_info[bioset]['issue'].lower() 
                self.snp_issue[snp] = self.nextbio_data[bioset]['snps'][snp]
                self.snp_issue[snp]['issue'] = self.bioset_info[bioset]['issue'].lower().split('|')[1].split(',')
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
                self.snp_issue = {} 
            elif 'data given for minor allele in hapmap' == self.bioset_info[bioset]['issue'].lower():
                print 'L data for minor hapmap',bioset, self.bioset_info[bioset]['issue'].lower() 
                self.snp_issue = {} # reinit 
                self.dataForHapMapMinorAllele(self.nextbio_data[bioset])
                self.snp_issue = {}
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
            elif 'fixed_effect' in self.bioset_info[bioset]['issue'].lower() :
                print 'M fixed_effect',bioset, self.bioset_info[bioset]['issue'].lower() 
                self.fixDataToFixedEffect(self.nextbio_data[bioset])
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset])) 
            elif 'p-value_(cmh) -> p-value' in self.bioset_info[bioset]['issue'].lower() :
                print 'N pvalue issue',bioset, self.bioset_info[bioset]['issue'].lower() 
                self.fixPvalue(self.nextbio_data[bioset])
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset])) 
            elif '->' in self.bioset_info[bioset]['issue'].lower() :
                print 'H NEW RULE',bioset, self.bioset_info[bioset]['issue'].lower() 
                print self.bioset_info[bioset]
            elif 'allele 1 based on hapmap and allele 1 from the publication are not consistent!' in self.bioset_info[bioset]['issue'].lower() :
                print 'O', bioset, self.bioset_info[bioset]['issue'].lower() 
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
            else:
                print 'G', bioset,  self.bioset_info[bioset]['issue'].lower()
                if flag == False:
                    self.snp_issue = {} # reinit 
                self.checkIfDataAllele1OrMinor(self.nextbio_data[bioset])
                self.cooked_data.append(self.util.keyToUpperCase(self.nextbio_data[bioset]))
            


    ### write in Nextbio style batch info ###
    def writeNextBioStyle(self, filename,data):
        logging.debug(' Function:  writeNextBioStyle %s' % filename)
        wbk = xlwt.Workbook()
        data_sheet = wbk.add_sheet("data")
        info_sheet = wbk.add_sheet("info")
   
        info_data_row = 0  
        data_row = 0
        data_row = self.excel_util.writeExcelRow(data_sheet,data_row,"=======================================================================================================================================================================================================================================" ,None)
        for d in data:
            bioset_id = d["BIOSET_ID"]
#             print bioset_id
            ### fill hand curated info for this bioset
            info_data_row, info_sheet = self.excel_util.writeInfoSheet(info_sheet,info_data_row, self.snp_table_header,self.bioset_info[bioset_id])
            
            ### fill nextbio data for this bioset
            data_row, data_sheet= self.excel_util.writeNextBioBiosetInfo(data_sheet,data_row,d)
  
            ### write SNP table snp_table_header
            data_row, data_sheet , self.column_key[bioset_id] = self.excel_util.writeExcelSnpTableHeader(data_sheet, data_row, self.column_key[bioset_id])
            ### write SNP data
            data_row, data_sheet  = self.excel_util.writeExcelSnpTableData(data_sheet, data_row, self.column_key[bioset_id], d['SNPS'])

        wbk.save("%s/%s.xls" % (config.OUTPUTPATH,filename))


    ### main for class ###
    def main(self):
        logging.debug(' Function:  FixHandCurratedData main %s' % self.disease_name)
        ''' get data for this disease '''
        ### read excel data
        self.nextbio_data, self.bioset_info, self.snp_table_header, self.column_key= self.excel_util.readExcelData("%s/%s_CB_QC.xls" % (config.DATAPATH,disease_name))

        ### fix data 
        self.fixData()
        ### record json for documentation
        self.util.writeOutput('cooked_%s' % self.disease_name,self.cooked_data)
        ### otput in xls
        self.writeNextBioStyle('cooked_%s' % self.disease_name,self.cooked_data )
        
        

""" MAIN function """    
if __name__ == "__main__":
    ## age_related_macular_degeneration, gallstone, venous_thrombosis, gout, atrial_fibrillation, allergic_rhinitis, gout,kidney_stone,alcohol_dependence,melanoma,rheumatoid_arthritis, age_related_macular_degeneration
    disease_name = "allergic_rhinitis"
    logging.basicConfig(filename='%s/%s%s_%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'fixData',disease_name), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean =  FixHandCurratedData(disease_name)
    clean.main()
    
    print'DONE with DMmain'
    
""" END OF main """ 