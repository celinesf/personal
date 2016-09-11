#!/usr/bin/env python
"""
   Check strands of data from PharmGKB as defined by dbsnp
   + complete the data with positions minor and major
   10/29/13 - 0.0.1 : creation
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Viragene Inc."
__maintainer__ = "Celine Becquet"
__email__ = "celine@vitagene.com"
__status__ = "dev" 


PATH = "/Users/Celine/vitagene/script_in_out/pharmGKB/check_strands"

import logging, re, copy
from SharedFunctions import SharedFunctions

class CheckStrand():
    def __init__(self):
        self.header = ["chromosome","position","alleles","reference","minor_allele","ancestral","functions","effect_allele","merged"]
        self.util = SharedFunctions(PATH)
        self.dbsnp_map = self.util.get_map('dbsnp_map')
        self.trait_supp_snp = self.util.get_map('trait_supp_snp')    
        self.function_map = self.util.get_map('function_map')
        self.new_header = self.util.get_map('header_map')
        self.trait_supp_dbsnp = []
    
    ''' setDbSnpFunction  + assign chr and position in proper format'''
    def setDbSnpFunction(self, snp_doc):
        logging.debug(' Function:  setDbSnpFunction - %s' % self.dbsnp) 
        
        if snp_doc['functions'] is not None:
            self.new_data['functions'] = None
            for doc in snp_doc['functions']:
                ### check gene name
                if doc['symbol'] is not None:
                    symbol = doc['symbol'].split('-AS1')[0]
                    if symbol.lower()  not in self.new_data['gene'].lower():
                        logging.info(' ADDING DBSNP GENE (setDbSnpFunction) - snp: %s, gene:%s, dbGene: %s' %(self.dbsnp,self.new_data['gene'], doc['symbol'] ))
                        self.new_data['gene'] = self.new_data['gene'] + '/' + symbol
                    ### add function
                if doc['fxnClass'] != 'reference':
                    if self.new_data['functions'] is None:
                        self.new_data['functions'] = doc['fxnClass']
                    else:
                        new_index = self.function_map.index(doc['fxnClass'])
                        old_index = self.function_map.index(self.new_data['functions'])
                        if new_index < old_index:
                            self.new_data['functions'] = doc['fxnClass']
                            self.new_data['effect_allele'] = doc['allele']
    
    def change_strand_text(self, text):
        logging.debug(' Function:  change_strand_text - text %s' % text) 
        old_geno = re.findall('[A,C,T,G]{2}' , text)

        for g in old_geno:
            new_geno = ''.join(self.util.reverseAlleles(g))
            text=text.replace(g,new_geno)
        return text
    
    def check_strand(self):
        logging.debug(' Function:  check_strand - snp %s' % self.dbsnp ) 
        if self.new_data['genotype'][0] not in self.new_data['alleles']:
            self.new_data['genotype'] = ''.join(self.util.reverseAlleles(self.new_data['genotype']))
            self.new_data['note'] = self.change_strand_text(self.new_data['note'])
            self.new_data['effect'] = self.change_strand_text(self.new_data['effect'])            
    
    def complete_snp_data(self):
        logging.debug(' Function:  complete_snp_data -') 
        for d in self.trait_supp_snp:
            self.new_data = copy.deepcopy(d)
            if d['snp'] in self.dbsnp_map:
                self.dbsnp = d['snp']
                for h in range(0,len(self.header)):
                    if self.header[h] in self.dbsnp_map[d['snp']]:
                        if self.header[h] == "functions":
                            # find the worth function and associated allele"
                            self.setDbSnpFunction(self.dbsnp_map[d['snp']])
                        else:
                            self.new_data[self.header[h]] = self.dbsnp_map[d['snp']][self.header[h]]
                    else:
                        if self.header[h] == "gene" or self.header[h] == "effect_allele":
                            # check symbol concatenation with gene
                            ''' '''
                        else:
                            self.util.warnMe('warning', " key not found %s (complete_snp_data)" %(self.header[h]))      
            else:
                self.util.warnMe('critical', " SNP not found %s (complete_snp_data)" %(d['snp']))
            self.check_strand()
            self.trait_supp_dbsnp.append(copy.deepcopy(self.new_data))

      
        self.util.output_file('trait_supp_dbsnp',self.trait_supp_dbsnp)
        self.util.write_output('trait_supp_dbsnp',self.new_header,self.trait_supp_dbsnp)
        
### main function
if __name__ == "__main__":
    logging.basicConfig(filename='%s/%s' % (PATH,'log'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__ input' )
    d = CheckStrand()
    d.complete_snp_data()

    print(' DONE -- Function: __main__' )
    logging.debug(' DONE -- Function: __main__' )
""" END OF main """
