#!/usr/bin/env python
"""
   Prepare input/ database like files for genetic algorithm
   10/29/13 - 0.0.1 : creation
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Viragene Inc."
__maintainer__ = "Celine Becquet"
__email__ = "celine@vitagene.com"
__status__ = "dev" 

"""
Prepare input/ database like files for genetic algorithm
"""

PATH = "/Users/Celine/vitagene/script_in_out/pharmGKB"
from SharedFunctions import SharedFunctions
import logging, copy

class TraitSuppSnpInput():
    def __init__(self):
        self.util = SharedFunctions(PATH)
        self.interactions_map = self.util.get_map('interactions_map') 
        self.traits_map = self.util.get_map('traits_map') 
        self.supplements_map = self.util.get_map('supplements_map') 
        self.output_data= []
        self.new_header = ['trait', 'supplement', 'interaction', 'gender', 'snp', 'genotype', 'score', 'effect', 'pmid', 'gene', 'note']
    
    def check_map(self,check_map, word):
        logging.debug(' Function: check_map word: %s' % word )
        
        key = word.lower().replace(' ',"_")
        if key in check_map and check_map[key] == word:
            logging.info(' found key: %s (check_map)' % key )
            word = key
        elif key in check_map and check_map[key] != word:
            logging.warning(' found key: %s but different name %s %s (check_map)' % (key, check_map[key] , word))
            print(' found key: %s but different name %s %s (check_map)' % (key, check_map[key] , word))
        elif '/' in word:
            words = word.split('/')
            new_word=[]
            for w in words:
                w = self.check_map(check_map, w)
                new_word.append(w)
            word = copy.deepcopy(new_word)
        elif key != 'null' and key != '':
            logging.warning(' NEW KEY: %s (check_map)' % (key))
            check_map[key] = word
            word = key
            print(' NEW KEY: %s (check_map)' % (key))
        return word

    def extract_genotype_data(self,key,value):
        logging.debug(' Function: extract_genotype_data key: %s, value: %s' % (key,value) )
        genotype = key.split('_')[1]
        if genotype not in self.genotype:
            self.genotype[genotype]={}
        self.genotype[genotype][key.split('_')[0]] = value

    def format_data(self,words):
        logging.debug(' Function: format_data line: %s' % words )
        self.data = {}
        self.genotype = {}
        for i in range(0,len(self.header)):
            if self.header[i] == 'trait' :
                words[i] = self.check_map(self.traits_map, words[i])
                self.data[self.header[i]] = words[i]
            elif self.header[i] == 'supplement' :
                words[i] = self.check_map(self.supplements_map, words[i])
                self.data[self.header[i]] = words[i]
            elif self.header[i] == 'interaction' :
                words[i] = self.check_map(self.interactions_map, words[i])
                self.data[self.header[i]] = words[i]
            elif 'a1' in self.header[i] or 'a2' in self.header[i]:
                self.extract_genotype_data(self.header[i],words[i])
            else: 
                self.data[self.header[i]] = words[i]
        for geno in self.genotype:
            self.genotype[geno].update(self.data)
            self.output_data.append(copy.deepcopy(self.genotype[geno]))


    def extract_data(self):
        logging.debug(' Function: extract_data' )
        ori_data = open("%s/map/%s" %(PATH, "SNP_supp_disease.txt"))
        lines = ori_data.readlines()
        nl = 0
        for line in lines:
            if nl == 0:
                self.header = self.util.get_header(line, None)
            elif nl>1:
                self.format_data(line.strip().split('\t'))
            nl += 1
        self.util.output_file('traits_map',self.traits_map)
        self.util.write_map('traits_map',self.traits_map)
        self.util.output_file('supplements_map',self.supplements_map)
        self.util.write_map('supplements_map',self.supplements_map)
        
        self.util.output_file('interactions_map',self.interactions_map)
        self.util.write_map('interactions_map',self.interactions_map)
        self.util.output_file('trait_supp_snp',self.output_data)
        self.util.write_output('trait_supp_snp',self.new_header,self.output_data)

### main function
if __name__ == "__main__":
    logging.basicConfig(filename='%s/%s' % (PATH,'log_TraitSuppSnpInput'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__ input' )
    d = TraitSuppSnpInput()
    d.extract_data()
    
    print(' DONE -- Function: __main__' )
    logging.debug(' DONE -- Function: __main__' )
""" END OF main """
