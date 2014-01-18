#!/usr/bin/env python
"""
   Prepare input/ database like files for Matt's Algorithm
   10/22/13 - 0.0.1 : creation
   10/29/13 - 0.0.2 : updated to accoutn for shared functions
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Viragene Inc."
__maintainer__ = "Celine Becquet"
__email__ = "celine@vitagene.com"
__status__ = "dev" 

"""
Functions to Prepare input/ database like files for Matt's Algorithm
"""
PATH = "/Users/Celine/vitagene/script_in_out/supp2disease"

import logging, re, copy
from SharedFunctions import SharedFunctions

class TraitSuppInput():
    def __init__(self):
        self.CATEGORY_MAP = {'A':'acute', 'B' : 'beauty', 'C' : 'chronic', 'P':'performance','W' : 'wellness'}
        self.traits_map ={}
        self.supplements_map ={}
        self.output_data = []
        self.header =None
        self.util = SharedFunctions(PATH)
        self.score_map = self.util.get_map('score_map')
        self.sub_trait_map = self.util.get_map('sub_trait_map')    
    
    def recover_trait_info(self,words):
        logging.debug(' Function: recover_trait_info words: %s' % words )
        self.new_trait = {}
        for i in range(0,3):
            if self.header[i] == 'category':
                categories = words[i].strip().split('/')
                words[i] = self.CATEGORY_MAP[categories[0]]
                for cat in range(1,len(categories)):
                    if categories[cat] in self.CATEGORY_MAP:
                        words[i] = words[i] + "/"+ self.CATEGORY_MAP[categories[cat]]
                    else:
                        print "WARNING I DONT KNOW WHAT THIS IS (recover_trait_info) - %s %s" % (cat,categories[cat])
                        logging.debug("WARNING I DONT KNOW WHAT THIS IS (recover_trait_info)- %s %s" % (cat, categories[cat]))
                
            self.new_trait[self.header[i]]= words[i].strip().lower().replace(' ','_')
            ### record disease/trait name
            if self.header[i] == 'trait' and words[i].strip().lower().replace(' ','_') not in self.traits_map:
                self.trait = words[i].strip().lower().replace(' ','_')
                self.new_trait[self.header[i]] = self.trait
                self.traits_map[self.trait] =  words[i].strip()

    def recover_score(self,words):
        logging.debug(' Function: recover_score - words: %s' %words)
        for i in range(3,len(words)):
            self.new_data = copy.deepcopy(self.new_trait)
            w = re.split('[(/:)]', words[i].strip())
            supp = w[0].strip()
            self.supp = supp.lower().replace(' ','_') 
            ### add supp to map
            if  self.supp not in self.supplements_map:
                self.supplements_map[ self.supp] = supp
            self.new_data["supplement"] =  self.supp
            for i in range(1,len(w)):
                if w[i] in self.score_map:
                    self.new_data[self.score_map[w[i]]['site']] = self.score_map[w[i]]['score']
                if w[i] not in self.score_map and w[i] != "":
                    sub_traits = w[i].strip().split(', ')
                    self.new_data['vitaganic_sub_trait'] = sub_traits[0].strip().replace(' ','_')
                    for sub in range(1,len(sub_traits)):
                        self.new_data['vitaganic_sub_trait'] = self.new_data['vitaganic_sub_trait']+ '/' + sub_traits[sub].strip().replace(' ','_')
                        if sub_traits[sub].strip().replace(' ','_') not in self.sub_trait_map:
                            print 'HERE', sub_traits[sub], words                
            self.output_data.append(copy.deepcopy(self.new_data))

    
    def extract_data(self):
        logging.debug(' Function: extract_data' )
        ori_data = open("%s/map/%s" %(PATH, "disease_supp_scores"))
        lines = ori_data.readlines()
        nl = 0
        for line in lines:
            if nl == 1:
                self.header = self.util.get_header(line,['dose','mayo_clinic','webmd','vitaganic','vitaganic_sub_trait'])
            elif nl>1:
                self.recover_trait_info(line.strip().split('\t'))
                self.recover_score(line.strip().split('\t'))
            nl += 1

        self.util.output_file('traits_map',self.traits_map)
        self.util.output_file('supplements_map',self.supplements_map)
        self.util.output_file('trait_supp_score',self.output_data)
        self.util.write_map('traits_map',self.traits_map)
        self.util.write_map('supplements_map',self.supplements_map)
        self.util.write_output('trait_supp_score',self.header,self.output_data)

        
### main function
if __name__ == "__main__":
    logging.basicConfig(filename='%s/%s' % (PATH,'log'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__ input' )
    d = TraitSuppInput()
    d.extract_data()

    print(' DONE -- Function: __main__' )
    logging.debug(' DONE -- Function: __main__' )
""" END OF main """
