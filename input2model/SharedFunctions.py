#!/usr/bin/env python
"""
   Functions used by multiple classes
   10/29/13 - 0.0.1 : 
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Viragene Inc."
__maintainer__ = "Celine Becquet"
__email__ = "celine@vitagene.com"
__status__ = "dev" 

"""
Functions used by multiple classes
"""


import logging, json, copy

class SharedFunctions():
    def __init__(self, path):
        ''' '''
        self.path = path
        
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
    
    def get_map(self, file_name):
        logging.debug(' Function: get_map - filename: %s' % (file_name) )
        score_map = open("%s/map/%s" %(self.path, file_name))
        return json.loads(score_map.read())
    
    def get_header(self,line, add):
        logging.debug(' Function: get_header line: %s' % line )
        header = line.lower().replace(' ','_').strip().split('\t')
        if add is not None:
            header.extend(add)
        return header
    
    def output_file(self,name,map_data):
        logging.debug(' Function: output_file - name :%s' %name )
        output = open('%s/%s' % (self.path,name),'w')
        output.write(json.dumps(map_data,indent=4))
        output.close()
    
    def write_map(self,name,data): 
        logging.debug(' Function: write_map - name :%s' %name )
        output = open('%s/%s.txt' % (self.path,name),'w')
        output.write('key\tvalue\n')
        for key in data:
            output.write('%s\t%s\n' % (key, data[key]))
        output.close()
    
    def write_output(self,name,header,data): 
        logging.debug(' Function: write_output - name :%s' %name )
        output = open('%s/%s.txt' % (self.path,name),'w')
        
        line = "#"
        for col in header:
            line += '\t' +col
        output.write(line +'\n')

        ndoc = 0
        for doc in data:
            line = str(ndoc)
            for col in header:
                if col in doc:
                    if type(doc[col]) is list:
                        doc[col]= '/'.join(doc[col])
                    line += '\t' + str(doc[col])
                else:
                    line += '\t' + "NULL"
            output.write(line +'\n')
            ndoc +=1     
        output.close()
    
    ''' warnMe '''   
    def warnMe(self, flag, comment):
        if flag == 'info':
            logging.info(comment)
            print comment
        elif flag== 'warning':
            logging.warning(comment)
            print "%s -%s" % (flag.upper(),comment)
        elif flag== 'critical':
            logging.critical(comment)
            print "%s -%s" % (flag.upper(),comment)
        elif flag== 'error':
            logging.error(comment)
            print "%s -%s" % (flag.upper(),comment)
        else:
            logging.info(comment)
            print "%s -%s" % (flag.upper(),comment) 

        ''' reverse
    '''
    def reverseAlleles(self, alleles):
        logging.debug(' Function: reverseAlleles - %s' % alleles)  
        reversed_allele = []
        for index in range(0,len(alleles)):
            allele= alleles[index ]
            if allele == 'G':
                reversed_allele.append('C')
            elif allele == 'C':
                reversed_allele.append('G')
            elif allele == 'A':
                reversed_allele.append('T')
            elif allele == 'T':
                reversed_allele.append('A')
            elif allele == '-':
                reversed_allele.append('-')
            else :
                self.util.warnMe('warning', ' CANT REVERSE ALLELE %s (reverseAlleles) - alleles:%s' % (allele,alleles))  
        if len(reversed_allele) == 1:
            reversed_allele = reversed_allele[0]
        return reversed_allele

