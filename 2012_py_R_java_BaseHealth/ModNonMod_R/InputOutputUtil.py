#!/usr/bin/env python

"""
    Utility functions used is every module
    06/1/13: version 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 
"""
    1) writeJsonFile
    2) warnMe
"""

import logging,json

class InputOutputUtil():
    def __init__(self,path , disease_name):
        self.path = path
        self.disease_name = disease_name
    
##################### utility functions  ################### 
    ''' writeJsonFile '''
    def writeJsonFile(self, to_change, document):
        logging.debug(' Function:  writeJsonFile - %s' % to_change)
        json_file = '%s/%s_%s.json' % (self.path,self.disease_name,to_change )
        output = open(json_file,'w')
        output.write ( json.dumps(document, indent=4))
        output.write('\n')
        output.close()
        
    ''' warnMe '''   
    def warnMe(self, flag, comment):
        if flag == 'info':
            logging.info(comment)
            print comment
        elif flag == "warning":
            logging.warning(comment)
            print "WARNING -%s" % comment
    
#####################    