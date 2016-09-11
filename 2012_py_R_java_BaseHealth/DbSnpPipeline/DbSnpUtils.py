#!/usr/bin/env python

"""
     Utility functions for DbSnpPipeline
     06/16/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging, json
import DbSnpConfig as config

class DbSnpUtils():
    def __init__(self):
        ''' nothing '''
    
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
      
    ''' get map for typo corrections '''
    def getTermMap(self,fname):
        logging.debug(' Function:  getTermMap' )
        f = open(fname,'r')
        typos = f.read()  
        f.close()
        return json.loads(typos)
        
            
    def writeOutput(self, filename,data):
        logging.debug(' Function:  writeOutput %s:' % filename)
        f = open("%s/%s.json" % (config.OUTPUTPATH,filename), "w")
        f.write(json.dumps(data, indent=4))
        f.close()
        
    ''' write unique key of list '''
    def writeList(self, filename,data):
        logging.debug(' Function:  writeList %s:' % filename)
        total = 0
        f = open("%s/%s.json" % (config.OUTPUTPATH,filename), "w")
        for d in data:
            f.write("%s\t%s\n" % (d, data[d]["count"])) 
            total +=  data[d]["count"]
        f.write("TOTAL\t%s\n" % (total))
    
  
 