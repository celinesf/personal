#!/usr/bin/env python

"""
     Utility functions  AlgoGene 
     06/24/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging,  json, sys
import AlgoGeneConfig as config

class AlgoGeneUtil():
    def __init__(self):
        ''' '''

######################### OUT FUNCTIONS ##############

    ### fecth_input_data ###
    # - open an input file
    # - convert the content from json to python dictionary
    # called by main
    ### 
    def fecth_input_data(self,filename, comment):
        logging.debug(' Function: fecth_input_data: filename: %s, comment %s' % (filename, comment))
        try:
            logging.info(' ***** Opening %s file: %s' % (comment, filename))
            inputfile = open(filename,'r')
        except:
            logging.critical(' COULD NOT OPEN VARIANT DATA FILE (main): %s' % filename)
            
        ### convert json
        try:
            inputdata = json.loads(inputfile.read())
            inputfile.close()
        except:
            logging.critical(' COULD NOT CONVERT %s FROM JSON (main): %s' % (comment.upper(),filename))
            sys.exit(' COULD NOT CONVERT %s FROM JSON (main): %s' % (comment.upper(),filename))
        return inputdata
    ### END OF fecth_input_data  

#     ''' writeJsonFile '''
#     def writeJsonFile(self, to_change, document):
#         logging.debug(' Function:  writeJsonFile - %s' % to_change)
#         json_file = '%s/%s_%s.json' % (config.OUTPUTPATH,self.disease_name,to_change )
#         output = open(json_file,'w')
#         output.write ( json.dumps(document, indent=4))
#         output.write('\n')
#         output.close()

    def writeOutput(self, filename,data):
        logging.debug(' Function:  writeOutput %s:' % filename)
        f = open("%s/%s.json" % (config.OUTPUTPATH,filename), "w")
        f.write(json.dumps(data, indent=4))
        f.close()
        
    ''' write unique key of list '''
    def writeList(self, filename,data):
        logging.debug(' Function:  writeList %s:' % filename)
        total = 0
        f = open("%s%s.json" % (config.OUTPUTPATH,filename), "w")
        for d in data:
            f.write("%s\t%s\n" % (d, data[d]["count"])) 
            total +=  data[d]["count"]
        f.write("TOTAL\t%s\n" % (total))

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



