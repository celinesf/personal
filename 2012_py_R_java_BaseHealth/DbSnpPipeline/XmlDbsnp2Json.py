#!/usr/bin/env python

"""
     Convert xml dbsnp data into json
     06/14/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


import logging, json, copy
import DbSnpConfig as config
from DbSnpUtils import DbSnpUtils
from SimpleDbSnpInfo import SimpleDbsnpInfo
from xmlloader import *


class XmlDbsnp2Json():
    def __init__(self, filename):
        self.util = DbSnpUtils()
        self.filename =filename
        self.simpleStructure = SimpleDbsnpInfo()
      

    ''' xmlToJson 
    BEFORE STARTING
        remove : from xmi key
        remove &amp; by  \u0026 signs
        &gt; to  \u003E
    '''
    def xmlToJson(self):
        logging.debug(' Function: xmlToJson') 
        self.json_data = {}
        
        input_data = file('%s/%s' %(config.DATAPATH,self.filename), 'r') 
#         input_data = file(config.INPUTFILE, 'r')   # dbsnp_list_ A document containing XML
        line = input_data.readline().replace('&amp;','\u0026').replace('&gt;','\u003E')
        info = None
        ns = 0
        while len(line)>0:
            line = input_data.readline()
            line = line.replace('    <', '<').replace('&amp;','\u0026').replace('&gt;','\u003E')
            if '<Rs rsId=' in line:
                info = line
            elif info is not None:
                info = "%s%s" % (info,line)
                ### focus on 1 snp info at a time
                if '</Rs>' in line:
                    out = open('%s/xml_temp_%s' % (config.OUTPUTPATH ,self.filename), 'w')
                    out.write(info)
                    out.close()
                    in_data = open('%s/xml_temp_%s' % (config.OUTPUTPATH ,self.filename), 'r')
                    xml_data =  StreamXMLLoader(in_data, 0)    
                    data = xml_data.expectXML()
                    dbsnp = data['Rs']['__attrs__']['rsId']  
                    logging.info('SNP : %s' % dbsnp)     
                    ### add simplify+ gather data functions
#                     self.util.warnMe('info', 'start for snp %s' %dbsnp)
                    self.json_data['rs%s' % dbsnp] = self.simpleStructure.simplifyStructure(data['Rs'])
                    
                    ### add to mongo
#                     print dbsnp
#                     nextbio_db['dbsnp'].remove({'dbsnp':'rs%s' % dbsnp}) 
#                     nextbio_db['dbsnp'].insert(copy.deepcopy(self.json_data['rs%s' % dbsnp]))                
                    in_data.close()
                    ns+=1
#                     if ns ==10: break
                    if ns %100 == 0:
                        print ns, 'rs%s' %dbsnp
        input_data.close()
#         print  json.dumps(self.json_data)    
        self.util.writeOutput('%s_xml2json' %self.filename, self.json_data)


""" MAIN function """    
if __name__ == "__main__":

    filename = '131029175606'
    ##################
    logging.basicConfig(filename='%s/%s%s' % (config.OUTPUTPATH,config.LOG_FILENAME ,'%s_XmlDbsnp2Json' % filename), filemode='w',
                        level=logging.INFO,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__' )

    clean = XmlDbsnp2Json(filename)
    ######## run to convert xml to json ####
#     clean.removeChar()
    clean.xmlToJson()
    ########
    
    print'DONE with ', filename
    
""" END OF main """ 