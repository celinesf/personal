"""
    Config file for for DbSnpPipeline
     06/16/12- 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

PATH = "/Users/Celine/vitagene/script_in_out/pharmGKB/dbsnp"
LOG_FILENAME = "logfile"
OUTPUTPATH = '%s' % PATH
DATAPATH = '%s' % PATH
INPUTFILE = '%s/maps/dbsnp_list_xml' % PATH
XMLMAP = '%s/maps/dbsnp_xml_map.json' % PATH


########### connect to mongo  ###############
# MONGOHOST= "54.244.22.12" ### tiger
# MONGO = "mongodb://adm1n:adm1npw0rd@%s:27017/admin" % MONGOHOST 
