"""
    Config file for Algogene pipeline for nextbio data
   06/20/13 - 1.0 
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

import pymongo

LOG_FILENAME = "logfile"
OUTPUTPATH = "/Users/celine/Box Documents/Celine/ScriptInOut/AlgoGene"
LD_LIMIT = 1e5 # default = 1e5
POP_LIMIT = 1000 # default = 1000
PVALUE = 0.05 # default = 0.05
SIGN = 3 # default = 3

########### connect to mongo  ###############
MONGOHOST= "54.245.239.40" ### tiger
MONGO = "mongodb://adm1n:adm1npw0rd@%s:27017/admin" % MONGOHOST 
CW = pymongo.Connection(MONGO, 27017)
GENOPHEN_DB = CW["genophen30"]
NEXTBIO_DB = CW["nextbio"]