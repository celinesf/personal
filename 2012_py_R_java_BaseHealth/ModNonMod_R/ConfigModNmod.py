"""
    Config for GetModNonModRiskFactors
   06/04/13: version 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 


LOG_FILENAME = "logfile"
OUTPUTPATH  = "/Users/celine/Box Documents/Genetic Scripts/Comprehensive/9-12-13_LC"
ROUTPUTPATH  = "/Users/celine/Box\ Documents/Genetic\ Scripts/Comprehensive/9-12-13_LC"
DATAPATH = "%s" % OUTPUTPATH

########### connect to mongo  ###############
MONGOHOST= "54.244.22.12" ### tiger
MONGO = "mongodb://adm1n:adm1npw0rd@%s:27017/admin" % MONGOHOST 
