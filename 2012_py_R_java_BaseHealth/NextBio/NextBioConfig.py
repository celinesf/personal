"""
    Config file for NextBio parser/ clean up scripts
   05/06/13 - 1.0 
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

################ connect to mongo  ###############
MONGOHOST= "54.244.22.12" ### tiger
MONGO = "mongodb://adm1n:adm1npw0rd@%s:27017/admin" % MONGOHOST 
### PATH TO CHANGE
PATH = "/Users/celine/Box Documents/Celine/ScriptInOut/NextBio/08-20-13"
LOG_FILENAME = "logfile"
################

################ step 0 cleanup_encoding ################
# DATAPATH_ORI = "%s/0_original_data" % PATH
# OUTPUTPATH = "%s/1_cleaned_2_curation" % PATH
# CODEMAP = "%s/maps/code_auto_correction.json" % PATH
################


################ step 1 GetData4HandCurration ################
# DATAPATH = "%s/1_cleaned_2_curation"% PATH
# OUTPUTPATH = "%s/2_curration_2_dbsnp" % PATH
# ### DISEASEMAP - change as needed
# # DISEASEMAP = "%s/maps/3.1_disease_lexicon_test.json" % PATH
# # DISEASEMAP = "%s/maps/3.1_disease_lexicon_v2.json" % PATH
# DISEASEMAP = "%s/maps/3.2_disease_lexicon_v4.json" % PATH
# ###
# KEYMAP = "%s/maps/key_map.json" % PATH
# ETHNICITYMAP = "%s/maps/ethnicity_map.json" % PATH
# DATATYPEMAP = "%s/maps/data_type.json" % PATH
# MODELMAP = "%s/maps/model_map.json" % PATH
# INFOHEADER = "%s/maps/bioset_info_header.json" % PATH
# METAMAP = "%s/maps/bioset_meta_map.json" % PATH
################

################ step 2 GetSnpList ################
# DATAPATH = "%s/2_curration_2_dbsnp" % PATH
# OUTPUTPATH = "%s/3_get_snp_2_dbsnp" % PATH
################


################ step 5 GetFormatedData ################
DATAPATH = "%s/5_currated_2_db"% PATH
OUTPUTPATH = "%s/6_format_output" % PATH
SNPINFOMAP = "%s/maps/snp_data_map.json" % PATH
FUNCTIONOMAP = "%s/maps/function_map.json" % PATH
INFOHEADER = "%s/maps/bioset_info_header.json" % PATH
METAMAP = "%s/maps/bioset_meta_map.json" % PATH
CURRATIONMAP ="%s/maps/hand_curration_map.json" % PATH 


################ FixGeneticsRepStrands ################
# OUTPUTPATH = "%s" % PATH
