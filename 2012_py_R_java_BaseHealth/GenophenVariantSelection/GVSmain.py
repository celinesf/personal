#!/usr/bin/python

'''
Created on Jul 13, 2012
@author: celine
'''

import json
import sys
import logging

### log file name
#from datetime import * 
LOG_FILENAME = "logfile" #_%s" % str(datetime.now())

### variant scoring class functions
import VariantScoring as VS
VS_class = VS.VariantScoring()

#from VariantScoring import VariantScoring
### LD variant selection class
import LDVariantSelection as LD
LD_class = LD.LDVariantSelection()

### shared functions and variables
import SharedFunctions as SF
STAT = SF.STATISTICS
SF_class = SF.SharedFunctions()

### summaries_disease ###
# clean dictionaries for output purposes
### 
def summaries_disease():
    logging.debug(' Function: summaries_disease')
    for et in STAT['variant_scoring']['different_ethnicity']:
        for key in STAT['variant_scoring'][et]:
            STAT['variant_scoring'][key] += STAT['variant_scoring'][et][key]
        for key in STAT['LD_selection'][et]:
            STAT['LD_selection'][key] += STAT['LD_selection'][et][key]
### END OF summaries_disease

### initialize_statistics ###
# Initialize the summary statistics dictionary
### 
def initialize_statistics(stat, num):
    logging.debug(' Function: initialize_statistics')
    if num == 0:
        stat['different_ethnicity'] =[] # IN
        stat['total_ethnicity'] = 0 # IN
        stat['total_sign_ethnicity'] = 0 # IN
        stat['total_nsign_ethnicity'] = 0 # IN
        stat['sum_variant'] = 0 # IN 2
        stat['sum_risk_allele'] = 0 # IN 
        stat['WARNING'] = {}
        stat['ERROR'] = {}

    stat['total_variant'] = 0 # IN 2
    stat['total_risk_allele'] = 0 # IN 
    
    ## sample info
    stat['total_sample'] = 0 # IN
    stat['sum_sample_size'] = 0
    stat['sum_sample_type'] = 0
    stat['total_sign_sample'] = 0 # IN
    stat['total_nsign_sample'] = 0 # IN
    
    ## model info
    stat['total_model'] = 0 # IN 00
    stat['total_sign_model'] = 0 # IN
    stat['total_nsign_model'] = 0 # IN
    stat['total_OK_OR_model'] = 0 # IN
    stat['total_NOTOK_OR_model'] = 0 # IN
    stat['total_OK_model'] = 0 # IN
    stat['total_NOTOK_model'] = 0 # IN
### END OF initialize_statistics

### input_filter ###
# filters out input for variant_scoring and LD_variant_selection
# reformat into dictionaries (i.e., no array)
# called by main
### 
def input_filter(indata):
    logging.debug(' Function: input_filter')

    ### Initialization 
    input_data = {} # gwas data -> variant_scoring
    input_info = {} # variant specific scores -> LD_variant_selection   
    
    STAT['variant_scoring'] = {}
    initialize_statistics(STAT['variant_scoring'], 0)

    ### Loop on variant documents/dictionaries
    for rs in indata:
        input_data[rs] = {}
        input_info[rs] = {}
        STAT['variant_scoring']['sum_variant'] += 1
        
        for key in indata[rs]:
            value = indata[rs][key]
            ## Records variant ID 
            if key.split("_")[1] == 'id':
                input_info[rs][key] = value
            
            ## Records nextbio score and pvalue
            if key.split("_")[0] == 'nextbio':
                input_info[rs][key] = value
            
            ##  records document ethnicity  
            elif key.split("_")[1] == 'ethnicity':
                input_data[rs][key] = value
                input_info[rs][key] = {}
                for key2 in indata[rs][key]:
                    STAT['variant_scoring']['total_ethnicity'] += 1
                    
                    ## update list of ethnicities
                    if key2 not in STAT['variant_scoring']['different_ethnicity']:
                        STAT['variant_scoring']['different_ethnicity'].append(key2)

                        ### INIT summary per ethnicity
                        STAT['variant_scoring'][key2] = {}
                        initialize_statistics(STAT['variant_scoring'][key2], 1 )
                    
                    ## record hapfreq in info dict for LD selection  
                    hapfreq = SF_class.convert_2_boolean(indata[rs][key][key2]['hapmap_freq'])
                    input_info[rs][key][key2] = hapfreq
                    
    return input_data, input_info
### END OF input_filter

### fecth_input_data ###
# - open an input file
# - convert the content from json to python dictionary
# called by main
### 
def fecth_input_data(filename, comment):
    logging.debug(' Function: fecth_input_data: filename: %s, comment %s' % (filename, comment))
    inputdata = {}
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

################ MAIN ###############
if __name__ == '__main__':
    
    ### default input files name
    file1 = 'variant_data.json'
    if len(sys.argv) >1:
        file1 =sys.argv[1]

    ### extract file name
    directory = file1.split('/')
    baseline = directory[len(directory)-1].split('.json')[0]
    
    ### otput file name
    file2 = baseline + '_LD.json'
    if len(sys.argv) >2:
        file2 =sys.argv[2]
    
    ### configure logs
    LOG_FILENAME=LOG_FILENAME+baseline
    logging.basicConfig(filename=LOG_FILENAME, filemode='w',level=logging.INFO,format='%(asctime)s - %(levelname)s -%(message)s')

    ### fetch input data
    input_data = fecth_input_data(file1,'variant data' )

    ### filter the input data into input for variant_scoring and LD_variant_selection
    (variant_data, variant_info) = input_filter(input_data)
    
    ### score variants
    variant_score = VS_class.variant_scoring(variant_data)
    
    ### intermediate output of variant scoring
    SF_class.output_data(variant_score, baseline+'_score.json','variant score' )
    
    ### fecth LD block information
    LD_data = fecth_input_data(file2,'LD data' )
        
    ### select best variants for each LD block
    variant_selection = LD_class.LD_variant_selection(variant_score,variant_info,LD_data)

    ### out best variant per LD block per ethnicity
    SF_class.output_data(variant_selection, baseline+'_selection.json','LD variant selection' )   

    ### output summary statistics
    summaries_disease()
    SF_class.output_data(STAT, baseline+'_summary.json','summaries' )  
    
    ########## TEST AND DEBUG ONLY #########
    SF_class.output_data(SF.WORKING, baseline+'_working.json','working dictionary' ) 

    print "I AM DONE WITH" , baseline

    pass
