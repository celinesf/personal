### @package config
# Config file for CoordinateTranslator.py 

"""
@file config.py
@author: Celine Becquet
@creation_date:  09/24/2016
"""

__author__ = "Celine Becquet"
__email__ = "cbecquet@ingenuity.com"
__status__ = "dev" 
__version__ = 1.0 

#---# libraries
import json
import logging
import sys
import os.path
import re

#---# custom modules
sys.path.append(os.path.abspath(os.path.dirname(__file__)))  
from config_functions import ConfigFunctions
utils = ConfigFunctions()

# from utils import ConfigVariables
# utils = ConfigVariables()

#---### Set Option for test purpose
# sys.argv.insert(1, "-h")
# sys.argv.insert(1, "--outdir")
# sys.argv.insert(2, "/Users/becquetc/ING_CODING/InvitaeCode/data/")
# #  
# sys.argv.insert(3, "-l")
# sys.argv.insert(4, "INFO")
#  
# sys.argv.insert(5, "-g")
# sys.argv.insert(6, "tete.txt")
#---### End set test options

#---### Default values for global variable 
### Name of file with transcript coordinates:\n 
# A four column (tab-separated) file containing the transcripts.
# The first column is the transcript name, 
# and the remaining three columns indicate it's genomic mapping: chromosome name, 
# 0-based starting position on the chromosome, and CIGAR string indicating the mapping.
###
TRANSCRIPT_FILENAME = "transcripts.txt" 

### Name of file with queries\n
# A two column (tab-separated) file indicating a set of queries.
# The first column is a transcript name, and the second column is a 0-based transcript coordinate.
###
QUERY_FILENAME = "queries.txt" 
 
### Name of the output file with the transcript coordinates translated to genome coordinates:\n
# A four column tab separated file with one row for each of the input queries.
# The first two columns are exactly the two columns from the second input file, 
# and the remaining two columns are the chromosome name and chromosome coordinate, respectively.
###
OUTPUT_FILENAME = "translation.txt" 
 
### Path to output file directory. 
# Default is "output_files" within the code package
###
OUTPUT_FILES_PATH = "" + (os.path.abspath(os.path.dirname(__file__)) + "/../output_files")
 
### Path to intput file directory. 
# Default is "input_files" within the code package
###
INPUT_FILES_PATH = "" + (os.path.abspath(os.path.dirname(__file__)) + "/../input_files")  

### Full path to output translation file
# Default is "/<PATH2PACKAGE>/output_files/translation.txt"
###
TRANSLATION_PATH = OUTPUT_FILES_PATH + "/" + OUTPUT_FILENAME

### Full path to input transcription file
# Default is "/<PATH2PACKAGE>/input_files/transcripts.txt"
###
TRANSCRIPT_FILE_PATH = INPUT_FILES_PATH + "/" + TRANSCRIPT_FILENAME

### Full path to input query file
# Default is "/<PATH2PACKAGE>/input_files/queries.txt"
###
QUERY_FILE_PATH = INPUT_FILES_PATH + "/" + QUERY_FILENAME 
 
### Logging file name
LOG_FILENAME = "logfile.log"

### Logging level: 'DEBUG', 'INFO', 'WARNING', 'ERROR' or 'CRITICAL'
LOG_LEVEL = "DEBUG"

### Verbosity level
# Default is VERBOSITY = 2: STDOUT format like log file: '%m/%d/%Y %H:%M:%S' - <file>.py:<class>.<function>:<line> - + <level> - <info>
# VERBOSITY = 1: STDOUT format (no time stamp): <file>.py:<class>.<function>:<line> - + <level> - <info>
# VERBOSITY = 0: no STDOUT
###
VERBOSITY = 2

#---### get and check options
utils.get_options(sys.argv[1:])

