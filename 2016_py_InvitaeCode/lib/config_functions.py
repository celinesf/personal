###@package config_functions
# Set of functions to configure global variable from command line options

"""
@file config_functions.py
@author: Celine Becquet
@creation_date:  09/24/2016
"""

__author__ = "Celine Becquet"
__email__ = "celine.becquet@gmail.com"
__status__ = "dev" 
__version__ = 1.0 

#---# libraries
import json
import logging
import doctest
import sys
import getopt
import os.path
import re
from pytz import timezone
from datetime import datetime  # @UnresolvedImport
#---# get custom modules
sys.path.append(os.path.abspath(os.path.dirname(__file__)))  
import config
from utility_functions import UtilityFunctions


### Set of functions to configure global variable from command line options
class ConfigFunctions:
     
    def __init__(self): 
        """ ConfigFunctions Class 
        Attributes:
            _utils (UtilityFunctions): Set if utility functions.
        """   
        self._utils = UtilityFunctions()   
    #---# END __init__
     
    """  _configure_logs """
    ### configures the log file
    def _configure_logs(self):
        """ Configure the logging file """
        logging.basicConfig(filename=config.OUTPUT_FILES_PATH + 
        "/" + config.LOG_FILENAME, filemode='w', level=config.LOG_LEVEL,
        format='%(asctime)s - %(levelname)s - %(filename)s:%(funcName)s:%(lineno)s - %(message)s')
        self._utils.warn_me("INFO", "Configured the logging file with level '" + config.LOG_LEVEL + "' here:" + config.OUTPUT_FILES_PATH + "/" + config.LOG_FILENAME)   
    #---# END _configure_logs 
     
    ### sanity check of the command line arguments 
    def _check_options(self, type, display, variable, nolog=None):
        """ Checks that options are sound 
         Args:
            type (str): Type of variable to check:  "directory" or "file"
            display (str): Display value in case variable is not sound
            variable (str): Variable to check.
            nolog (:obj:`int`, optional): 1 if need to not log the info.
            
        Examples: # attempt at pytest/ doctest.testmod     
            >>> import config
            >>> from utility_functions import UtilityFunctions
            >>> from config_functions import ConfigFunctions
            >>> utils = ConfigFunctions()
            
            The string arguments 'type', 'display', and 'variable' are required
            >>> utils._check_options()
            Traceback (most recent call last):
              ...
            TypeError: _check_options() missing 3 required positional arguments: 'type', 'display', and 'variable'
            
            Type='directory' or 'file'
            >>> utils._check_options('type', 'display', 'variable')
            Traceback (most recent call last):
            ...
            SystemExit
            
            Correct type and existing directory
            >>> config.VERBOSITY=0
            >>> utils._check_options('directory', 'display', config.OUTPUT_FILES_PATH)           
        """
        self._utils.warn_me("INFO", "Checking options type: '" + type + "' display: '" + display + "' variable: '" + variable + "'", nolog) 
         
        found = 0
        #---# Check if directory exists
        if(type == "directory"): 
            if (os.path.isdir(variable)): 
                found = 1
        elif(type == "file"): 
            if (os.path.isfile(variable)): 
                found = 1
        else:
            self._utils.warn_me("CRITICAL", "The type '" + type + "' doesn't exist - value should be 'directory' or 'file'", nolog)
            self._get_help()
        if(not found):
            self._utils.warn_me("CRITICAL", "The '" + display + "' '" + type + "' specified is incorrect or doesn't exist: '" + variable + "'", nolog)
            self._get_help()
    #---# END _check_options 
        

    ### Extract options from the command line arguments
    def get_options(self, argv):
        """ Extracts options from the command line 
            and calls _configure_logs 
            and check soundness of variables
        Args:
            argv (list): List of options from command line

        Examples: attempt at pytest/ doctest.testmod TODO     
            >>> import config
            >>> from utility_functions import UtilityFunctions
            >>> from config_functions import ConfigFunctions
            >>> utils = ConfigFunctions()
            
            Parameter argv is required
            >>> utils.get_options()
            Traceback (most recent call last):
              ...
            TypeError: get_options() missing 1 required positional argument: 'argv'
            
            No options is OK
            >>> config.VERBOSITY=0
            >>> utils.get_options([])
            
            Get Manual and exit
            >>> utils.get_options(["-h"])
            Traceback (most recent call last):
            ...
            SystemExit
        """
        self._utils.warn_me("INFO", str(argv), 1)
         
        #---# get options
        try:
            opts, args = getopt.getopt(argv, "t:q:f:i:o:hg:l:v:", 
                                       ["transcripts=", "queries=", "outfile=", "indir=", "outdir=", "help", "logfile", "loglevel","verbosity"])
        except getopt.GetoptError:
            self.warn_me("CRITICAL", "wrong option", 1)
            self._get_help()
    
        #---# Start option loop
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self._get_help()
            elif opt in ("-1", "--transcripts"):
                config.INPUT1_FILENAME = arg
            elif opt in ("-2", "--queries"):
                config.INPUT2_FILENAME = arg
            elif opt in ("-o", "--output"):
                config.OUTPUT_FILENAME = arg
            elif opt in ("-i", "--indir"):
                config.INPUT_FILES_PATH = arg 
                config.INPUT_FILES_PATH = re.sub(r'/$', "", config.INPUT_FILES_PATH)
            elif opt in ("-d", "--outdir"):
                config.OUTPUT_FILES_PATH = arg
                config.OUTPUT_FILES_PATH = re.sub(r'/$', "", config.OUTPUT_FILES_PATH)
                self._check_options("directory", "output", config.OUTPUT_FILES_PATH, 1)
            elif opt in ("-g", "--logfile"):
                config.LOG_FILENAME = arg
            elif opt in ("-l", "--loglevel"):
                config.LOG_LEVEL = arg.upper()
                if(config.LOG_LEVEL not in ("DEBUG","INFO","WARNING","ERROR","CRITICAL")):
                    self._utils.warn_me("CRITICAL", "for option -l or --loglevel LOGLEVEL, accepted values are 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL' only",1)
                    self._get_help()
            elif opt in ("-l", "--verbosity"):
                config.VERBOSITY = arg
        #---# END option loop
        self._configure_logs()
         
        #---# Configure full paths to input and output files
        config.TRANSLATION_PATH = config.OUTPUT_FILES_PATH + "/" + config.OUTPUT_FILENAME
        config.TRANSCRIPT_FILE_PATH = config.INPUT_FILES_PATH + "/" + config.TRANSCRIPT_FILENAME
        config.QUERY_FILE_PATH = config.INPUT_FILES_PATH + "/" + config.QUERY_FILENAME 
         
        #---# Check input directory and files
        self._check_options("directory", "input", config.INPUT_FILES_PATH)
        self._check_options("file", "transcript", config.TRANSCRIPT_FILE_PATH)
        self._check_options("file", "query", config.QUERY_FILE_PATH)       
 
        self._utils.warn_me("INFO", "DONE getting/checking options")
    #--- END get_options
     
    ### output the manual for coordinate_translator.py
    def _get_help(self):
        """ Output help document """ 
        print("USAGE")
        print("\tcoordinate_translator.py [OPTIONS] ")
        print("DESCRIPTION") 
        print("\tTranslates (0-based) transcript coordinates to (0-based) genome coordinates")
          
        print("OPTIONS")
        print("\t\t-t FILENAME, --transcripts FILENAME")
        print("\t\t\tFILENAME specifies the name of a four column (tab-separated) file containing the transcripts.")
        print("\t\t\tThe first column is the transcript name, and the remaining three columns indicate it's genomic mapping: ")
        print("\t\t\tchromosome name, 0-based starting position on the chromosome, and CIGAR string indicating the mapping. ")
        print("\t\t\tBy default, FILENAME='transcripts.txt'")
          
        print("\t\t-q FILENAME, --queries FILENAME")
        print("\t\t\tFILENAME specifies the name of a two column (tab-separated) file indicating a set of queries.")
        print("\t\t\tThe first column is a transcript name, and the second column is a 0-based transcript coordinate.")
        print("\t\t\tBy default, FILENAME='queries.txt'")
    
        print("\t\t-i PATH, --indir PATH")
        print("\t\t\tPATH specifies the full path to the input files. By default, PATH=<PATH2PACKAGE>/input_files")
    
        print("\t\t-f FILENAME, --outfile FILENAME")
        print("\t\t\tFILENAME specifies the name of a four column tab separated file with one row for each of the input queries.")
        print("\t\t\tThe first two columns are exactly the two columns from the query file.")
        print("\t\t\tThe remaining two columns are the chromosome name and chromosome coordinate, respectively.")
        print("\t\t\tBy default, FILENAME='translation.txt'")
        
        print("\t\t-o PATH, --outdir PATH")
        print("\t\t\tPATH specifies the full path to the output and log files.")
        print("\t\t\tBy default, PATH=<PATH2PACKAGE/output_files>")
          
        print("\t\t-r FILENAME, --logfile FILENAME")
        print("\t\t\tFILENAME specifies the name of the log file.")
        print("\t\t\tBy default, FILENAME=logfile.log")
          
        print("\t\t-l LOGLEVEL, --loglevel LOGLEVEL")
        print("\t\t\tLOGLEVEL specifies the logging level (i.e. controls the verbosity) in the logging file.")
        print("\t\t\tLOGLEVEL='DEBUG', 'INFO', 'WARNING', 'ERROR' or 'CRITICAL' (by default LOGLEVEL=DEBUG)")
  
        print("\t\t-v VERBOSITY, --verbosity VERBOSITY")
        print("\t\t\tVERBOSITY= 0, 1 or 2 (by default VERBOSITY=2)")
        print("\t\t\tDefault is VERBOSITY = 2: STDOUT format like log file: '%m/%d/%Y %H:%M:%S' - <file>.py:<class>.<function>:<line> - + <level> - <info>.")
        print("\t\t\tVERBOSITY = 0: STDOUT format (no time stamp): <file>.py:<class>.<function>:<line> - + <level> - <info>.")
        print("\t\t\tIf VERBOSITY = None: no STDOUT.")
 
        print("OTHER OPTIONS")
        print("\t\t-h, --help")
        print("\t\t\tOutput this help manual")
        exit()
#         sys.exit()
    #--- END _get_help 
#--- END ConfigFunctions 

    
### main function
def main():
    """ test the functions """
    doctest.testmod(verbose=True)
    CF = ConfigFunctions()
#---# END main
 
if __name__ == "__main__":
    main()      

#---# END main    

