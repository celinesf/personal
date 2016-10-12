###@package utility_functions
# Set of global utility functions

"""
@file utility_functions.py
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
import inspect
import re
import pickle
import gzip
import subprocess
import shlex
from pytz import timezone
from datetime import datetime  # @UnresolvedImport
#---# get custom modules
sys.path.append(os.path.abspath(os.path.dirname(__file__)))  
import config

### Set of global utility functions
class  UtilityFunctions:
    
    def __init__(self):
        """ UtilityFunctions class """
    #---# END __init__   
    
    ### print log to stdout 
    def warn_me(self, flag, comment, nolog=None):
        """ writes info to stdout and log file 
        Args:
            flag (str): Level of logging (flag=DEBUG, INFO, WARNING or ERROR).
            comment (str): Information to convey.
            nolog (:obj:`int`, optional): 1 if need to not log the info - used whenever warn_me is called prior to logging configuration.
        """
        file = inspect.currentframe().f_back.f_code.co_filename
        classname = inspect.currentframe().f_back.f_locals["self"].__class__.__name__
        line = inspect.currentframe().f_back.f_lineno
        function = inspect.currentframe().f_back.f_code.co_name
        flag = flag.lower() 
        file = re.sub(r'.*(InvitaeCode/)', '', file)
#         file = re.sub(r'.*(/lib/|bin/|InvitaeCode/)', '', file)
#         file = re.sub(r'bin/../', '', file)
        info = file + ":" + classname + "." + function + ":" + str(line) + " - " + flag.upper() + " - " + comment
        #---# start if output logs
        if(nolog == None):
            if flag == 'info':
                logging.info(info)  
            elif flag == 'warning':
                logging.warning(info)     
            elif flag == 'critical':
                logging.critical(info)  
            elif flag == 'error':
                logging.error(info)   
            else:
                logging.info(info)
        #---# END if output logs
        
        if ( flag != 'debug' and config.VERBOSITY != 0):
            """ STDOUT same info as in log file """
            if(config.VERBOSITY == 2):
                info = datetime.now(timezone('US/Pacific')).strftime('%m/%d/%Y %H:%M:%S') + " - " + info
            print(info)
        #---# END if print STDOUT log
    #---# END warn_me
         
    """ store_hash - no need? """   
    def store_hash(self, hash, filename, FILE_STREAM=None):
        """ pickle dump data into a binary/gziped storage file 
        Args:
            hash (dict): Dictionary to output.
            filename (str): Storage file full path.
            FILE_STREAM (:obj:`file`, optional): File Object to append hash to.
        """
        logging.debug("filename:'" + filename)
        FILE = FILE_STREAM
        if not FILE_STREAM:
            if  os.path.isfile(filename):
                command = "rm -f \"" + filename + "\""
                self.run_command_line(command, filename + "_rm_store_hash.txt")
            FILE = open_file(filename , "w", "gzip")
     
        if hash:
            try:
                pickle.dump(hash, FILE)
            except Exception:
                self.warn_me("CRITICAL", "Not able to pickle dump filename:'" + filename + "'")
                raise Exception ("Not able to to pickle dump filename:'" + filename + "'")
                sys.exit()
        else:
              self.warn_me("CRITICAL", "Hash is empty - Not able to pickle dump filename:'" + filename + "'")
              sys.exit()
        if (not  FILE_STREAM):
            FILE.close()
    #---# END store_hash
    
    """ retrieve_hash - no need?"""  
    def retrieve_hash(self, filename, FILE_STREAM=None) :
        """ load pickle binary/gziped data from a storage file 
        Args:
            filename (str): Storage file full path.
            FILE_STREAM (:obj:`file`, optional): File Object to extract hash to.
        Returns:
            A dictionary if exists, dies otherwise 
        """
        logging.debug("filename:'" + filename)

        FILE = FILE_STREAM
        if not FILE_STREAM:
            FILE = open_file(filename , "r", "gzip")
    
        hash = {}
        try:
            hash = pickle.load(FILE)
        except Exception:
           self.warn_me("CRITICAL", "Not able to pickle load filename:'" + filename + "'")
           raise Exception ("Not able to pickle load filename:'" + filename + "'")
           sys.exit()
        if (not FILE_STREAM):
            FILE.close()
           
        return hash;
    #---# END retrieve_hash
     
    ### run internal command line
    def run_command_line(self, command, filename):
        """ run shell command line 
        Args:
            command (str): command line to run 
            filename (str): full path of file to log stdout and stderr.
        Returns:
            String output by the command line if successful, dies otherwise
        """
        args = shlex.split(command)
    
        self.warn_me("info", "command: '" + command + "' filename:'" + filename + "'\n" + json.dumps(args, sort_keys=True, indent=4))
    
        proc = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (output, error) = proc.communicate()
        exitcode = proc.returncode
    
        if exitcode != 0 :
            self.warn_me("CRITICAL", " TRY AGAIN Failed: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
            proc = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (output, error) = proc.communicate()
            exitcode = proc.returncode
        if exitcode != 0 :   
            self.warn_me("CRITICAL", " Failed: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
            if (output):
                FILE = open_file(filename + "_FAILED_OUT.txt" , "w", "utf8")
                FILE.write(command + "\n")
                FILE.write(output.decode('utf8'))
                FILE.close()
            if (error):            
                FILE = open_file(filename + "_FAILED_ERROR.txt" , "w", "utf8")
                FILE.write(command + "\n")
                FILE.write(error.decode('utf8'))
                FILE.close()
            self.warn_me("CRITICAL", " Failed: exit code= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
            sys.exit()
        if (output):
            self.warn_me("INFO", "SUCCESS- I had stdout: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
            FILE = open_file(filename + "_success_OUT.txt" , "w", "utf8")
            FILE.write(command + "\n")
            FILE.write(output.decode('utf8'))
            FILE.close()
        if (error):
            self.warn_me("INFO", "SUCCESS- I had stderr: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
            FILE = open_file(filename + "_success_ERROR.txt" , "w", "utf8")
            FILE.write(command + "\n")
            FILE.write(error.decode('utf8'))
            FILE.close()
        return output
    #---# END run_command_line
    
    ### open a file in a specific mode 
    def open_file(self, filename, read_write, encoding=None) :
        """ open a file - remove it if already exists and in writing mode 
        Args:
            filename (str): Full path to file to open
            read_write (str): How to open the file: "r", "w" or "a" 
            encoding (str) optional: file encoding: "gzip" or "utf-8".
        Returns:
            A File object if successful, dies otherwise 
       
        Examples: attempt at pytest/ doctest.testmod TODO     
            >>> import config
            >>> from utility_functions import UtilityFunctions
            >>> utils = UtilityFunctions()
            
            Parameters 'filename', 'read_write', and 'encoding' is required
            >>> utils.open_file()
            Traceback (most recent call last):
              ...
            TypeError: open_file() missing 3 required positional arguments: 'filename', 'read_write', and 'encoding'

            read_write needs to be 'w','r' or'a'
            >>> config.VERBOSITY=0
            >>> utils.open_file( 'filename', 'read_write', 'encoding')
            Traceback (most recent call last):
            ...
            SystemExit

            >>> utils.open_file( config.TRANSLATION_PATH, 'r', 'utf8')
            <_io.TextIO ... encoding='utf-8'>
        """
        logging.debug("filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
        
        #---# check correct read_write
        if(read_write not in ('w','r','a')):
            self.warn_me("CRITICAL", "read_write needs to be 'w','r' or'a':'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
            sys.exit()
               
        #---# remove existing write file
        if re.search("w", read_write) and (os.path.islink(filename) or os.path.isfile(filename)):
            command = "rm -f \"" + filename + "\""
            self.run_command_line(command, filename + "_rm_open_file.txt")
        
        #---# open file with specific encoding and read_write
        if(encoding == "gzip"):
            try:
                f = gzip.open(filename, read_write)
            except Exception:
                self.warn_me("warning", "Not able to open gzip filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
                raise Exception ("Not able to open gzip filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
        elif(encoding == "utf8"):
            try:
                f = open(filename, read_write,encoding='utf-8')
            except Exception:
                self.warn_me("warning", "Not able to open utf8 filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
                raise Exception ("Not able to open utf8 filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
        else:
            try:  
                 f = open(filename, read_write)
            except Exception:
                self.warn_me("warning", "Not able to open text filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
                raise Exception ("Not able to open text filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
        #---# return file if possible
        if f:
            return f
    #---# END open_file
#---# END UtilityFunctions  


### main function
def main():
    """ test the functions """
    UF = UtilityFunctions()
#---# END main

if __name__ == "__main__":
    main()  
    
    """ test the functions """
    doctest.testmod(verbose=True)
#---# END main 