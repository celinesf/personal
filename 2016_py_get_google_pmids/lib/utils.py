'''
Created on Apr 5, 2016

@author: becquetc
'''
__author__ = "Celine Becquet"
__email__ = "cbecquet@ingenuity.com"
__status__ = "dev" 
import json
import logging
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
# import pytz
from datetime import datetime  # @UnresolvedImport
### get custom modules ###
sys.path.append( os.path.abspath(os.path.dirname(__file__)))  
import config


##################### START FUNCTIONS ####################
''' WarnMe '''   
def WarnMe(flag, comment):
    """ write info to stdput and log file """
    file = inspect.currentframe().f_back.f_code.co_filename
    line = inspect.currentframe().f_back.f_lineno
    function = inspect.currentframe().f_back.f_code.co_name
    flag = flag.lower() 
    logging.debug("file:'" + file + "' function:'" + function + "' line:'" + str(line) + "' - flag:'" + flag + "' - comment:'" + comment + "'")
    
    file = re.sub(r'/Users/becquetc/ING_CODING/GetGooglePmids/get_google_pmids/', '', file)
    file = re.sub(r'bin/../', '', file)
    
    if flag == 'info':
        logging.info("file:'" + file + "' function:'" + function + "' line:'" + str(line) + "' - flag:'" + flag + "' - comment:'" + comment + "'")  
    elif flag == 'warning':
        logging.warning("file:'" + file + "' function:'" + function + "' line:'" + str(line) + "' - flag:'" + flag + "' - comment:'" + comment + "'")     
    elif flag == 'critical':
        logging.critical("file:'" + file + "' function:'" + function + "' line:'" + str(line) + "' - flag:'" + flag + "' - comment:'" + comment + "'")  
    elif flag == 'error':
        logging.error("file:'" + file + "' function:'" + function + "' line:'" + str(line) + "' - flag:'" + flag + "' - comment:'" + comment + "'")   
    else:
        logging.info("file:'" + file + "' function:'" + function + "' line:'" + str(line) + "' - flag:'" + flag + "' - comment:'" + comment + "'")
    print(datetime.now(timezone('US/Pacific')).strftime('%m/%d/%Y %H:%M:%S') + "-" + file + ":" + function + ":" + str(line) + " - " + flag.upper() + " - " + comment)
#### END WarnMe

 
''' StoreHash '''   
def StoreHash(hash, filename, FILE_STREAM=None):
    logging.debug("filename:'" + filename)
    """ pickle dump data into a binary/gziped storage file """
    FILE = FILE_STREAM
    if not FILE_STREAM:
        if  os.path.isfile(filename):
            command = "rm -f \"" + filename + "\""
            RunCommandLine(command, filename + "_rm_StoreHash.txt" )
        FILE = OpenFile(filename , "w", "gzip")
 
    if hash:
        try:
            pickle.dump(hash, FILE)
        except Exception:
            WarnMe("CRITICAL", "Not able to pickle dump filename:'" + filename + "'")
            raise Exception ("Not able to to pickle dump filename:'" + filename + "'")
            sys.exit()
    else:
          WarnMe("CRITICAL", "Hash is empty - Not able to pickle dump filename:'" + filename + "'")
          sys.exit()
    if (not  FILE_STREAM):
        FILE.close()
#### END StoreHash

''' RetrieveHash '''  
def RetrieveHash(filename, FILE_STREAM=None):
    logging.debug("filename:'" + filename)
    """ load pickle binary/gziped data from a storage file """
     
    FILE = FILE_STREAM
    if not FILE_STREAM:
        FILE = OpenFile(filename , "r", "gzip")

    hash = {}
    try:
        hash = pickle.load(FILE)
    except Exception:
       WarnMe("CRITICAL", "Not able to pickle load filename:'" + filename + "'")
       raise Exception ("Not able to pickle load filename:'" + filename + "'")
       sys.exit()
    if (not FILE_STREAM):
        FILE.close()
        
    return hash;
#### END RetrieveHash
 
''' RunCommandLine ''' 
def RunCommandLine(command, filename):
    """ run shell command line """
    args = shlex.split(command)

    WarnMe("info", "command: '" + command + "' filename:'" + filename + "'\n" + json.dumps(args, sort_keys=True, indent=4))

    proc = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (output, error) = proc.communicate()
    exitcode = proc.returncode

    if exitcode != 0 :
        WarnMe("CRITICAL", " TRY AGAIN Failed: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
        proc = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (output, error) = proc.communicate()
        exitcode = proc.returncode
    if exitcode != 0 :   
        WarnMe("CRITICAL", " Failed: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
        if (output):
            FILE = OpenFile(filename + "_FAILED_OUT.txt" , "w", "utf8")
            FILE.write(command+"\n")
            FILE.write(output.decode('utf8'))
            FILE.close()
        if (error):            
            FILE = OpenFile(filename + "_FAILED_ERROR.txt" , "w", "utf8")
            FILE.write(command+"\n")
            FILE.write(error.decode('utf8'))
            FILE.close()
        WarnMe("CRITICAL", " Failed: exit code= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
        sys.exit()
    if (output):
        WarnMe("INFO", "SUCCESS- I had stdout: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
        FILE = OpenFile(filename + "_success_OUT.txt" , "w", "utf8")
        FILE.write(command+"\n")
        FILE.write(output.decode('utf8'))
        FILE.close()
    if (error):
        WarnMe("INFO", "SUCCESS- I had stderr: exitcode= '" + str(exitcode) + "' \nstdout: " + output.decode('utf8') + "\nstderr: " + error.decode('utf8') + "\n")
        FILE = OpenFile(filename + "_success_ERROR.txt" , "w", "utf8")
        FILE.write(command+"\n")
        FILE.write(error.decode('utf8'))
        FILE.close()
    return output
#### END RunCommandLine

''' OpenFile ''' 
def OpenFile(filename, read_write, encoding) :
    logging.debug("filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
    """ open a file - remove it if already exists and in writing mode """
    if re.search("w", read_write) and (os.path.islink(filename) or os.path.isfile(filename)):
        command = "rm -f \"" + filename + "\""
        RunCommandLine(command, filename + "_rm_OpenFile.txt" )
    if(encoding == "gzip"):
        try:
            f = gzip.open(filename, read_write)
        except Exception:
            WarnMe("warning", "Not able to open gzip filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
            raise Exception ("Not able to open gzip filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
    elif(encoding == "utf8"):
        try:
            f = open(filename, read_write)
        except Exception:
            WarnMe("warning", "Not able to open utf8 filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
            raise Exception ("Not able to open utf8 filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
    else:
        try:  
             f = open(filename, read_write, encoding=encoding)
        except Exception:
            WarnMe("warning", "Not able to open text filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
            raise Exception ("Not able to open text filename:'" + filename + "' read_write:'" + read_write + "' encoding:'" + encoding + "'")
    if f:
        return f
#### END OpenFile


"""GetOptions"""
def GetOptions(argv):
    logging.debug(str(argv))
#     print (str(argv))
    """ Get options from command line """
    try:
       opts, args = getopt.getopt(argv,"hu:t:d:v:",["uat_or_qc=","topic=","directory=","version=","help"])
    except getopt.GetoptError:
        WarnMe("warning", "wrong option")
        GetHelp()
 
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            GetHelp()
        elif opt in ("-u", "--uat_or_qc"):
            config.UATorQC = arg.upper()
            if(not config.UATorQC in ("UAT","QC")):
                WarnMe("CRITICAL", "'"+config.UATorQC+"' is not a valid option. You need to specify whether you want to extract 'UAT' or 'QC' pmids using option -u UAT_OR_QC, --uat_or_qc UAT_OR_QC")
                GetHelp()
        elif opt in ("-t", "--topic"):
            config.topic = arg
        elif opt in ("-f", "--directory"):
            config.folder_id = arg
        elif opt in ("-v", "--version"):
            config.version = arg
#     print ("uat_or_qc " +config.UATorQC+ " topic " +config.topic+" directory " +config.folder_id+ " version " +config.version)
    #### CHECK ALL REQUIRED OPTIONS HERE
    info="Extraction PMIDs for: "
    if(config.UATorQC == ""):
        WarnMe("CRITICAL", "You need to specify whether you want to extract UAT or QC pmids using option -u UAT_OR_QC, --uat_or_qc UAT_OR_QC")
        GetHelp()
    else:
        info=info+"'"+config.UATorQC+"'"
        if(config.topic == "" ):
            if( config.UATorQC == "UAT"):
                WarnMe("CRITICAL", "for UAT, You need to specify a topic to extract pmids for using option -t TOPIC, --topic TOPIC")
                GetHelp()
            elif(config.UATorQC == "QC"):
                if(config.folder_id == ""):
                    WarnMe("CRITICAL", "You need to specify the either the folder ID for the all QC batch to extract pmids for, using option -d DIR_ID, --directory DIR_ID")
                    GetHelp() 
                elif (config.version != ""):
                    WarnMe("CRITICAL", "You need to specify the topic of the QC batch to extract pmids for, using option -t TOPIC, --topic TOPIC")
                    GetHelp()  
                else:
                    info=info+" for all google sheets in"  
                    info=info+" Folder: '"+config.folder_id+"'" 
        else:
            if(config.UATorQC == "UAT"):
                if( config.folder_id == "" ):
                    WarnMe("CRITICAL", "You need to specify the folder ID for the UAT round you want to extract pmids for, using option -d DIR_ID, --directory DIR_ID")
                    GetHelp()
                else:
                    info=info+", Topic: '"+config.topic+"'"
                    info=info+", Folder: '"+config.folder_id+"'" 
            elif(config.UATorQC == "QC"):        
                if( config.folder_id == "" ):
                    if (config.version == ""):
                        WarnMe("CRITICAL", "You need to specify the QC the version 'vxxx' of the QC file to extract pmids for, using option -v VERSION, --version VERSION")
                        GetHelp()
                    else:
                        config.folder_id="0B__wQSJJ3xYTfkduQkFSLU9Idk5HRDlWOTNQN0tUYnJoR2xlZHZySFZ3akdORzdzNFkwMTg"
                        info=info+", Folder: '"+config.folder_id+"'" 
                        info=info+", Topic: '"+config.topic+"'"
                        info=info+", Version: '"+config.version+"'"
                elif( config.folder_id != "" ):
                    if (config.folder_id != "0B__wQSJJ3xYTfkduQkFSLU9Idk5HRDlWOTNQN0tUYnJoR2xlZHZySFZ3akdORzdzNFkwMTg"):
                        WarnMe("CRITICAL", "The folder_id you specified isn't the one expected for Clinical QC batchs +'"+config.folder_id+"'")
                        GetHelp()
                    else:
                        if (config.version == "" ):
                            config.topic=""
                            info=info+". I will ignore, Topic: '"+config.topic+"'"
                            info=info+", and run for all QC batches in Folder: '"+config.folder_id+"'" 
                        else:
                             info=info+", Folder: '"+config.folder_id+"'" 
                             info=info+", Topic: '"+config.topic+"'"
                             info=info+", version: '"+config.version+"'" 
                
    WarnMe("INFO", info)
      
#### END GetOptions

"""GetHelp"""
def GetHelp():
    logging.debug("")
    """ Get help """
    print("USAGE")
    print("\tget_google_pmids.py --uat_or_qc UAT_OR_QC --topic TOPIC [OPTIONS] ")
    print("DESCRIPTION") 
    print("\tExtracts pmids from the specified google QC or UAT rounds")
    
    print("OPTIONS")
    print("\t\t-u UAT_OR_QC, --uat_or_qc UAT_OR_QC")
    print("\t\t\tUAT_OR_QC the type of QC to extract pmids for. UAT_OR_QC can be 'UAT' or 'QC' only. This option is required ")
    
    print("\t\t-t TOPIC, --topic TOPIC")
    print("\t\t\tTOPIC specifies the topic for which to extract pmids. This option is required")
    
    print("OPTIONS required if UAT_OR_QC='UAT'")
    print("\t\t-d DIR_ID, --directory DIR_ID")
    print("\t\t\tDIR_ID specifies the folder id for the UAT round for which to extract pmids. This option is required if UAT_OR_QC='UAT'")
    
    print("OPTIONS required if UAT_OR_QC='QC'")
    print("\t\t-v VERSION, --version VERSION")
    print("\t\t\tFVERSION specifies the version i.e., '001' in the file name for the QC round for which to extract pmids. This option is required if UAT_OR_QC='QC'")
    
    print("OTHER OPTIONS")
    print("\t\t-h, --help")
    print("\t\t\tOutput this help manual")
    sys.exit(2)
#### END GetHelp

"""InitPmidData"""   
def InitPmidData():
    logging.debug("")
    """ retrieve or init the PMID list """  
    ### retrieve/init pmid list
    if (os.path.isfile(config.pmid_list_store_path)): 
        PMID_LISTS=RetrieveHash(config.pmid_list_store_path)
    else:
        PMID_LISTS = {"PMIDS":{}}
    return PMID_LISTS
#### END InitPmidData

"""InitTopicData"""   
def InitTopicData(topic,PMID_LISTS):
    logging.debug("")
    """ init the PMID list for this folder/topic """  
    ### Init PMID_LISTS for this topic
    if(topic != "" and not topic in PMID_LISTS):
        PMID_LISTS[topic]={}
    if(topic != "" and not config.UATorQC in PMID_LISTS[topic]):
        PMID_LISTS[topic][config.UATorQC]={}
    return PMID_LISTS
#### END InitTopicData

"""GetPmidsQCCode"""   
def GetPmidsQCCode(nrow,pmid, info):
    logging.debug("nrow:'"+ str(nrow) + "' pmid:'" + pmid + "' info:'" + info + "'")
    int_pmid=0
    int_info=0
    pattern_content=re.compile(r'\#|\*')
    pmid=pattern_content.sub("",pmid)
    try: 
        int_pmid =int(pmid)
    except Exception: 
        if(pmid != "" and not re.search("total", pmid, re.IGNORECASE)):
            WarnMe("CRITICAL", "At row: '"+str(nrow)+"' NOT A PMID: '"+pmid+ "' and info is '"+info+ "'")  
            sys.exit()
    if(int_pmid>0 ):
        try: 
            int_info =int(info)
        except Exception: 
            if(info !="NC" and info !="" and info !="1b"):
                WarnMe("CRITICAL", "At row: '"+str(nrow)+"' pmid '"+pmid+ "' has no info '"+info+ "'")  
                sys.exit()
        if(int_pmid>0 and int_info in config.QC_FALSE_POSITVE_CODE):
            return int_info
        else:
            return 0
#### END GetPmidsQCCode