'''
get_google_pmids gets lists of pmids from QC and uat google sheets
Created on Apr 5, 2016

@author: becquetc
'''
__author__ = "Celine Becquet"
__email__ = "cbecquet@ingenuity.com"
__status__ = "dev" 

import json
import logging
import sys
import os.path
# import inspect
import re
import gc
### get custom modules ###
sys.path.append(re.sub("bin", "", os.path.abspath(os.path.dirname(__file__))) + "lib")  
import config  
import utils 
import qc_utils
import uat_utils
# print(json.dumps(sys.path , sort_keys=True, indent=4))
#     proc = psutil.Process(os.getpid())
#     gc.collect()
#     gc.collect()
#     print("C:" + str(proc.memory_info().rss))

  
"""OutputPmidList"""
def OutputPmidList(PMID_LISTS,DRIVE_TREE):
    logging.debug("")
    """ output text file for pmid list """  
    FILE = utils.OpenFile(config.pmid_list_path, 'w',"utf8")
    FILE.write("pmid\ttopic\tUAT or QC\tfolder_name\tfolder_id\tfile_name\tfile_id\n")
    for topic in (PMID_LISTS):
        if(topic !="PMIDS"):
            for UAT_or_QC in (PMID_LISTS[topic]):
                for folder_id in (PMID_LISTS[topic][UAT_or_QC]):
                    folder_name= DRIVE_TREE[folder_id]["title"]
                    print(topic+" "+UAT_or_QC+" "+folder_id+" "+folder_name)
                    for pmid in (PMID_LISTS[topic][UAT_or_QC][folder_id]):
                        for file_id in (PMID_LISTS[topic][UAT_or_QC][folder_id][pmid]):
                            file_name = PMID_LISTS[topic][UAT_or_QC][folder_id][pmid][file_id]
                            FILE.write(pmid+"\t")
                            FILE.write(topic+"\t")
                            FILE.write(UAT_or_QC+"\t")
                            FILE.write(folder_name+"\t")
                            FILE.write(folder_id+"\t")
                            FILE.write(file_name+"\t")
                            FILE.write(file_id+"\n")              
    FILE.close()  
#### END OutputPmidList
   
def main(argv):
    """ Main function """

    ##################
    ################## MAIN START #######################
    logging.debug('starting')

    utils.GetOptions(argv[1:])

    """ loop on topics """
    if(config.UATorQC == "UAT"):
        (PMID_LISTS,DRIVE_TREE) = uat_utils.GetUATPmids()
    elif(config.UATorQC == "QC"):
        (PMID_LISTS,DRIVE_TREE) = qc_utils.GetQCPmids()
        
    OutputPmidList(PMID_LISTS,DRIVE_TREE)

if __name__ == "__main__":
   
    main(sys.argv)  
