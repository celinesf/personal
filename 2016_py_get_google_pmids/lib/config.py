'''
Config files for get_google_pmids 
Created on Apr 5, 2016

@author: Celine Becquet
'''
__author__ = "Celine Becquet"
__email__ = "cbecquet@ingenuity.com"
__status__ = "dev" 

import json
import logging
import sys
import os.path
### end import ###
# sys.path.append(os.path.abspath(os.path.dirname(__file__)))  

print("Current working dir : %s" % os.getcwd())

#### list of options ###
topic = ""
folder_id = ""
UATorQC = ""
version = ""

sys.argv.insert(1, "--uat_or_qc")

# sys.argv.insert(2, "UAT")
# sys.argv.insert(3, "--topic")

# sys.argv.insert(4, "BRCA")
# sys.argv.insert(5, "--directory")
# sys.argv.insert(6, "0B0WrZ2XZYv2ncEgzMlBRV2RCUGs") #  DONE 140  "0B0WrZ2XZYv2ncEgzMlBRV2RCUGs":"BRCA QC" ,05/20

# sys.argv.insert(4, "UCS")
# sys.argv.insert(5, "--directory")
# sys.argv.insert(6, "0B8R_OxurEtHFUThPbVZSZHpEZGs") #TODO WAIT  57 "0B8R_OxurEtHFUThPbVZSZHpEZGs" :"UCS UAT Scoring Sheets",  

# sys.argv.insert(4, "Hereditary Cancer v2")
# sys.argv.insert(5, "--directory")
# sys.argv.insert(6, "0B5vkaWupAx_vTHJBMzR6M19ZcDg") #DONE  50  0B5vkaWupAx_vTHJBMzR6M19ZcDg": "QC round 1", 
# sys.argv.insert(6, "0B5vkaWupAx_vM0F5MmYyVS1LVFk") #DONE 49   "0B5vkaWupAx_vM0F5MmYyVS1LVFk": "QC round 2", 
# sys.argv.insert(6, "0B5vkaWupAx_vSmIxdFRhSTllZE0")  #DONE 50  0B5vkaWupAx_vSmIxdFRhSTllZE0: "QC round 3", 
# sys.argv.insert(6, "0B5vkaWupAx_vVWxlaC1CYlNmTFU")  #DONE 32  0B5vkaWupAx_vVWxlaC1CYlNmTFU: "QC round 4"
# sys.argv.insert(6, "0B5vkaWupAx_vaW42ZHRURG9IRVE")  #TODO WAIT 52  0B5vkaWupAx_vaW42ZHRURG9IRVE: "QC round 5",
# sys.argv.insert(6, "0B5vkaWupAx_vRmRsbThnWm8xdms")  #TODO WAIT 50  0B5vkaWupAx_vRmRsbThnWm8xdms: "QC round 6",

# sys.argv.insert(4, "CGS") 
# sys.argv.insert(5, "--directory")
# sys.argv.insert(6, "0B8R_OxurEtHFRzZNZmh4R05UUHM") #DONE 74 "0B8R_OxurEtHFRzZNZmh4R05UUHM": "CGS UAT Scoring Sheets",
# sys.argv.insert(6, "0B8quviPZBKo7cU94b0d6Q3pVN1U") #DONE 50 "0B8quviPZBKo7cU94b0d6Q3pVN1U": "CGS2 Bibliography UAT round 2"
# sys.argv.insert(6, "0B8quviPZBKo7Qm5RTlJuZlYtQ0k") #DONE 50  "0B8quviPZBKo7Qm5RTlJuZlYtQ0k":"CGS2 Bibliography UAT round 3",

# sys.argv.insert(4, "Hereditary Cancer v1") DONE
# sys.argv.insert(5, "--directory")
# sys.argv.insert(6, "0B8R_OxurEtHFfkVxTnltLUhQQTRRSUp1M1pSNkZlelEzd0VaRk45YVg2SlVodHVCSEFCYkE") #DONE   11 "0B8R_OxurEtHFfkVxTnltLUhQQTRRSUp1M1pSNkZlelEzd0VaRk45YVg2SlVodHVCSEFCYkE":"HerCan QC Previous Scoring Sheets",Hereditary Cancer QC Summary
# sys.argv.insert(6, "0B8R_OxurEtHFfnV0S0ZNT3dtTi14bkw5SjdOcWNqQW1JaF9YZ3FPRHpEdF83WkRZTTV4RUE") #DONE  39 "0B8R_OxurEtHFfnV0S0ZNT3dtTi14bkw5SjdOcWNqQW1JaF9YZ3FPRHpEdF83WkRZTTV4RUE":"HerCan QC Round 1 Scoring Sheets",
# sys.argv.insert(6, "0B8R_OxurEtHFfkNkNU12VlFtTzRkYmp5RmpScUV6QWVpaFQ2Y1BlVGlLWnZKNDQxLXNFM2c") #DONE 34   "0B8R_OxurEtHFfkNkNU12VlFtTzRkYmp5RmpScUV6QWVpaFQ2Y1BlVGlLWnZKNDQxLXNFM2c":"HerCan QC Round 2 Scoring Sheets",
# sys.argv.insert(6, "0B8R_OxurEtHFaWVaMDB4WDJCTUk") #DONE  20  "0B8R_OxurEtHFaWVaMDB4WDJCTUk":"Her Can QC Round 3 Scoring Sheets",
# sys.argv.insert(6, "0B8R_OxurEtHFNE9ucXY2YWFoeEE") #DONE   25 "0B8R_OxurEtHFNE9ucXY2YWFoeEE":"Her Can QC Round 4 Scoring Sheets",                        


###### QC
sys.argv.insert(2, "QC")
sys.argv.insert(3, "--directory")
sys.argv.insert(4, "0B__wQSJJ3xYTfkduQkFSLU9Idk5HRDlWOTNQN0tUYnJoR2xlZHZySFZ3akdORzdzNFkwMTg")
#    
sys.argv.insert(3, "--topic")
# sys.argv.insert(4, "Hereditary Cancer v1")
# sys.argv.insert(4, "Trusight cancer")
sys.argv.insert(4, "ACMG")
# sys.argv.insert(4, "CGS")
# sys.argv.insert(4, "CGS")
#    
sys.argv.insert(5, "--version")
sys.argv.insert(6, "v026")

# sys.argv.insert(7, "--help")

# HERE Clinical QC-coverage-v001-lynch-100
# HERE Clinical QC-coverage-v002-non-brca-50
# HERE Clinical QC-coverage-v003-non-brca-50
# HERE Clinical QC-coverage-v004-non-brca-50
# HERE Clinical QC-coverage-v005-non-brca-50
# HERE Clinical QC-coverage-v006-ucs-50
# HERE Clinical QC-coverage-v007-ucs-50
# HERE Clinical QC-coverage-v008-ucs-50
# HERE Clinical QC-coverage-v009-ucs-50
# HERE Clinical QC-coverage-v010-cgs2-50
# HERE Clinical QC-coverage-v011-cgs2-50
# HERE Clinical QC-coverage-v012-cgs2-50
# HERE Clinical QC-coverage-v013-ucs-50
# HERE Clinical QC-coverage-v014-ucs-50
# HERE Clinical QC-coverage-v015-ucs-50
# HERE Clinical QC-coverage-v016-ucs-hb-50
# HERE Clinical QC-coverage-v017-ucs-50
# HERE Clinical QC-coverage-v018-cgs2-50
# HERE Clinical QC-coverage-v019-cgs2-50
# HERE Clinical QC-coverage-v020-cgs2-50
# HERE Clinical QC-coverage-v021-heredv2-50
# HERE Clinical QC-coverage-v022-heredv2-50
# HERE Clinical QC-coverage-v023-heredv2-50
# HERE Clinical QC-coverage-v024-heredv2-50-abandoned
# HERE Clinical QC-coverage-v024-heredv2-50-updated
# HERE Clinical QC-coverage-v025-HCv2-50
# HERE Clinical QC-coverage-v026-HCv2-50-deprecate
# HERE Clinical QC-coverage-vXXX-template-50


OUTPUTPATH = (os.path.abspath(os.path.dirname(__file__)) + "/../data")  # TODO to change

LOG_FILENAME = "logfile.log"
credentials_drive_filename = "CB_google_drive_credential_file.json"
credentials_sheet_filename = credentials_drive_filename  
full_drive_filename = "drivefull_store"
tree_drive_filename = "drivetree"
pmid_list_filename = "pmid_list"
pmid_issues_filename = "pmid_issues"
qcdate_filename = "QC_date"

full_drive_path = OUTPUTPATH + "/" + full_drive_filename + ".txt"
tree_drive_path = OUTPUTPATH + "/" + tree_drive_filename + ".txt"
# tree_drive_store_path = OUTPUTPATH + "/" + tree_drive_filename + ".pklz"

pmid_list_path = OUTPUTPATH + "/" + pmid_list_filename + ".txt"
pmid_list_store_path = OUTPUTPATH + "/" + pmid_list_filename + ".pklz"

pmid_issues_path = OUTPUTPATH + "/" + pmid_issues_filename + ".txt"
qcdate_file_path = OUTPUTPATH + "/" + qcdate_filename + ".txt"
qc_file_date = ""

UAT_SUMMARY_NAMES = {
"BRCA":{
        "10LiTFULqap2dxzOjsvUzsZxuepc9VIU0zE6eng2VHlA": "Final BRCA UAT Summary", 
       },
"UCS":{
       "1CW3KWTkyPkj3060EqFa7qRLlPngpIKpCxuwkv-6uelc": "UCS UAT Coverage Summary",
       },
"Hereditary Cancer v2":{          
                        "1o1EFr-YJOyoSz_-a7A_XTOw-w5p196XIvDUbSMig4Kw": "Hereditary Cancer v2 bibliography summary",
                        },
"Hereditary Cancer v1":{                     
                       "1sDT1blJlLf4xQY_jmp983B5-pw9mXMOSQc-ReogcdT4":"HerCan nonLynch nonBRCA UAT Summary", 
                       "1y8L-ePp7bf_9TKst8UDukeVh7Lojmzdu9M-wEt5aW5k":"Hereditary Cancer QC Summary", 
                       },
"CGS":{
       "1DGYDKVtvmMC7XdQ5aGRARS6UdARQV6iy4Dw3irCRxO8":"CGS UAT Summary", 
       }
}

UAT_TOPIC_NAMES = {"BRCA":"BRCA",
            "Hereditary Cancer v2":"Hereditary Cancer v2",
            "Hereditary Cancer v1":"HerCan|Her Can",
            "UCS":"UCS",
            "CGS":"CGS|CGS2",
            }

QC_TOPIC_NAMES = [
            {"topic":"Hereditary Cancer v1",
             "regex":"non-brca|lynch"},
            {"topic":"Hereditary Cancer v2",
             "regex":"heredv2|HCv2"},
           {"topic": "UCS",
             "regex":"UCS"},
           {"topic": "CGS",
             "regex":"CGS|CGS2"},
           {"topic": "BRCA",
             "regex":"BRCA"},
            {"topic": "Trusight cancer",
             "regex":"TrusightACMGnewonly"},
        {"topic": "ACMG",
             "regex":"TrusightACMGnewonly"},
            ]
QC_REGEX = ""
for t in QC_TOPIC_NAMES:
    regex = t["regex"]
    if(QC_REGEX == ""):
        QC_REGEX = regex
    else:
        QC_REGEX = QC_REGEX + "|" + regex
print(QC_REGEX)

QC_REPORT_TAB = ["report", "report-after IKRS update"]
QC_CONTENT_TAB = ["Content report - final", "content report-final", "Content report-final",
                  "Content report"]
QC_REPORT_COL = {
               9:{
                "ncol":4,
                "header": {"Batch Coverage PMIDs":1,
                     "Batch 2 PMIDs":1,
                     "All PMIDs from Online resources":1}
                },
                8:{
                    "ncol":4,
                    "header":{"Batch 1 PMIDs":1} 
                    }
               }  

QC_PMID_COL = {
             2:{"PMID missing in KB":1,
                "QC-only PMIDs":1,  # remove * and # in front of pmids
                },
            3:{"PMIDs missing in KB":1,
               "PMID missing in KB":1,
                "PMIDs missing in KB (unique)":1,
                },
            6:{"PMIDs missed in KB (Batch 1 specific)":1,
                }
             }
QC_FALSE_POSITIVE_COL = {
                         8:{"Error category":1, "Code":1, "error category":1
                            } ,
                         6:{"Error category":1} ,
                         5:{"Error category":1} ,
                         4:{"Error type":1} ,
                      }

QC_FALSE_POSITVE_CODE = {
                  5,  # #falseP
6,  # # other language
13,  # # out of scope
15,  # #out of scope for future prioritization (nice-to-have right now)
}


#### logging info
logging.basicConfig(filename=OUTPUTPATH + "/" + LOG_FILENAME, filemode='w',
level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(filename)s:%(funcName)s:%(lineno)s - %(message)s')

