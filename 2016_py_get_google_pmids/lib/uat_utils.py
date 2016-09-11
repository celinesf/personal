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

import copy
import re
### get custom modules ###
sys.path.append( os.path.abspath(os.path.dirname(__file__)))  
import config  
import utils 
import google_utils 


##################### START FUNCTIONS ####################
"""GetUATPmids"""
def GetUATPmids():
    logging.debug("")
    """GetUATPmids: extracts pmids from UAT files for each topics """
    print(json.dumps(config.UAT_SUMMARY_NAMES, sort_keys=True, indent=4))

    ### get folder file list   
    (SHEETS,DRIVE_TREE) = google_utils.GetDriveSheets() 

    ### sanity check that you are in correct topic
    folder_name= DRIVE_TREE[config.folder_id]["title"]
    CheckCorrectTopic(config.topic,folder_name,config.folder_id,DRIVE_TREE)
    uat_round_number= GetUATRoundNumber(folder_name)
    
    ### get OOS pmids from summary sheet
    REJECTED={}
    for file_id in config.UAT_SUMMARY_NAMES[config.topic]:
        print(config.UAT_SUMMARY_NAMES[config.topic][file_id])
        if(not REJECTED and (( folder_name != "HerCan QC Previous Scoring Sheets" and config.UAT_SUMMARY_NAMES[config.topic][file_id] != "Hereditary Cancer QC Summary") or 
           (folder_name == "HerCan QC Previous Scoring Sheets" and config.UAT_SUMMARY_NAMES[config.topic][file_id] == "Hereditary Cancer QC Summary")) ):
            REJECTED=GetRejectedFromUatSummary( uat_round_number,file_id,config.UAT_SUMMARY_NAMES[config.topic][file_id] ,SHEETS)
        if(REJECTED):
            break
    print("HERE")
    print(json.dumps(REJECTED, sort_keys=True, indent=4))
   
    ### retrieve/init pmid list
    PMID_LISTS= utils.InitPmidData()
    ### Init PMID_LISTS for UAT and this topic
    PMID_LISTS= utils.InitTopicData(config.topic,PMID_LISTS);
    ### clean up if already have this folder data
    PMID_LISTS = CleanUpUAT(PMID_LISTS)

    ### simple log the list of files found for this folder
    utils.WarnMe("INFO", "Ready to extract UAT pmids in folder_name: '"+folder_name+"' folder_id: '"+ config.folder_id+ "'" )
    FS = utils.OpenFile(config.OUTPUTPATH+"/SHEETS_"+config.UATorQC+"_"+config.topic+"_"+folder_name+".txt", 'w',"utf8")
    for i in SHEETS.openall():
        if(i.id in (DRIVE_TREE[config.folder_id]["children"])):
           FS.write(i.title+"\t"+i.id+"\n")
    FS.close()

    ###  get data for each files in this folder
    nchild = 0
    nuat = 0
    for child_id in sorted(DRIVE_TREE[config.folder_id]["children"], key=lambda x:DRIVE_TREE[config.folder_id]["children"][x]["title"]):
        child_name = DRIVE_TREE[config.folder_id]["children"][child_id]["title"]
        nchild += 1
        if(DRIVE_TREE[config.folder_id]["children"][child_id]["type"] == "spreadsheet"):
            content = google_utils.GetSheetContent(child_name, child_id, "Publications", SHEETS)                 
            if content:
                nuat += 1
                print("folder name '" + folder_name+"' child_name:'"+str(child_name) + "' nchild:'" + str(nchild) + "' nuat " + str(nuat))
                PMID_LISTS = ExtractPmidsFromUatContent(child_name,child_id,config.topic,config.folder_id, content,PMID_LISTS,REJECTED)
        else:
            print("I am a folder:'" + child_name + "'")   
    ### store / output hash              
    utils.StoreHash(PMID_LISTS, config.pmid_list_store_path)
 
    FILE = utils.OpenFile(config.pmid_list_path+".json", 'w',"utf8")
    FILE.write(json.dumps(PMID_LISTS, sort_keys=True, indent=4))
    FILE.close()  
    return(PMID_LISTS,DRIVE_TREE)
#### END GetUATPmids



"""GetRejectedFromUatSummary"""
def GetRejectedFromUatSummary(uat_round_number,file_id,file_title , SHEETS):
    logging.debug("")  
    utils.WarnMe("info", "google sheet file_title: '"+file_title+"', uat_round_number: '"+str(uat_round_number)+"', file_id: '"+file_id+"'" )
    REJECTED={}
    sheet = SHEETS.open_by_key(file_id)
    worksheet_list = sheet.worksheets()
    for wk in worksheet_list:
        if(not REJECTED and re.search(r'Article Coverage Stats',wk.title, re.IGNORECASE)  ):
            utils.WarnMe("info", "google tab sheet title: '"+wk.title+"', file_title: '"+file_title+"', uat_round_number: '"+str(uat_round_number)+"', file_id: '"+file_id+"'" )
            content = google_utils.GetSheetContent(file_title, file_id, wk.title, SHEETS)      
            headerOOS=0  
            ncolOOS=0 
            npmid=0       
            current_round=0 
            nrow=0
            if content:
                for row in content:
                    ### get header
                    if(nrow==0):
                       (headerOOS, ncolOOS,ncolPMID)=GetOOSHeaderInfo(row,file_title, file_id)
                       utils.WarnMe("info", "Using  headerOOS: '"+headerOOS+"' column number '"+str(ncolOOS)+"' in google sheet file_title '"+file_title+"', file_id'"+file_id+"'" )
                    ### get rejected from old UAT summary format
                    elif(ncolOOS >0):
                        ### find 1st round number round
                        if(re.search(r'(.* )?round (?P<round>\d+)( .*)?', row[0], re.IGNORECASE)):
                            if(current_round==uat_round_number):
                                utils.WarnMe("warning", "I am done recovering rejected pmids for UAT round: '"+str(uat_round_number)+"' in '"+file_title+"', file_id'"+file_id+"'" )
                                #break
                            current_round=GetUATRoundNumber(row[0])
                        elif(nrow==1):
                             current_round=1
                        ### get rejected
                        if(current_round==uat_round_number and  row[ncolOOS] !="" ):
                            if( headerOOS == "OOS/FP pmids" ):
                                pmids_list=re.split(' ', row[ncolOOS])
                                for pmid in pmids_list:
                                    int_pmid=0
                                    try: 
                                        int_pmid =int(pmid)
                                    except Exception: 
                                        if(pmid != "" and pmid != "?" ):
                                            utils.WarnMe("CRITICAL", "At row: '"+str(nrow)+"' NOT A PMID: '"+pmid+ "'")  
                                            sys.exit()
                                    if(int_pmid>0):
                                        REJECTED[int_pmid]=1                         
                            ### get rejected from old UAT summary format
                            elif( headerOOS == "Content Error category" and ncolPMID>0):
                                pmid = row[ncolPMID]
                                int_pmid=0
                                try: 
                                    int_pmid =int(pmid)
                                except Exception: 
                                    if(pmid != "" and pmid != "?" ):
                                        utils.WarnMe("CRITICAL", "At row: '"+str(nrow)+"' NOT A PMID: '"+pmid+ "'")  
                                        sys.exit()
                                if(int_pmid>0):
                                    info=row[ncolOOS]
                                    print(pmid)
                                    print(info)
                                    int_info = utils.GetPmidsQCCode(nrow,pmid, info)
                                    if(int_info>0):
                                        REJECTED[pmid]=int_info
                            else:
                                utils.WarnMe("CRITICAL", "SHOULD NOT BE HERE - headerOOS: '"+headerOOS+"' column number '"+str(ncolOOS)+"' in google sheet file_title '"+file_title+"', file_id'"+file_id+"'" )
                                sys.exit()
                    else:
                        utils.WarnMe("warning", "I didn't find any headerOOS in google sheet file_title '"+file_title+"', file_id'"+file_id+"'" )
                        break
                    ### find if has "OOS/FP pmids" or "CQE Error category"
                    nrow+=1
                ### end loop on rows
                #print(content)
            ### end if content
            else:
                utils.WarnMe("critical", "I could not find an 'Article Coverage Stats' like tab in google sheet file_title '"+file_title+"', file_id'"+file_id+"'" )            
    return REJECTED
#### END GetRejectedFromUatSummary

"""GetOOSHeaderInfo"""
def GetOOSHeaderInfo(row,file_title, file_id,):
    logging.debug("")  
    header=row
    ncol=0
    headerOOS=0  
    ncolOOS=0
    ncolPMID=0
    for c in header:
        if(c=='PMID (1 PMID per row)'):
            ncolPMID=ncol
        if(c=="OOS/FP pmids" or c == "Content Error category"):
            headerOOS=c
            ncolOOS=ncol
            utils.WarnMe("info", "Using  headerOOS: '"+headerOOS+"' column number '"+str(ncolOOS)+"' in google sheet file_title '"+file_title+"', file_id'"+file_id+"'" )
            break
        ncol+=1
    return  (headerOOS, ncolOOS,ncolPMID)
#### END GetOOSHeaderInfo

"""GetUATRoundNumber"""
def GetUATRoundNumber(folder_name):
    logging.debug("folder_name: '"+folder_name+"'")  
    regex=re.compile(r'(.* )?round (?P<round>\d+)( .*)?', re.IGNORECASE)
    result = regex.match(folder_name)
    nround=0
    ### get round number
    if(result):
        nround=result.groupdict()["round"]
        utils.WarnMe("info", "The "+config.UATorQC+" folder_name: '"+folder_name+"' is for round number '"+ str(nround)+ "'" )
    else:
        nround=1
        utils.WarnMe("warning", "The "+config.UATorQC+" folder_name: '"+folder_name+"' has DEFAULT for round number '"+ str(nround)+ "'" )
    return int(nround)
#### END GetUATRoundNumber

"""CleanUpUAT"""
def CleanUpUAT(PMID_LISTS):
    logging.debug("")  
    PMID_LISTS[config.topic][config.UATorQC][config.folder_id]={}
    NEW_LIST= copy.deepcopy(PMID_LISTS)
    for pmid in PMID_LISTS["PMIDS"]:
        if(config.topic in PMID_LISTS["PMIDS"][pmid]):
           if( config.UATorQC in PMID_LISTS["PMIDS"][pmid][config.topic]):
                if config.folder_id in PMID_LISTS["PMIDS"][pmid][config.topic][config.UATorQC]:
                    del NEW_LIST["PMIDS"][pmid][config.topic][config.UATorQC][config.folder_id]  
                if(len(NEW_LIST["PMIDS"][pmid][config.topic][config.UATorQC]) ==0):
                    del NEW_LIST["PMIDS"][pmid][config.topic][config.UATorQC]
                if(len(NEW_LIST["PMIDS"][pmid][config.topic]) ==0)     :
                    del NEW_LIST["PMIDS"][pmid][config.topic]
                if(len(NEW_LIST["PMIDS"][pmid]) ==0)     :
                    del NEW_LIST["PMIDS"][pmid]     
    return NEW_LIST
#### END CleanUpUAT

"""CheckCorrectTopic"""
def CheckCorrectTopic(topic,folder_name,folder_id, DRIVE_TREE):
    logging.debug("folder_name: '"+ folder_name+ "' - folder_id: '"+folder_id+ "' - topic:'"+topic + "'")
    if(re.search(r''+config.UAT_TOPIC_NAMES[topic]+'',folder_name, re.IGNORECASE)):
        utils.WarnMe("info", "The "+config.UATorQC+" folder_name: '"+folder_name+"' is IN the correct topic: '"+ topic+ "'" )
    else:
        utils.WarnMe("info", "The "+config.UATorQC+" folder_name: '"+folder_name+"' is NOT the correct topic: '"+ topic+ "' - I will attempt to rescue using the parent folder" )
        if("parents" in DRIVE_TREE[folder_id] ):
            for key in DRIVE_TREE[folder_id]["parents"]:
                if key in DRIVE_TREE:
                    if(re.search(r''+config.UAT_TOPIC_NAMES[topic]+'', DRIVE_TREE[key]["title"], re.IGNORECASE)):
                        utils.WarnMe("info", "The "+config.UATorQC+" folder_name: '"+folder_name+"' in the parent folder '"+DRIVE_TREE[key]["title"]+"' is IN the correct topic: '"+ topic+ "'" )
                    else:
                        utils.WarnMe("CRITICAL", "The "+config.UATorQC+" folder_name: '"+folder_name+"' in the parent folder '"+DRIVE_TREE[key]["title"]+"' is NOT in the correct topic: '"+ topic+ "'" )
                        sys.exit()
                else:
                    utils.WarnMe("CRITICAL", "The "+config.UATorQC+" folder_name: '"+folder_name+"' has the parent folder_id: '"+key+"' not found in google drive tree and is NOT in the correct topic: '"+ topic+ "'" )
                    sys.exit()
        else:
            utils.WarnMe("CRITICAL", "The "+config.UATorQC+" folder_name: '"+folder_name+"' has the NO parent folder and is NOT in the correct topic: '"+ topic+ "'" )
            sys.exit()
    """ check that uat folder or QC file is for the correct topic """
#### END CheckCorrectTopic

"""ExtractPmidsFromUatContent"""
def ExtractPmidsFromUatContent(child_name,child_id,topic, folder, content,PMID_LISTS,REJECT):
    logging.debug("child_name: '"+ child_name+ "' - child_id: '"+child_id+ "' - topic:'"+topic + "'")
    """extract pmids from a UAT file content """
    nrow=0
    ncol=0
    for row in content:
        if(nrow==0):
            header=row
            for c in header:
                if(c=="Pubmed ID"):
                    break
                ncol+=1
        elif ncol != len(row) :
            pmid_info= row[ncol]
            pmids = re.findall(r'\d+[\.]?\d*',pmid_info)
            for pmid in pmids:
                try:
                    int(pmid)
                except Exception: 
                    utils.WarnMe("critical", "pmid: '"+pmid+"' is not a valid PMID in child_name: '"+ child_name+ "' - child_id: '"+child_id+ "' - topic:'"+topic + "'" )
                    FILE_ISSUES = utils.OpenFile(config.pmid_issues_path,"a","utf8")
                    FILE_ISSUES.write(topic+"\t"+pmid +"\t"+child_id+"\t"+child_name+"\n")
                    FILE_ISSUES.close()
                if(not pmid in REJECT):
                    if(not pmid in PMID_LISTS[topic]["UAT"][folder] ):
                        PMID_LISTS[topic]["UAT"][folder][pmid]={}
                    PMID_LISTS[topic]["UAT"][folder][pmid][child_id]=child_name
                    if(not pmid in PMID_LISTS["PMIDS"] ):
                        PMID_LISTS["PMIDS"][pmid]={}
                    if(not topic in PMID_LISTS["PMIDS"][pmid] ):
                        PMID_LISTS["PMIDS"][pmid][topic]={}
                    if(not "UAT" in PMID_LISTS["PMIDS"][pmid][topic] ):
                        PMID_LISTS["PMIDS"][pmid][topic]["UAT"]={}
                    if(not folder in PMID_LISTS["PMIDS"][pmid][topic]["UAT"] ):
                        PMID_LISTS["PMIDS"][pmid][topic]["UAT"][folder]={}
                    PMID_LISTS["PMIDS"][pmid][topic]["UAT"][folder][child_id]=child_name
        else:
            utils.WarnMe("critical", "could not find Pubmed ID in header for file child_name: '"+ child_name+ "' - child_id: '"+child_id+ "' - topic:'"+topic + "'")
        nrow+=1
    return PMID_LISTS
#### END ExtractPmidsFromUatContent