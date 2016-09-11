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
import os.path
import re
import copy
import pytz
from pytz import timezone
from datetime import datetime  # @UnresolvedImport
### get custom modules ###
sys.path.append( os.path.abspath(os.path.dirname(__file__)))  
import config  
import utils 
import google_utils  


##################### START FUNCTIONS ####################
"""GetQCPmids"""
def GetQCPmids():
    logging.debug("")

    GetLatestDate()

    """extract pmids from QC files if was updated since last run """
    PMID_LISTS = utils.InitPmidData()

    ### get folder file list      
    (SHEETS,DRIVE_TREE) = google_utils.GetDriveSheets() 
    
    folder_name= DRIVE_TREE[config.folder_id]["title"]
    ### simple log the list of files found for this folder
    utils.WarnMe("INFO", "Ready to extract QC pmids in folder_name: '"+folder_name+"' folder_id: '"+ config.folder_id+ "'" )
    FS = utils.OpenFile(config.OUTPUTPATH+"/SHEETS_"+config.UATorQC+"_"+folder_name+".txt", 'w',"utf8")
    for i in SHEETS.openall():
        if(i.id in (DRIVE_TREE[config.folder_id]["children"])):
           FS.write(i.title+"\t"+i.id+"\n")
    FS.close()
    
    regex=re.compile(r'Clinical QC-coverage-(?P<version>v\d\d\d)-(?P<topic>'+config.QC_REGEX+')?(-\D+)?-(?P<size>\d+|\d+-updated)$', re.IGNORECASE)
    nchild = 0
    nuat = 0
    for child_id in sorted(DRIVE_TREE[config.folder_id]["children"], key=lambda x:DRIVE_TREE[config.folder_id]["children"][x]["title"]):
        child_name = DRIVE_TREE[config.folder_id]["children"][child_id]["title"]
        nchild += 1
        if(DRIVE_TREE[config.folder_id]["children"][child_id]["type"] == "spreadsheet"):
            result = regex.match(child_name)
            ### get topic for this file
            print(result)
            if(result):
                (topic,version) = GetQcTopic(result)
                need_to_run=1
                if (config.version=="" and config.topic==""):
                    ### get topic for this file
                    need_to_run= CheckQcBatchDate(config.qc_file_date,DRIVE_TREE[config.folder_id]["children"][child_id]["modifiedDate"])
                elif (config.version==version and config.topic==topic):
                    need_to_run= -1.0;
                elif(config.version==version and config.topic!=topic):
                    utils.WarnMe("CRITICAL", " Found correct version:'"+version+"' - but wrong topic. I found '"+topic+"' when I expected:'"+config.topic+"' child_name:'"+str(child_name) + "' nchild:'" + str(nchild)+"'" )
                    sys.exit()
                ### extract QC data
                if(need_to_run <0.0):
                    utils.WarnMe("info", " Updating version:'"+version+"', topic:'"+topic+"' child_name:'"+str(child_name) + "' nchild:'" + str(nchild)+"'" )
                 
                    ### INIT
                    PMID_LISTS=utils.InitTopicData(topic,PMID_LISTS)
                    ### clean up if already have this folder data
                    PMID_LISTS=CleanUpQC(topic,child_id,child_name,PMID_LISTS)
            
                    ### get content
                    report =  GetContentOrReportSheet(config.QC_REPORT_TAB,child_name,child_id,version,topic,config.folder_id, SHEETS)
                    content = GetContentOrReportSheet(config.QC_CONTENT_TAB,child_name,child_id,version,topic,config.folder_id, SHEETS)
           
                    if(report and content):
                        PMID_LISTS = ExtractPmidsFromQcContent(child_name,child_id,version,topic,config.folder_id,report, content,PMID_LISTS)
                    else:
                        utils.WarnMe("CRITICAL", "could not get the content and report from sheet child_name:'"+str(child_name) + "' nchild:'" + str(nchild)+"'" )
                        sys.exit()
                elif(config.version=="" and config.topic==""):
                    utils.WarnMe("info", " No need to update version:'"+version+"', topic:'"+topic+"' child_name:'"+str(child_name) + "' nchild:'" + str(nchild)+"'" )
                if (config.version==version and config.topic==topic):
                    utils.WarnMe("info", " DONE with version:'"+version+"', topic:'"+topic+"' child_name:'"+str(child_name) + "' nchild:'" + str(nchild)+"'" )
                    break

            elif(child_name != "Clinical QC-coverage-v026-HCv2-50-deprecate" 
                   and child_name != "Clinical QC-coverage-vXXX-template-50"
                   and child_name != "Clinical QC-coverage-v024-heredv2-50-abandoned"):
                utils.WarnMe("info", "File names '"+child_name+"' doesn't match QC file name" )
                sys.exit()
        else:
            print("I am a folder:'" + child_name + "'")
                                    
    utils.StoreHash(PMID_LISTS, config.pmid_list_store_path)
 
    FILE = utils.OpenFile(config.pmid_list_path+".json", 'w',"utf8")
    FILE.write(json.dumps(PMID_LISTS, sort_keys=True, indent=4))
    FILE.close()  
    
    ### record latest date 
    FILE = utils.OpenFile(config.qcdate_file_path, 'w',"utf8")# TODO change date to datetime.now(timezone('US/Pacific'))
#     FILE.write(config.qc_file_date)
    FILE.write(datetime.now(timezone('US/Pacific')).strftime('%Y-%m-%dT%H:%M:%S.%f%Z'))
    FILE.close() 
    return(PMID_LISTS,DRIVE_TREE)
#### END GetQCPmids

def GetLatestDate():
    logging.debug("")
    if (os.path.isfile(config.qcdate_file_path)):
        FILE = utils.OpenFile(config.qcdate_file_path, 'r',"utf8")
        config.qc_file_date=FILE.read()
        FILE.close()  
        utils.WarnMe("info", "I found the file '"+config.qcdate_file_path+"', the date was '"+config.qc_file_date+"'" )
    else:
        config.qc_file_date="1900-01-01T00:00:44.522872PDT" 
#### END GetLatestDate

"""GetQcTopic"""
def GetQcTopic(regex_result):
    re_topic=regex_result.groupdict()["topic"]
    re_version=regex_result.groupdict()["version"]
    logging.debug("topic: '"+re_topic+"' version: '"+re_version+"'")
    """GetQcTopic: extracts the appropriate topics for a QC google sheet"""
    found_topic=""
    print(re_topic)
    print(re_version)
    
    for topic_reg in config.QC_TOPIC_NAMES:
        regex=topic_reg["regex"]
        if (re.search(regex,re_topic,re.IGNORECASE)):
            found_topic=topic_reg["topic"]
            break
    return (found_topic,re_version)
#### END GetQcTopic

def CheckQcBatchDate(latest_date, batch_date):
    logging.debug("latest_date: '"+latest_date+"', batch_date: '"+batch_date+"'")

    ### convert run time to standart time
    local = pytz.timezone ("America/Los_Angeles")
    naive = datetime.strptime(config.qc_file_date, '%Y-%m-%dT%H:%M:%S.%f%Z')
    local_dt = local.localize(naive, is_dst=None)
    last = local_dt.astimezone (pytz.utc)
    last=last.strftime ("%Y-%m-%dT%H:%M:%S.%f%Z")
    
    ### add zone to batch date
    batch_date=re.sub(r'Z$',"000UTC",batch_date)

    d1 = datetime.strptime(last, '%Y-%m-%dT%H:%M:%S.%f%Z')
    d2 = datetime.strptime(batch_date, '%Y-%m-%dT%H:%M:%S.%f%Z')

    delta=(d1-d2).total_seconds()
    utils.WarnMe("info", "latest_date: '"+last+"', batch_date: '"+batch_date+"', d1-d2: '"+str(d1-d2)+"', delta: '"+ str(delta)+"'" )
    return(delta)
#### END CheckQcBatchDate

"""CleanUpQC"""
def CleanUpQC(topic,child_id,child_name,PMID_LISTS):
    logging.debug("topic: '"+topic+"',child_id: '"+child_id+"'child_name: '"+child_name+"'")  
    ### ini folder data
    if not config.UATorQC in PMID_LISTS[topic]:
        PMID_LISTS[topic][config.UATorQC]={}
    if not config.folder_id in PMID_LISTS[topic][config.UATorQC]:
        PMID_LISTS[topic][config.UATorQC][config.folder_id]={}
    
    ### clean up topic pmid data for this file
    NEW_LIST= copy.deepcopy(PMID_LISTS)   
    for pmid in PMID_LISTS[topic][config.UATorQC][config.folder_id]:
        if child_id in PMID_LISTS[topic][config.UATorQC][config.folder_id][pmid]:
            del NEW_LIST[topic][config.UATorQC][config.folder_id][pmid][child_id]
        if(len(NEW_LIST[topic][config.UATorQC][config.folder_id][pmid]) == 0):
            del NEW_LIST[topic][config.UATorQC][config.folder_id][pmid]
    for pmid in PMID_LISTS["PMIDS"]:
        if(topic in PMID_LISTS["PMIDS"][pmid]):
           if( config.UATorQC in PMID_LISTS["PMIDS"][pmid][topic]):
                if config.folder_id in PMID_LISTS["PMIDS"][pmid][topic][config.UATorQC]:
                    if child_id in PMID_LISTS["PMIDS"][pmid][topic][config.UATorQC][config.folder_id]:
                        del NEW_LIST["PMIDS"][pmid][topic][config.UATorQC][config.folder_id][child_id]  
                        if(len(NEW_LIST["PMIDS"][pmid][topic][config.UATorQC][config.folder_id]) ==0):
                            del NEW_LIST["PMIDS"][pmid][topic][config.UATorQC][config.folder_id]
                        if(len(NEW_LIST["PMIDS"][pmid][topic][config.UATorQC]) ==0):
                            del NEW_LIST["PMIDS"][pmid][topic][config.UATorQC]
                        if(len(NEW_LIST["PMIDS"][pmid][topic]) ==0)     :
                            del NEW_LIST["PMIDS"][pmid][topic]
                        if(len(NEW_LIST["PMIDS"][pmid]) ==0)     :
                            del NEW_LIST["PMIDS"][pmid]     
    return NEW_LIST
#### END CleanUpQC

"""GetContentOrReportSheet"""
def GetContentOrReportSheet(TABS,child_name,child_id,version,topic,folder_id,SHEETS):
    logging.debug("child_name: '"+child_name+"' child_id: '"+child_id+"' topic: '"+topic+"' version: '"+version+"' folder_id: '"+folder_id+"'")
#     print("child_name: '"+child_name+"' child_id: '"+child_id+"' topic: '"+topic+"' version: '"+version+"' folder_id: '"+folder_id+"'")
    """ try to open report tab or alternative name"""
    for tab_name in TABS:
        try:
            content = google_utils.GetSheetContent(child_name, child_id, tab_name, SHEETS) 
        except Exception:
            utils.WarnMe("warning", "Worksheet: '" + tab_name + "' was not found in the supplied Google Spreadsheet'" + child_name + "'. Please check the worksheet name.")
            # raise Exception ( "Worksheet: '" +tab_name +"' was not found in the supplied Google Spreadsheet: '" +child_name +"'. Please check the worksheet name.")
        if content:
           utils.WarnMe("INFO", "Opened Worksheet: '" + tab_name + "' in Google Spreadsheet'" + child_name + "'.") 
           return content 
#### END GetContentOrReportSheet  

"""ExtractPmidsFromQcContent"""
def ExtractPmidsFromQcContent(child_name,child_id,version,topic,folder,report, content,PMID_LISTS):
    logging.debug("child_name: '"+child_name+"' child_id: '"+child_id+"' topic: '"+topic+"' version: '"+version+"' folder: '"+folder+"'")
    """ExtractPmidsFromQcContent: extracts pmids from report, and remove false positives found in content QC google sheet"""
    REPORT = GetQcReportPmids(child_name,child_id,version,topic,folder,report )
    CONTENT = GetQcContentPmids(child_name,child_id,version,topic,folder,content )
    
    RFILE = utils.OpenFile(config. OUTPUTPATH + "/content_" +child_name+".txt" , 'w',"utf8")
    for pmid in CONTENT:
        RFILE.write(pmid+"\n")
    RFILE.close() 
    
    RFILE = utils.OpenFile(config. OUTPUTPATH + "/report_" +child_name+".txt" , 'w',"utf8")
    FILE = utils.OpenFile(config. OUTPUTPATH + "/final_" +child_name+".txt" , 'w',"utf8")
    for pmid in REPORT:
        RFILE.write(pmid+"\n")
        if(not pmid in CONTENT):
            FILE.write(pmid+"\n")
            if(not pmid in PMID_LISTS[topic]["QC"][folder] ):
                PMID_LISTS[topic]["QC"][folder][pmid]={}
            PMID_LISTS[topic]["QC"][folder][pmid][child_id]=child_name
            if(not pmid in PMID_LISTS["PMIDS"] ):
                PMID_LISTS["PMIDS"][pmid]={}
            if(not topic in PMID_LISTS["PMIDS"][pmid] ):
                PMID_LISTS["PMIDS"][pmid][topic]={}
            if(not "QC" in PMID_LISTS["PMIDS"][pmid][topic] ):
                PMID_LISTS["PMIDS"][pmid][topic]["QC"]={}
            if(not folder in PMID_LISTS["PMIDS"][pmid][topic]["QC"] ):
                PMID_LISTS["PMIDS"][pmid][topic]["QC"][folder]={}
            PMID_LISTS["PMIDS"][pmid][topic]["QC"][folder][child_id]=child_name 
    RFILE.close() 
    FILE.close() 
    return PMID_LISTS
#### END ExtractPmidsFromQcContent  

"""GetQcContentPmids"""
def GetQcContentPmids(child_name,child_id,version,topic,folder_id,content ):
    logging.debug("child_name: '"+child_name+"' child_id: '"+child_id+"' topic: '"+topic+"' version: '"+version+"' folder_id: '"+folder_id+"'")
    """ extracts pmids from report tab in QC google sheet"""
    
    nrow=0
    nhead=-1
    npmid=-1
    ninfo=-1
    PMIDS={}
    for row in content:
        ### find header
        if(nrow == 0 and nhead<0):
            ncol=0
            for c in row:
                if(ncol in config.QC_PMID_COL and c in config.QC_PMID_COL[ncol]):
                    utils.WarnMe("Info", "Content pmid header of child_name: '"+child_name+ " at nrow:'" + str(nrow) + "' ncol'" + str(ncol) +  "' name'" + c + "'")          
                    npmid=ncol
                elif(ncol in config.QC_FALSE_POSITIVE_COL and c in config.QC_FALSE_POSITIVE_COL[ncol]):
                    utils.WarnMe("Info", "Content false positive header of child_name: '"+child_name+ " at nrow:'" + str(nrow) + "' ncol'" + str(ncol) + "' name'" + c + "'")          
                    ninfo=ncol
                    break
                ncol+=1
        ### recover pmids
        elif(ninfo>=0 and npmid>=0 ):
            pmid=row[npmid]
            info=row[ninfo]
            int_info = utils.GetPmidsQCCode(nrow,pmid, info)
            if(int_info>0):
                PMIDS[pmid]=int_info
        elif((ninfo>=0 and npmid<0 ) or(ninfo<0 and npmid>0 )):
            utils.WarnMe("CRITICAL", "Missing a header column: npmid='"+str(npmid)+"' ninfo='"+str(ninfo)+ "'")  
            sys.exit()
        nrow+=1
    print(PMIDS)
    return PMIDS
#### END GetQcContentPmids

"""GetQcReportPmids"""
def GetQcReportPmids(child_name,child_id,version,topic,folder_id,content ):
    logging.debug("child_name: '"+child_name+"' child_id: '"+child_id+"' topic: '"+topic+"' version: '"+version+"' folder_id: '"+folder_id+"'")
    """ extracts pmids from report tab in QC google sheet"""
    pattern_report=re.compile(', |,')
    nrow=0
    nhead=-1
    ncol=0
    PMIDS={}
    for row in content:
        ### find header
        if(nrow in config.QC_REPORT_COL and nhead<0):
            ncol=0
            for c in row:
                if(ncol == config.QC_REPORT_COL[nrow]["ncol"] and c in config.QC_REPORT_COL[nrow]["header"]):
                    utils.WarnMe("Info", "Report header of child_name: '"+child_name+ " at nrow:'" + str(nrow) + "' ncol'" + str(ncol) +  "' name'" + c + "'")          
                    nhead=nrow
                    break
                ncol+=1
        ### recover pmids
        elif(nhead>=0 and ncol == config.QC_REPORT_COL[nhead]["ncol"] and c in config.QC_REPORT_COL[nhead]["header"]):
            pmids_phrase = row[ncol]
            pmids_list=pattern_report.split(pmids_phrase)
            int_pmid=0;
            for pmid in pmids_list:
                try: 
                    int_pmid=int(pmid)
                except Exception: 
                    if(pmid != "" and pmid !="sur1 deletion TTC Phe 1388"):
                        utils.WarnMe("CRITICAL", "At row: '"+str(nrow)+"' NOT A PMID: '"+pmid+ "'")  
                        sys.exit()
                if(int_pmid>0):
                    PMIDS[pmid]=1
        nrow+=1
    return PMIDS
#### END GetQcReportPmids