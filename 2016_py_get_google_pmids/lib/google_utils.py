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
import inspect
import re
from pydrive.auth import GoogleAuth 
from pydrive.drive import GoogleDrive  
import gspread
# from oauth2client.client import flow_from_clientsecrets
from oauth2client.file import Storage
from oauth2client.service_account import ServiceAccountCredentials

### get custom modules ###
sys.path.append(os.path.abspath(os.path.dirname(__file__)))  
import config
import utils

######################### GOOGLE SHEETS RELATED functions #########################
"""GetSheets"""
def GetDriveSheets():
    logging.debug("")
    """ get new credential ad reload the google drive and drive tree """
#     (drive) = GetGoogleDriveAuthorization()
#       
#     utils.WarnMe("info", "getting drive")
#     file_list = drive.ListFile({'q': "trashed=false"}).GetList()
# #     utils.StoreHash(file_list, config.full_drive_store_path)     
#     FILE = utils.OpenFile(config.full_drive_path, 'w',"utf8")
#     FILE.write(json.dumps(file_list , sort_keys=True, indent=4))
#     FILE.close()    
#       
#     """ organize drive hierarchy """
#     DRIVE_FILE_TREE = GetDriveTree(file_list)
#       
#     """ CheckParentsChildren completeness"""
#     DRIVE_FILE_TREE = CheckParentsChildren(DRIVE_FILE_TREE)
#           
#     FILE = utils.OpenFile(config.tree_drive_path, 'w',"utf8")
#     FILE.write(json.dumps(DRIVE_FILE_TREE , sort_keys=True, indent=4))
#     FILE.close()
#           
#     utils.StoreHash(DRIVE_FILE_TREE, config.OUTPUTPATH + "/TREE.gz")
    DRIVE_FILE_TREE = utils.RetrieveHash(config.OUTPUTPATH + "/TREE.gz")
    """ loop on topics """
    (SHEETS) = getGoogleSheetsAuthorization()
    return(SHEETS,DRIVE_FILE_TREE) 
#### END GetSheets

# This module checks for and creates the credential file from Google
def getGoogleSheetsAuthorization():
    logging.debug("credential file name:'" + config.credentials_sheet_filename + "'")
    """ get google sheet credential from json keyfile """
    
#     scope = ['https://spreadsheets.google.com/feeds']
#  
#     credentials = ServiceAccountCredentials.from_json_keyfile_name('pydrive-ec7004f48761.json', scope)
#     return(gspread.authorize(credentials))
    
    """ find existing credential or get drive client keys from settings.yalm """
    #### First check to see if a credential file is present in the current folder, if so use it
    if (os.path.isfile(config.credentials_sheet_filename)): 
        logging.info("Credential file (" + config.credentials_sheet_filename + ") found.")
        storage = Storage(config.credentials_sheet_filename)
        credentials = storage.get() 
        print(credentials)
# 
#     else:
#         # If it's not present, get access
#         logging.info("Credential file (" + config.credentials_sheet_filename + ") not found. Generating a new one.\n")
#         # Define the flow objectb
#         flow = flow_from_clientsecrets('client_secrets.json',
#             scope='https://spreadsheets.google.com/feeds/',
#             redirect_uri='urn:ietf:wg:oauth:2.0:oob')
#         # Get the authorization uri
#         auth_uri = flow.step1_get_authorize_url()
# 
#         # Prompt the user to follow the link and paste the code
#         sys.stdout.write("REQUEST: Open the following URL and paste the returned value: \n" + auth_uri + "\n")
#         sys.stdout.write("REQUEST: Paste the generated value here: ")
#         authorizationCode = input() #raw_input()
#         credentials = flow.step2_exchange(authorizationCode)
# 
#         # Save the HTTP object
#         sys.stdout.write("REQUEST: Storing the HTTP credential object in the current directory.\n")
#         storage = Storage(config.credentials_sheet_filename)
#         storage.put(credentials)
# 
#         # Retrieve the credential object. This will allow access to the required Google Docs
#         logging.info("Retrieving the credential object.\n")
#         credentials = storage.get()

    utils.WarnMe("info", "Done using/saving GOOGLE SHEET credentials to '" + config.credentials_sheet_filename + "' and accessing google drive")
    # Return the credential
    return(gspread.authorize(credentials));
    #### End

"""
GetSheetContent
"""
def GetSheetContent(child_name, child_id, tab_name, SHEETS):
    logging.debug("child_name:'" + child_name + "' - child_id:'" + child_id + "' - :'" + tab_name + "'")
    """ Get content from a google sheet """
#     utils.WarnMe("info", "child_name:'" + child_name + "' - child_id:'" + child_id + "' - :'" + tab_name + "'")
#     sheet = SHEETS.open_by_key("1sYUz1SwH9UV62odA74Ohb7xE5KSg3Ro-GozmwqhN0BQ")#child_id)
#     sheet = SHEETS.open(child_name)
    
    sheet = SHEETS.open_by_key(child_id)

    content = []
    if sheet:
        logging.info("Spreadsheet: '" + child_name + "' found and opened")
        tab = {}
        try:
            tab = sheet.worksheet(tab_name)
        except Exception:
            utils.WarnMe("warning", "\tWorksheet: '" + tab_name + "' was not found in the supplied Google Spreadsheet'" + child_name + "'. Please check the worksheet name.")
            # raise Exception ( "Worksheet: '" +tab_name +"' was not found in the supplied Google Spreadsheet: '" +child_name +"'. Please check the worksheet name.")
        if tab:
            logging.info( "Worksheet: '" + tab_name + "' found and opened")
            content = tab.get_all_values()
        else:
            utils.WarnMe("error", "Worksheet: '" + tab_name + "' NOT found in Spreadsheet: '" + child_name )
    else:
        utils.WarnMe("info", "Spreadsheet: '" + child_name + "' NOT found")
     
#     print("here"+json.dumps(content,  sort_keys=True, indent=4)) 
    
    return content
######################### END GOOGLE SHEETS RELATED functions #########################


######################### START FUNCTION #########################
######################### GOOGLE DRIVE RELATED functions #########################
"""
GetGoogleDriveAuthorization
"""
def GetGoogleDriveAuthorization():
    logging.debug("credential file name:'" + config.credentials_drive_filename + "'")
    """ find existing credential or get drive client keys from settings.yalm """

    gauth = GoogleAuth()

    # Try to load saved client credentials
    gauth.LoadCredentialsFile(config.credentials_drive_filename)
    if gauth.credentials is None:
        utils.WarnMe("info", "Could not find '" + config.credentials_drive_filename + "' - will ask to get token from url")

        # Authenticate if they're not there
        gauth = GoogleAuth()
#         gauth.oauth_scope(['https://www.googleapis.com/auth/drive','https://spreadsheets.google.com/feeds'])
        auth_url = gauth.GetAuthUrl()  # Create authentication url user needs to visit
        # Prompt the user to follow the link and paste the code
        sys.stdout.write("REQUEST: Open the following URL and paste the returned value: \n" + auth_url + "\n")
        sys.stdout.write("REQUEST: Paste the generated value here: ")
        code = input()  # raw_input()
        gauth.Auth(code)
         
        gauth.LocalWebserverAuth()
        utils.WarnMe("info", "Done asking token from url")

    elif gauth.access_token_expired:
        # Refresh them if expired
        utils.WarnMe("info", "Refreshing expired token")
        gauth.Refresh()
    else:
        # Initialize the saved creds
        utils.WarnMe("info", "Initialize the saved credentials")
        gauth.Authorize()
    # Save the current credentials to a file
    gauth.SaveCredentialsFile(config.credentials_drive_filename)
    drive = GoogleDrive(gauth)
    utils.WarnMe("info", "Done using/saving GOOGLE DRIVE credentials to '" + config.credentials_drive_filename + "' and accessing google drive")
    return (drive)


######################### GetDriveTree functions #########################
""" GetDriveTree(FILE_LIST) """
def GetDriveTree(FILE_LIST):
    """ constructs the tree of folder/children/files from my google drive """
    logging.debug("")
    DRIVE_TREE = {}
    for file_info in FILE_LIST:
        """ init file/folder info """
        filename = file_info['title']
        fileid = file_info['id']
        #### init
        if not fileid in DRIVE_TREE:
            DRIVE_TREE[fileid] = {}       
        DRIVE_TREE = GetFileFolderInfo(filename, fileid, DRIVE_TREE, file_info)    
        #### clean up
        del file_info, filename, fileid
    return DRIVE_TREE
#### END GetFileFolderInfo

""" GetFileFolderInfo """
def GetFileFolderInfo(filename, fileid, DRIVE_TREE, INFO):
    """ recover the file/folder information """
    logging.debug("filename:'" + filename + "' - fileid:'" + fileid + "'")
    for key in ('title', "id", 'alternateLink',"modifiedDate", "parents"):
        if key == "parents" and INFO[key]:
            """ assign parents and children names and ids """
            DRIVE_TREE = GetParentChildrenInfo(filename, fileid, DRIVE_TREE, INFO[key])
            
        elif key == "alternateLink" :
            """ Assign type - folder vs spreadsheet """   
            DRIVE_TREE = GetFolderFileType(filename, fileid, DRIVE_TREE, INFO[key])
            if not fileid in DRIVE_TREE: 
                break
        elif key != "parents":
            DRIVE_TREE[fileid][key] = INFO[key]
        #### clean up
        del key
    return DRIVE_TREE
#### END GetFileFolderInfo


###### GetFileFolderInfo functions #######
""" GetFolderFileType """
def GetFolderFileType(filename, fileid, DRIVE_TREE, link):
    """ assign the Folder/File Type information """
    logging.debug("filename:'" + filename + "' - fileid:'" + fileid + "'")  
    
    DRIVE_TREE[fileid]["alternateLink"] = link
    if re.search("folderview", link):
        DRIVE_TREE[fileid]["type"] = "folder"
    elif re.search("spreadsheets", link):
        DRIVE_TREE[fileid]["type"] = "spreadsheet"
    else:
        """ ignore non folder/ sheet files """
        DRIVE_TREE.pop(fileid)
    return DRIVE_TREE
#### END GetFolderFileType

""" GetParentChildrenInfo """
def GetParentChildrenInfo(child_name, child_id, DRIVE_TREE, PARENTS):
    """ recover the parent/children information """
    logging.debug("child_name:'" + child_name + "' - child_id:'" + child_id + "'")    
   
    for PARENT in PARENTS:        
        parent_id = PARENT["id"]
        DRIVE_TREE[child_id]["parents"] = {}  # init
        DRIVE_TREE[child_id]["parents"][parent_id] = {}
        DRIVE_TREE[child_id]["parents"][parent_id]["id"] = parent_id
        
        """ update children info in parent """
        if parent_id in DRIVE_TREE:
            parent_name = DRIVE_TREE[parent_id]["title"] 
            DRIVE_TREE[child_id]["parents"][parent_id]["title"] = parent_name
            
            #### init children info
            if not "children" in DRIVE_TREE[parent_id]:
                DRIVE_TREE[parent_id]["children"] = {}
            DRIVE_TREE[parent_id]["children"][child_id] = {"title":child_name, "id":child_id, "type" :DRIVE_TREE[child_id]["type"],"modifiedDate" :DRIVE_TREE[child_id]["modifiedDate"] }
            #### clean up
            del parent_name
        #### clean up
        del PARENT, parent_id
    return DRIVE_TREE
#### END GetParentChildrenInfo
###### END GetFileFolderInfo functions ######
######################### END GetDriveTree functions ######################### 


#########################  CheckParentsChildren functions ######################### 
""" CheckParentsChildren """
def CheckParentsChildren(DRIVE_TREE):
    """ check have all the parent/children information """
    logging.debug("")    

    for fileid in sorted(DRIVE_TREE):
        FILEINFO = DRIVE_TREE[fileid]
        filename = FILEINFO["title"]
        if "parents" in FILEINFO:
            (FILEINFO, DRIVE_TREE) = CheckParentsInfo(filename, fileid, FILEINFO, DRIVE_TREE)
        if "children" in FILEINFO:
            (FILEINFO, DRIVE_TREE) = CheckChildrenInfo(filename, fileid, FILEINFO, DRIVE_TREE)
        #### clean up          
        del filename, fileid, FILEINFO

    return DRIVE_TREE
#### END CheckParentsChildren

""" CheckChildrenInfo """
def CheckChildrenInfo(parent_name, parent_id, FILEINFO, DRIVE_TREE):
    """ check have all the children information """
    logging.debug("parent_name:'" + parent_name + "' - parent_id:'" + parent_id + "'") 
    
    for child_id in sorted(FILEINFO["children"]):
        child_name = FILEINFO["children"][child_id]["title"]
        """ check children already has parent name and info"""
        if "parents" in DRIVE_TREE[child_id]:
            if parent_id in DRIVE_TREE[child_id]["parents"]:
                if not "title" in DRIVE_TREE[child_id]["parents"]:
                    DRIVE_TREE[child_id]["parents"][parent_id]["title"] = parent_name
            else:
                utils.WarnMe("CRITICAL", "no parent ID found in children file '" + child_name + "' - type:'" + FILEINFO["type"] + "' - parent_name:'" + parent_name + "'")
#                 print(json.dumps(FILEINFO, sort_keys=True, indent=4))
#                 print(json.dumps(DRIVE_TREE[child_id], sort_keys=True, indent=4))
#                 sys.exit()
        else:
            utils.WarnMe("CRITICAL", "no parent DEFINED found in children file '" + child_name + "' - type:'" + FILEINFO["type"] + "' - parent_name:'" + parent_name + "'")
            sys.exit()
        #### clean up
        del child_name     
    return (FILEINFO, DRIVE_TREE) 
#### END CheckParentsInfo

""" CheckParentsInfo """
def CheckParentsInfo(child_name, child_id, FILEINFO, DRIVE_TREE):
    """ check have all the parent information """
    logging.debug("child_name:'" + child_name + "' - child_id:'" + child_id + "'") 

    for parent_id in sorted(FILEINFO["parents"]):
        """ check already have parent name and info"""
        parent_name = ""
        if "title" in FILEINFO["parents"][parent_id]:
            parent_name = FILEINFO["parents"][parent_id]["title"]
        elif parent_id in DRIVE_TREE:
            logging.info("Had to update parent info for child_name:'" + child_name + "' - type:'" + FILEINFO["type"] + "' - parent_name:'" + parent_name + "' - parent_id:'" + parent_id + "'")
            parent_name = DRIVE_TREE[parent_id]["title"]
            FILEINFO["parents"][parent_id]["title"] = parent_name
        elif FILEINFO["type"] == "spreadsheet":
            utils.WarnMe("info", "no parent info found for child_name:'" + child_name + "' - type:'" + FILEINFO["type"] + "' - parent_id:'" + parent_id + "'")

        """ check children in parent """
        if parent_name != "":
            if "children" in DRIVE_TREE[parent_id] and child_id in DRIVE_TREE[parent_id]["children"]:
                pass
            else:
                if not "children" in DRIVE_TREE[parent_id] :
                    DRIVE_TREE[parent_id]["children"] = {}
                DRIVE_TREE[parent_id]["children"][child_id] = {"title":child_name, "id":child_id, "type" :FILEINFO["type"],'modifiedDate':FILEINFO["modifiedDate"] }
        #### clean up
        del parent_name, parent_id     
    return (FILEINFO, DRIVE_TREE) 
#### END CheckParentsInfo
######################### END CheckParentsChildren functions ######################### 
######################### END GOOGLE DRIVE RELATED functions #########################
