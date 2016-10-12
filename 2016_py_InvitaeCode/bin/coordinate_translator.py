### @package coordinate_translator
# Translates a (0-based) transcript coordinate to a (0-based) genome coordinate.\n
# Assumes that the file with transcript information is small enough to keep information in a tuple.\n
# The transcripts names are assumed case-sensitive.\n
# <a href="InvitaeBioinformaticsExerciseB.pdf" target="_blank"><b>Invitae exercise</b></a>
###
# \mainpage CoordinateTranslator
# Translates a (0-based) transcript coordinate to a (0-based) genome coordinate.\n
# Assumes that the file with transcript information is small enough to keep information in a tuple.\n
# The transcripts names are assumed case-sensitive.\n
# See Invitae exercise <a href="InvitaeBioinformaticsExerciseB.pdf" target="_blank"><b>here</b></a>\n
#<PRE> 
# USAGE
#     coordinate_translator.py [OPTIONS] 
# DESCRIPTION
#     Translates (0-based) transcript coordinates to (0-based) genome coordinates
# OPTIONS
#         -t FILENAME, --transcripts FILENAME
#             FILENAME specifies the name of a four column (tab-separated) file containing the transcripts.
#             The first column is the transcript name, and the remaining three columns indicate it's genomic mapping: 
#             chromosome name, 0-based starting position on the chromosome, and CIGAR string indicating the mapping. 
#             By default, FILENAME='transcripts.txt'
#         -q FILENAME, --queries FILENAME
#             FILENAME specifies the name of a two column (tab-separated) file indicating a set of queries.
#             The first column is a transcript name, and the second column is a 0-based transcript coordinate.
#             By default, FILENAME='queries.txt'
#         -i PATH, --indir PATH
#             PATH specifies the full path to the input files. By default, PATH=<PATH2PACKAGE>/input_files
#         -f FILENAME, --outfile FILENAME
#             FILENAME specifies the name of a four column tab separated file with one row for each of the input queries.
#             The first two columns are exactly the two columns from the query file.
#             The remaining two columns are the chromosome name and chromosome coordinate, respectively.
#             By default, FILENAME='translation.txt'
#         -o PATH, --outdir PATH
#             PATH specifies the full path to the output and log files.
#             By default, PATH=<PATH2PACKAGE/output_files>
#         -r FILENAME, --logfile FILENAME
#             FILENAME specifies the name of the log file.
#             By default, FILENAME=logfile.log
#         -l LOGLEVEL, --loglevel LOGLEVEL
#             LOGLEVEL specifies the logging level (i.e. controls the verbosity) in the logging file.
#             LOGLEVEL='DEBUG', 'INFO', 'WARNING', 'ERROR' or 'CRITICAL' (by default LOGLEVEL=DEBUG)
#         -v VERBOSITY, --verbosity VERBOSITY
#             VERBOSITY= 0, 1 or 2 (by default VERBOSITY=2)
#             Default is VERBOSITY = 2: STDOUT format like log file: '%m/%d/%Y %H:%M:%S' - <file>.py:<class>.<function>:<line> - + <level> - <info>.
#             VERBOSITY = 0: STDOUT format (no time stamp): <file>.py:<class>.<function>:<line> - + <level> - <info>.
#             If VERBOSITY = None: no STDOUT.
# OTHER OPTIONS
#         -h, --help
#             Output this help manual
""" 
@file coordinate_translator.py
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
import os.path
import inspect
import re
import copy

#---# custom modules 
sys.path.append(re.sub("bin", "", os.path.abspath(os.path.dirname(__file__))) + "lib")  
import config
from utility_functions import UtilityFunctions

### Translates (0-based) transcript coordinates to (0-based) genome coordinates\n 
# USAGE: coordinate_translator.py [OPTIONS] \n
### use option -h or --help to output the manual
class CoordinateTranslator():   
    
    def __init__(self):
        """ CoordinateTranslator Class
        Attributes:
            _utils (Utils): Set if utility functions.
        """
        self._utils = UtilityFunctions();
        self.mapping_info = {"transcript":{}, "genome":{}}
        self._utils.warn_me("INFO", "Starting CoordinateTranslator")
        self._cigar_pat = re.compile(r"\d+[MIDNSHP=X]{1}")
        self._cigar_check = re.compile(r"[a-z]|[A-C]|[E-G]|[J-L]|O|[Q-R]|[T-W]|[Y-Z]|(\d+[MIDNSHP=X]{2,10})")
    #---# END __init__ 
  
    ### Extract transcripts information
    def _open_transcript_file(self):
        """ Read the file with transcript's information
            record transcripts information in a Dictionary
            -> Assumes memory is not an issue
        """
        self._utils.warn_me("INFO", " opening transcript file " + config.TRANSCRIPT_FILE_PATH)
        F = self._utils.open_file(config.TRANSCRIPT_FILE_PATH, "r", "utf8")   
        for line in F:
            line = line.rstrip('\r\n')
            data = line.rsplit('\t')
            transcript=data[0]
            ### add transcript to dictionary 
            if(transcript not in self.mapping_info["transcript"]):
                self.mapping_info["transcript"][transcript] = {}
                self.mapping_info["transcript"][transcript]["name"] = data[0]
                self.mapping_info["transcript"][transcript]["chr"] = data[1]
                self.mapping_info["transcript"][transcript]["start_position"] =  data[2]
                self.mapping_info["transcript"][transcript]["cigar"] = data[3]
                self._check_transcript_data(transcript)
                
                ### get more info from cigar
                if(transcript in self.mapping_info["transcript"] ):
                    self._get_transcript_info_from_cigar(transcript)
            ### ignore duplicates
            else:
                self._utils.warn_me("WARNING", " I already found a transcript names: '" + transcript + "' with information:\n" + json.dumps(self.mapping_info["transcript"][data[0]], sort_keys=True, indent=4) + "I will ignore the current information: '" + line)
        F.close()
    #---# END _open_transcript_file 
    
    ### Sanity check of a transcript data
    def _check_transcript_data(self, read_name):
        """ Sanity check of the transcript data 
            TODO add test
        """
        self._utils.warn_me("DEBUG", "Data sanity check for transcript '" + read_name + "'")
        ### check chr
        if(not (re.match(r'^chr(X|Y|[1-9][0-9]?)$', self.mapping_info["transcript"][read_name]["chr"], flags=re.IGNORECASE)) or  
           (re.match(r'^chr([1-9][0-9]?)$', self.mapping_info["transcript"][read_name]["chr"], flags=re.IGNORECASE) 
            and int(re.split(r'^chr', self.mapping_info["transcript"][read_name]["chr"], flags=re.IGNORECASE)[1]) > 22)):  # check chr
            self._utils.warn_me("WARNING", " I will ignore '" + read_name + "' because the chromosome information: '"+self.mapping_info["transcript"][read_name]["chr"]+"' was not of the form 'chr<num>' with 0 < num <=22:\n" + json.dumps(self.mapping_info["transcript"][read_name], sort_keys=True, indent=4))
            self.mapping_info["transcript"].pop(read_name)
        ### check start position
        elif(not re.match(r'^\d+$', str(self.mapping_info["transcript"][read_name]["start_position"]))):  # check position
            self._utils.warn_me("WARNING", " I will ignore '" + read_name + "' because the start position: '"+self.mapping_info["transcript"][read_name]["start_position"]+"' was not an integer >= 0:\n" + json.dumps(self.mapping_info["transcript"][read_name], sort_keys=True, indent=4))
            self.mapping_info["transcript"].pop(read_name)
        ### check cigar format
        elif(re.search(self._cigar_check, self.mapping_info["transcript"][read_name]["cigar"]) 
             or not re.match(self._cigar_pat, self.mapping_info["transcript"][read_name]["cigar"]) 
             or len(re.findall(r"\d+" , self.mapping_info["transcript"][read_name]["cigar"])) != len(re.findall(r"[MIDNSHP=X]" , self.mapping_info["transcript"][read_name]["cigar"]))
             or len(re.findall(self._cigar_pat , self.mapping_info["transcript"][read_name]["cigar"])) != len(re.findall(r"[MIDNSHP=X]" , self.mapping_info["transcript"][read_name]["cigar"]))
             or len(re.findall(self._cigar_pat , self.mapping_info["transcript"][read_name]["cigar"])) != len(re.findall(r"\d+" , self.mapping_info["transcript"][read_name]["cigar"]))          
             or(
                 len(re.findall(r"\d+" , self.mapping_info["transcript"][read_name]["cigar"])) == len(re.findall(r"[MIDNSHP=X]" , self.mapping_info["transcript"][read_name]["cigar"]))
                 and len(re.findall(self._cigar_pat , self.mapping_info["transcript"][read_name]["cigar"])) == len(re.findall(r"[MIDNSHP=X]" , self.mapping_info["transcript"][read_name]["cigar"]))
                 and len(re.findall(self._cigar_pat , self.mapping_info["transcript"][read_name]["cigar"])) == len(re.findall(r"\d+" , self.mapping_info["transcript"][read_name]["cigar"]))          
                 and re.search("\d+[H|S]", self.mapping_info["transcript"][read_name]["cigar"]) 
                 and not (re.search("^\d+[H|S]", self.mapping_info["transcript"][read_name]["cigar"]) or re.search("\d+[H|S]$", self.mapping_info["transcript"][read_name]["cigar"]))
                )
             ):  # check cigar
            self._utils.warn_me("WARNING", " I will ignore '" + read_name + "' because cigar: '" + self.mapping_info["transcript"][read_name]["cigar"] + "' doesn't follow cigar format\n" + json.dumps(self.mapping_info["transcript"][read_name], sort_keys=True, indent=4))
            self.mapping_info["transcript"].pop(read_name)
        else:
            self.mapping_info["transcript"][read_name]["start_position"]=int(self.mapping_info["transcript"][read_name]["start_position"])
        #---# END check chr, position, cigar        
    #---# END _check_transcript_data 
        
    ### Extracts info from cigar
    def _get_transcript_info_from_cigar(self,  transcript):
        """Extracts info from cigar"""
        self._utils.warn_me("DEBUG", "Extracts info for transcript '" + transcript + "' from cigar: "+ self.mapping_info["transcript"][transcript]["cigar"]+"\n"+json.dumps(self.mapping_info["transcript"][transcript], sort_keys=True, indent=4) )
        self.mapping_info["transcript"][transcript]["tpos"]={}
        chr=self.mapping_info["transcript"][transcript]["chr"]
        if(chr not in self.mapping_info["genome"]):
            self.mapping_info["genome"][chr]={}
       
        total=0
        tpos=0
        gpos=self.mapping_info["transcript"][transcript]["start_position"]
        nsymbol=0
        ### split cigar entries
        for centry in self._cigar_pat.findall(self.mapping_info["transcript"][transcript]["cigar"]):
            ccount = int(centry[:-1])
            csymbol = centry[-1]         
            tend=tpos+ccount-1
            gend=gpos+ccount-1
            nexttpos=tend+1
            nextgpos=gend+1

            # M 0 alignment match (can be a sequence match or mismatch). Adds positions to transcript and genomic coordinates
            # = 7 sequence match. Adds positions to transcript and genomic coordinates
            # X 8 sequence mismatch. Adds positions to transcript and genomic coordinates
            if( re.findall("^[M|=|X]$", csymbol) ):
                """ Adds positions to transcript and genomic coordinates """
                total += ccount

            # I 1 insertion to the reference. Adds positions to transcript coordinates
            elif( re.findall("^[I]$", csymbol) ):
                """ Adds positions to transcript coordinates """
                nextgpos=gpos
                gpos-=1
                gend=gpos   
                total += ccount
            # S 4 soft clipping (clipped sequences present in SEQ). Adds positions to transcript coordinates
            elif( re.findall("^S$", csymbol) ):
                """ Adds positions to transcript coordinates """
                nextgpos=gpos
                gpos=-999
                gend=-999   
                total += ccount
            # D 2 deletion from the reference. Adds positions to genomic coordinates
            # N 3 skipped region from the reference. Adds positions to genomic coordinates
            elif( re.findall("^[D|N]$", csymbol) ):
                """ Adds positions to genomic coordinates """
                nexttpos=tpos
                tpos-=1
                tend=tpos      
            # H 5 hard clipping (clipped sequences NOT present in SEQ). no position added
            elif( re.findall("^[H]$", csymbol) and nsymbol==0):
                """ no position added """
                nexttpos=tpos
                tpos-=ccount
                tend=nexttpos-1            
                nextgpos=gpos
                gpos-=ccount
                gend=nextgpos-1  
            elif( re.findall("^[H]$", csymbol) and nsymbol>0):
                """ no position added """
                nexttpos=tpos
#                 tpos=tpos
#                 tend=tpos+count-1            
                nextgpos=gpos
#                 gpos-=ccount
#                 gend=nextgpos-1 
            # P 6 padding (silent deletion from padded reference). no position added - ignore
            elif( re.findall("^[P]$", csymbol) ):
                """ no position added """
                nexttpos=tpos
                tpos-=1
                tend=tpos            
                nextgpos=gpos
                gpos-=1
                gend=gpos              

            ### initiate transcript info
            if(tpos not in self.mapping_info["transcript"][transcript]["tpos"]):
                self.mapping_info["transcript"][transcript]["tpos"][tpos]={}
            if(tend not in self.mapping_info["transcript"][transcript]["tpos"][tpos]):
                self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]={}
            self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]["type"]=csymbol
            self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]["count"]=ccount
            self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]["start_tpos"]=tpos 
            self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]["end_tpos"]=tend
            self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]["start_gpos"]=gpos 
            self.mapping_info["transcript"][transcript]["tpos"][tpos][tend]["end_gpos"]=gend

            ### initiate genomic info
            if(gpos not in self.mapping_info["genome"][chr]):
                self.mapping_info["genome"][chr][gpos]={}
            if(str(gend) not in self.mapping_info["genome"][chr][gpos]):
                self.mapping_info["genome"][chr][gpos][gend]={transcript:{}}
            self.mapping_info["genome"][chr][gpos][gend][transcript]["type"]=csymbol
            self.mapping_info["genome"][chr][gpos][gend][transcript]["count"]=ccount
            self.mapping_info["genome"][chr][gpos][gend][transcript]["start_gpos"]=gpos
            self.mapping_info["genome"][chr][gpos][gend][transcript]["end_gpos"]=gend
            self.mapping_info["genome"][chr][gpos][gend][transcript]["start_tpos"]=tpos
            self.mapping_info["genome"][chr][gpos][gend][transcript]["end_tpos"]=tend
            
            ### change current genomic/transcript positions
            tpos=nexttpos
            gpos=nextgpos
            nsymbol+=1
        #---# END split cigar loop

        self.mapping_info["transcript"][transcript]["length"]=total

    #---# END _get_transcript_info_from_cigar 

    ### Reads queries and outputs associated genomic coordinates
    def _answer_queries(self):
        """ Reads the file with queries
            outputs genomic coordinate corresponding at the transcript coordinate
        """
        self._utils.warn_me("INFO", " opening query file " + config.QUERY_FILE_PATH)
 
        OUT = self._utils.open_file(config.TRANSLATION_PATH, "w", "utf8")   
        F = self._utils.open_file(config.QUERY_FILE_PATH, "r", "utf8")   
        nl=0
        for line in F:
            nl+=1
            line = line.rstrip('\r\n')
            data = line.rsplit('\t')
            print(data)
            if(len(data)==2):
                transcript = data[0] 
                position = int(data[1])
                if(transcript not in self.mapping_info["transcript"]):
                    self._utils.warn_me("WARNING", " I could not find a transcript names: '" + transcript + "' in the file: " + config.TRANSCRIPT_FILE_PATH+". no data wil be given for " + line)
                    OUT.write(line+"\tNA\tNO DATA FOUND\n")           
                else:
#                     transcript_lenght=self.mapping_info["transcript"][transcript]["length"]
#                     if(position < 0 or position > transcript_lenght):
#                         self._utils.warn_me("WARNING","QUERY ISSUE: POSITION - At line# " + str(nl)+", query for transcript names: '" + transcript + ": QUERIED POSITION: "+str(position)+" OUTSIDE THE TRANSCRIPT OF LENGHT "+str(transcript_lenght))
#                         OUT.write(line+"\tNA\tout of bound\n")
#                     else:
                    (gpos,comment) = self._get_genomic_coordinate(transcript, position)
                    OUT.write(line+"\t"+str(gpos)+ "\t"+comment+"\n")
                    print(line+"\t"+str(gpos)+ "\t"+comment+"\n")
            else:
                self._utils.warn_me("WARNING","QUERY ISSUE: FORMAT - At line# " + str(nl)+", query : '" + line + " len(data)="+str(len(data))+". Make sure the query file is tab delimited")
                OUT.write(line+"\tQUERY ISSUE: FORMAT - len(data)= "+str(len(data))+". Make sure the query file is tab delimited\n")

        F.close()
        OUT.close()
    #---# END _open_transcript_file 
    
    
        
    ### Extracts genomic coordinate
    def _get_genomic_coordinate(self,  transcript,position):
        """ Returns the genomic coordinate for position 'position' in transcript 'transcript' 
            as well as the information associated with this coordinate
            TODO: add test
        """
        self._utils.warn_me("info", "Extracts genomic coordinate for for transcript '" + transcript + "' at position: "+ str(position))
        gpos=None
        comment=""
        add=""
        for tstart in sorted(self.mapping_info["transcript"][transcript]["tpos"]):
            print("start"+str(tstart))
            if(tstart <= position):
                
                for tend in sorted(self.mapping_info["transcript"][transcript]["tpos"][tstart]):
                    if( position <=tend ):
                        print(str(tstart)+" "+str(tend)+" "+json.dumps(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend], sort_keys=True, indent=4))
                        # M 0 alignment match (can be a sequence match or mismatch). Adds positions to transcript and genomic coordinates
                        if(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"] == "M"):
                            comment+=add+"aligned"
                            
                        # = 7 sequence match. Adds positions to transcript and genomic coordinates
                        if(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"] == "="):
                            comment+=add+"match"
                            
                        # X 8 sequence mismatch. Adds positions to transcript and genomic coordinates
                        if(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"] == "X"):
                            comment+=add+"mismatch"

                        # I 1 insertion to the reference. Adds positions to transcript coordinates
                        # S 4 soft clipping (clipped sequences present in SEQ). Adds positions to transcript coordinates
                        elif( re.findall("^[I|S]$", self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"]) ):
                            comment+=add 
                            comment+= str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["count"])
                            comment+=self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"]
                            comment+=": "
                            comment+=str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["start_tpos"])
                            comment+="-"
                            comment+=str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["end_tpos"])
                        # D 2 deletion from the reference. Adds positions to genomic coordinates
                        # N 3 skipped region from the reference. Adds positions to genomic coordinates
                        elif( re.findall("^[D|N]$", self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"])):
                            comment+=add
                            comment+=str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["count"])
                            comment+=self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"]+": "
                            comment+=str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["start_gpos"])
                            comment+="-"
                            comment+=str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["end_gpos"])   
                        # H 5 hard clipping (clipped sequences NOT present in SEQ). no position added
                        # P 6 padding (silent deletion from padded reference). no position added - ignore
                        elif( re.findall("^[H|P]$", self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"])):
                            comment+=add
                            comment+=str(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["count"])
                            comment+=self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"]
                            if(re.findall("^[H]$", self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"])):
                                 comment += ", out of bound"
                        if(add == ""):
                            add=", "
                            if(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"] != "H" and self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["start_gpos"] == self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["end_gpos"]):
                                gpos= self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["start_gpos"] 
                            elif(self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["type"] != "H" ):
                                gpos= self.mapping_info["transcript"][transcript]["tpos"][tstart][tend]["start_gpos"] + ( position-tstart)
                            
            #---# END loop tend
        #---# END loop tstart
        if(comment == ""):
            comment = "out of bound"

        if(gpos == -999 or gpos== None):
            gpos="NA"
        
        return (gpos, comment)
    #---# END _get_genomic_coordinate 
#---#END CoordinateTranslator
    
### main function
def main():
    """ Runs the transcript translator given the options """
    
    CT = CoordinateTranslator()
    CT._open_transcript_file()
#     print(CT.mapping_info)
    print(json.dumps(CT.mapping_info, sort_keys=True, indent=4))
    CT._answer_queries()
#---# END main
 
if __name__ == "__main__":
    main()  
    import doctest
    doctest.testmod()
