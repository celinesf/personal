#!/usr/bin/env python

"""
    split GWAS_standard9.txt
   05/16/13 - 1.0
   06/06/13: rerun all for CB2
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 

"""
    1) clean NextBio data from non unicode characters       
"""

import  os, json, re, copy
import NextBioConfig as config

''' get map for typo corrections '''
def getTermMap(fname):
    f = open(fname,'r')
    typos = []
    nl = 1
    for line in f:
        s = line.strip()
        words = s.split("*\t*")
        if len(words)>1:
            #typo is the encoded character.
            #rep is the character that we replace the type character
            #num is the number of first line it was found.
            #count is the count of times the typo character was found.
            #values is the text where the character was found.
            toadd =  {"typo":words[0].split("*")[1], "rep":words[1].split("*")[0], 'num':nl,"count":0, "values":[]}
            typos.append(toadd)
        nl +=1
    f.close()
    return typos #

def identify_errors(batchname, code_map):
    input_data = open("%s/%s.txt" % (config.DATAPATH_ORI,batchname),'r')

    line = input_data.readline() ### first line ========
    nl=1
    while len(line)>0:
        line = input_data.readline() 
        new_line = line.split("\n")[0]
        nc = 1
        for typos in code_map:
            if typos["typo"]  in new_line:
                typos["count"] += new_line.count(typos["typo"])
                if new_line not in typos["values"]:
                    typos["values"].append(new_line)
            nc+=1
        nl+=1
        if (nl % 50000) == 0 :
            print nl
    input_data.close()  

    code_map = find_new_errors(code_map)
    
    output_file = open("%s/typos2_%s.json" % (config.OUTPUTPATH,batchname) ,'w') 
    output_file.write(json.dumps(code_map, indent=4)) 
    output_file.close()
    
def find_new_errors(code_map):
    new_data = copy.deepcopy(code_map)
    for data in new_data:
        new_values = []
        for line in data['values']:
            new_line = line
            for typo in new_data:
                if typo['typo'] in line:
                    new_line = re.sub("\%s" % typo["typo"],typo["rep"], new_line)
            new_values.append(new_line)     
        data['values'] = copy.deepcopy(new_values) 
    return new_data

    
def correct_errors(batchname, code_map):
    input_data = open("%s/%s.txt" % (config.DATAPATH_ORI,batchname),'r')
    output_file = open("%s/%s_v8.txt" % (config.OUTPUTPATH,batchname),'w')

    line = input_data.readline() ### first line ========
    nl=1
    nt =0 
    output_file.write(line)
    while len(line)>0:
        line = input_data.readline() 
        nc = 1
        for typos in code_map:
            if typos["typo"]  in line:
                typos["count"] += line.count(typos["typo"])
                nt += line.count(typos["typo"])
                if line not in typos["values"]:
                    typos["values"].append(line)
                line = line.replace(typos["typo"], typos["rep"])
            nc+=1
        nl+=1
        if (nl % 50000) == 0 :
            print nl
        output_file.write(line)
    print nt
    input_data.close()  
    output_file.close()

""" MAIN function """    
if __name__ == "__main__":

    batchname ="gwasbiosets_aug2013_3items_update"
    code_map = getTermMap("%s" % (config.CODEMAP))
    print batchname
    identify_errors(batchname, code_map)         
    correct_errors(batchname, code_map)
    print'DONE with DMmain'
    
""" END OF main """ 