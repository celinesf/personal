#!/usr/bin/python
#'''
#Created on Jan 27, 2012
#
#@author: celine
#'''
# command line Can have the following arguments:
# -i path_to_Illumina_HD (current directory by default)
# -b path_to_Amazon_HD (current directory by default)
# -g path_to_file_GID/IlluminaID (file named 'GID_Illumina.txt' in current directory by default)
# -o path_to_output_file  (file named 'GID_ILU_Gender.txt' in current directory by default)
#-e end_of_filenames_to_copy ('A01' by default)
# 
# The script reads the tab delimited text file with the corresponding GID and illumina ID 
# It then assesses whether this Illumina ID data corresponds to a male or a female and output that information in the file 'path_to_output_file'
# It creates a folder names GID in 'path_to_Amazon_HD' for each Illumina ID  ending with 'end_of_filenames_to_copy' in the HD 'path_to_Illumina_HD'
# It will copy the Illumina ID folders ending with 'end_of_filenames_to_copy' from the HD 'path_to_Illumina_HD'
# into the folders of corresponding GID in 'path_to_Amazon_HD'


import os
import sys
import re
import linecache
import time

ID=''
dirHD = os.getcwd( )
DIRbckp =  os.getcwd( ) 
fout = 'GID_ILU_Gender.txt' #GID_ILU_Gender1test'   # list with gender after copying         
fin= 'GID Illumina link.txt' #GID_Illumina.txt' # table of GID Illu id from malekeh    

# get arguments from command line
if len(sys.argv) >1:
    i = 0
    for arg in sys.argv:
        if arg == '-i':  dirHD = sys.argv[i+1]
        elif arg == '-b':  DIRbckp = sys.argv[i+1]  
        elif arg == '-g': fin = sys.argv[i+1]  
        elif arg == '-o':  fout = sys.argv[i+1]
        elif arg == '-e': ID = sys.argv[i+1]
        i += 1
        
print "File with GID/illu ID:", fin
print "Path to Illumina HD:",dirHD
print "Path to Amazon HD:", DIRbckp
print "Output file with genders:", fout
print "End of file concidered:", ID 
fileIN =    open(fin)
data = fileIN.read()
fileIN.close() 

ilu = [i for i in os.listdir(dirHD) if i.endswith(ID) ]
print "I will analyze MEMBERS ID:" ,ilu
fileOUT = open(fout,'w')
for iluid in ilu:
    if re.search('%s' % iluid.lower(), data.lower()) :
        nf = 0 # find the iluid once
        for m in re.finditer('%s' % iluid.lower(), data.lower()):
            t0=time.clock()
            nf += 1
            if nf > 1: 
                print "DUPLICATE iluID %s" % iluid
            start = m.start()
            lineno = data.count('\r', 0, start) + 1
            word = m.group(0)
            line = linecache.getline(fin, lineno)
            word =(line.strip().split()) 
            GID = word[3]
            dirdbsnp= '%s/%s/Variations/dbSNP/' % (dirHD,iluid)
            startf = 'chr'
            endf = '_custom_dbSNP.txt'
            sex = 0
            chrfile = '%s%s%s' % (startf, 'Y', endf)
            fnin = '%s%s' % (dirdbsnp,chrfile)
            try:
                fileY =    open(fnin)
                fileY.close()
                sex = 1
                fileOUT.write("%s\t%s\tMALE\n" % (GID,iluid))
            except :
                fileOUT.write("%s\t%s\tFEMALE\n" % (GID,iluid))
                sex = 0
            fileOUT.flush()
            
            sys.stdout.flush() 
            #if not os.path.exists("%s/%s" % (DIRbckp,GID)):
            #    print "I am creating the folder  %s into %s" % (GID,DIRbckp)
            #    os.system("mkdir %s/%s/" % (DIRbckp,GID))
            #else:
            #    print "I already have %s/%s" % (DIRbckp,GID)
            #print "I am copying %s into %s" % (iluid, GID)
            #os.system("cp -r %s/%s %s/%s/" % (dirHD,iluid,DIRbckp,GID))  
            t1=time.clock()
            print "Done copying %s into %s. Time %d" % (iluid, GID, (t1-t0))   
            sys.stdout.flush()           
fileOUT.close()
         
