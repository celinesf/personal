#! /usr/bin/env python
import logging,   io, operator, re

OUTPUTPATH = "/Users/celine/Box Documents/Celine/ScriptInOut/NextBio"
DATAPATH = "%s/08-14-13_cleanAt" % OUTPUTPATH
batchname = "gwasbiosets_aug2013_CB4.txt"

input_data = io.open("%s/%s" % (DATAPATH,batchname),'r',encoding='utf-8-sig')
output = open('%s/bioset_4270.txt' % OUTPUTPATH,'w')
print '%s/bioset_4270.txt' % OUTPUTPATH
line = input_data.readline() ### first line ========
output.write(line)
while len(line)>0:
    line = input_data.readline() 
    if '4270\t22064162' in line:
        output.write(line)
        
input_data.close()
output.close()
        