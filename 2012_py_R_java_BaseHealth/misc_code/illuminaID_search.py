#!/usr/bin/python
# Lists the Illumina ID folder available in a folder abd run L1L2prime script for those
#'''
#Created on Jan 26, 2012
#@author: celine
#'''

import os
import sys
import time
IND=''
IND='LP000-DNA_'
dirHD =os.getcwd( )
end='A01'
if len(sys.argv) >1:
    i=0
    for arg in sys.argv:
        if sys.argv[i] =='-d':
            dirHD = sys.argv[i+1]
        elif sys.argv[i] =='-s':
            IND = sys.argv[i+1]
        elif sys.argv[i] =='-e':
            end = sys.argv[i+1]
        i+=1

ilu = [i for i in os.listdir(dirHD) if i.startswith(IND) and ( i.endswith(end) ) ]
print "I will analyze MEMBERS ID:" ,IND, end,ilu
sys.stdout.flush()


for iluid in ilu:
    t0=time.clock()
    print "Analyzing MEMBER ID:" ,iluid
    sys.stdout.flush()
    from l1_l2_PRIME import main
    main(dirHD,iluid)
    t1=time.clock()
    print "DONE Analyzing MEMBERS ID: %s\t%d" ,(iluid,t1-t0)
    sys.stdout.flush()



