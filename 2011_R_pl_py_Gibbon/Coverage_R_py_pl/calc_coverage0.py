#! /usr/bin/env python

# calculate length of BACs, then estimate coverage for each bac based on the reads.

import os,sys
from numpy import *

name='Jo'
f=[1,2]

print '%s %s\n' % (f[0], f[1])

gibbondir = '/Users/celine/0_data/0Work_2010/Dropbox/Gibbon2/ChrisCode/'

# start by calculating the length of each bac, then initialize each into a vector of 0's

# parse fasta file
fastafile = gibbondir + 'nomascus2.fasta'
bacseqdict = {}
for line in file(fastafile):
	# skip missing lines
	if len(line.strip()) > 0:		
		if line.startswith('>'):
			# then BAC name
			bacname = line.strip().replace('>','')
			bacseqdict[bacname] = []
		else:
			bacseqdict[bacname] += line.strip()
		
# now, initialize vector for each bac name
bacdict = {}
infobac = {}
#bacfile = file('BACsizes'  , 'w')
for bacname in bacseqdict:
	bacdict[bacname] = zeros((len(bacseqdict[bacname])),dtype=int)
	infobac[bacname] = zeros(2,dtype=int)
        infobac[bacname][0]=len(bacseqdict[bacname])
     #   bacfile.write('%s\t%d\n' % (bacname, len(bacdict[bacname])))
#bacfile.close()

filelist = [i for i in os.listdir(gibbondir) if i.endswith('-b') and (i.startswith('t_%s'% (f[0])) or i.startswith('s_%s'% (f[1])))]

print '%s files total...' % (len(filelist))

for i,filename in enumerate(filelist):
	print '%s %s\n' % (filename, len(filelist) - i),
	sys.stdout.flush()
	# loop through each line in the file here, find mapped bac and bpstart
	for line in file(gibbondir + filename):
		line = line.strip().split()
		bacname = line[2]
		bpstart = int(line[3])
		bacdict[bacname][bpstart:bpstart + 100] += 1
                if infobac[bacname][0]>bpstart:  infobac[bacname][0]=bpstart
                if infobac[bacname][1]<bpstart+100:  infobac[bacname][1]=bpstart+100

bacfile = file('BACsizes_'+name  , 'w')
for bacname in bacseqdict:
        bacfile.write('%s\t%d\t%d %d\n' % (bacname, len(bacdict[bacname]),infobac[bacname][0],infobac[bacname][0] ))
bacfile.close()                

toto='gibbon_coverage_'+name 
print '%s\n' % (toto)
outfile = file('gibbon_coverage_'+name  , 'w')
for bacname in bacdict:
	print '%s\n' % (bacname),
	sys.stdout.flush()
	outfile.write( '%s\t%s\n' % (bacname,'\t'.join([str(val) for val in bacdict[bacname]])))
outfile.close()
print 'done!'
