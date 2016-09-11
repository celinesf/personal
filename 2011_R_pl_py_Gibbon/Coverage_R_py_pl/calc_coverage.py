# calculate length of BACs, then estimate coverage for each bac based on the reads.

import os,sys
from numpy import *

gibbondir = '/home/celine/Gibbon2/'

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
for bacname in bacseqdict:
	bacdict[bacname] = zeros((len(bacseqdict[bacname])),dtype=int)

# now, parse files
filelist = [i for i in os.listdir(gibbondir) if i.endswith('-b')]
#filelist = ['s_8_008_qseq.txt-b'] # small test list

print '%s files total...' % (len(filelist))
for i,filename in enumerate(filelist):
	print '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s %s' % (filename, len(filelist) - i),
	sys.stdout.flush()
	# loop through each line in the file here, find mapped bac and bpstart
	for line in file(gibbondir + filename):
		line = line.strip().split()
		bacname = line[2]
		bpstart = int(line[3])
		bacdict[bacname][bpstart:bpstart + 100] += 1

print '\n###\nwriting output...'
outfile = file('gibbon_coverage_by_bac.txt'  , 'w')
for bacname in bacdict:
	print '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s' % (bacname),
	sys.stdout.flush()
	outfile.write( '%s\t%s\n' % (bacname,'\t'.join([str(val) for val in bacdict[bacname]])))

outfile.close()
print 'done!'