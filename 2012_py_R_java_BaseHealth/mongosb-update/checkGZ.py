#!/usr/bin/python
'''
Created on October 9, 2012
@author: celine

update 2.1 genomics documents (in json format)
- add real and relative position
- add gene name
- update/clean effect
'''
import gzip

path ="/Users/celine/Box Documents/Celine/ScriptInOut/mongodb-update/"
filename = "PG0000802-BLD.snps.vcf.gz"

f = gzip.open('%s%s' % (path,filename),'r')

print f.readline()