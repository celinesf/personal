#! /usr/bin/env python
# recover all info from csv (after clean up spaces issues) and update the Authors, Year, Jorunal, Pages, Title, Abstract, Mesh, URL from .bib reference table
# Make things consistant: 
# gender: mixed-Mixed, 
# Broad Ethnicty Asian
import os,sys, re, gzip,string
from numpy import *

DISEASE="HeartDiseaseArticle.bib"
DB="Publications_1.csv"
out='Heart_Curated_1.csv'
IN=open(DB)
sh=IN.readline()
IN.close
sh=sh.strip().split(",")
header={}
for s in sh:
        header[s] = s

curdir = os.getcwd() + '/' 

new=0;
end=0;
paper={} #is this a new paper?
BIB=open(DISEASE)


OUT=open(out,'w')
OUT.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,\"\"\n" % (header['\"PMID\"'],header['\"Authors\"'],header['\"Year\"'],header['\"Journal\"'],header['\"Volume\"'],header['\"Num\"'],header['\"Pages\"'],header['\"Title\"'],header['"Abstract"'],header['\"Mesh\"'],header['\"URL\"']))
for line in BIB:
       
        article = line.strip().split('@article{')
        art=article[0].replace('\"','\'') # remove double quotes
        tp = art.strip().split('},')
        end= len(tp)
        if new==1  :
                
                word =art.strip().split(' = ')
                value=word[1].strip().split('}')    
                tp= len(value) 
                value=value[0].strip().split('{')
                tp1=len(value) 
                
                if tp1>2 and tp>2:
                        value[1]=''
                        v=word[1].strip().split('}')
                        for s in v:
                                p=s.strip().split('{')
                                for t in p:
                                        if t!='':
                                                value[1]="%s%s" % (value[1],t)
                                                
                paper[word[0]]=value[1]
       
        if end==1 and len(article)==1 and new==1:                
                st= [s for s in file(DB) if paper['Pmid'] in s]
                if len(st)==1:
                        s=st[0].strip().split(',')
                        if 'Au' in paper:
                                s[1]='\"%s\"' % paper['Au']
                        elif 'Author' in paper:
                                s[1]='\"%s\"' % paper['Author']
                        
                        if 'Year' in paper:
                                s[2]= paper['Year']
                 
                        if 'Journal' in paper:
                                s[3]= '\"%s\"' % paper['Journal']
        
                        if 'Volume' in paper:
                                s[4]= '\"%s\"' % paper['Volume']
                        else:
                                s[4]= '\"\"'   
                        if 'Number' in paper:
                                s[5]= '\"%s\"' % paper['Number']
                        else:
                                s[5]= '\"\"' 
                        if 'Pages' in paper:
                                s[6]= '\"%s\"' % paper['Pages']
 
                        if 'Title' in paper:
                                s[7]= '\"%s\"' % paper['Title']

                        if 'Abstract' in paper:
                                s[8]= '\"%s\"' % paper['Abstract']
                        
                        if 'Keywords' in paper:
                                s[9]= '\"%s\"' % paper['Keywords']
                        elif 'Mesh' in paper:
                                s[9]= '\"%s\"' % paper['Mesh']

                        if 'Url' in paper:
                                s[10]= '\"%s\"' % paper['Url']
                        elif 'Doi' in paper:
                                s[10]= '\"%s\"' % paper['Doi']

                        j = 0
                        for w in s:
                                j += 1
                                OUT.write("%s," % w)
                                if j > 10: 
                                        break
                        OUT.write("\"\"\n")
                else:
                        print "PROBLEM, The paper pmid %s was not in %s" % (paper['Pmid'],DB)
                new=0
                paper={}

  
             
                  
        if art=='' and len(article)>1 and new==0:
                new=1
               
       

     
OUT.close              
BIB.close
