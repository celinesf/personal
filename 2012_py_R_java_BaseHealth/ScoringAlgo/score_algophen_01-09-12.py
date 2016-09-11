#!/usr/bin/python
############ ALgo
# to is the total
# pm number pmid /where CI doesnt include 1 and or PValue<0.05 /total pmid
# et number BroadEthnicity rep, es ratio sign/ns
# st Is gwas from SNPInfo 'yes' ->1, if Comments contain 'meta' (SNPInfo), or if SNPOR Comemnts contain '/' -> 2. otherwise 0 
# na Caseno/500, /1000, /200
# sa SampleSize/1000, /2000, /400
# di NextBioDisease countain disease name (make list of possible match)


# re Comments contain ns and '[' count ',' -> -n-1 can '[' count , only ->n+1 + 
# mt ratio from metanalyses else 0
###### later l Comments contain 'LD' and not 'NOLD' [check in SNPLD how many snp are replicated?]
# gp number SNP for GenePaper if not NULL or NA
# gn number of SNP for GeneNextBio is not NBNA or NBNA

import MySQLdb
import string
import re
import sys

disease=1
if len(sys.argv) >1:
        disease = string.atoi(sys.argv[1])

NextBio ={}
if disease == 1:
        NextBio = ['Diabetes mellitus', 'Insulin','endocrine pancreas','Fasting plasma glucose level']
elif disease == 4:
        NextBio = ['intestine']
elif disease == 5:
        NextBio = ['Dementia','Alzheimer']
elif disease == 6  :
        NextBio = ['cholesterol','Triglyceride','Body mass index','Cardiomyopathy','cardiac function','Heart diseases','Cardiovascular','Coronary artery']
elif disease == 7  :
        NextBio = ['Cerebrovascular','Cardiovascular','Vascular']
elif disease == 8:
        NextBio = ['Mood disorder']
elif disease == 9:
        NextBio = ['blood pressure','Hypertensive']
elif disease == 103:
        NextBio = ['Breast']
elif disease == 105:
        NextBio = ['Alopecia']
elif disease == 106:
        NextBio = ['Lung cancer','airway']
elif disease == 107:
        NextBio = ['Migraine']
elif disease == 108:
        NextBio = ['Body mass index', 'Obesity', 'Weight','Waist']
elif disease == 109:
        NextBio = ['mineral density','Arthropathies','musculoskeletal','Mineral Density']
elif disease == 110:
        NextBio = ['Prostate']
elif disease == 112 or disease == 104:
        NextBio = ['Colitis','Inflammatory bowel','Crohn']


# Open database connection
db = MySQLdb.connect(host="192.168.2.104",user="celin",passwd="celin",db="AlgoPhen")
# prepare a cursor object using cursor() method
cursor = db.cursor()

template_score = "to= pm= et= st= ca= sa= mt= ld= gp= gn= di= re= fr= "
def signOR(ntot,ns,curs):
        ORv={}
        snpor = curs.fetchall()
        desc = curs.description
        ntot += 1
        founds=0;
        for OR in snpor: ### check all CI
                for (name, value) in zip(desc, OR) :
                        ORv[name[0]] = value   
                pval = ORv['Pvalue'] 
                ci = ORv['CI']
               # print pval
                #print ci
                if founds == 0:
                        if  ci is not None:
                               # print 'y',ci
                                c0=string.atof(ci.split("-", 1)[0])
                                c1=string.atof(ci.split("-", 1)[1])
                                if(c1 >1 and  c0> 1)  or (c0 <1 and c1 < 1): 
                                        founds +=1
                        elif pval is not None and  pval <=0.05:
                               #print 'z', pval
                                founds +=1
                       # elif  pval is None and ci is not None:
                      
                else:
                        break  
        ns += founds
        return (ntot,ns)


############################## SNP specific ######################.
sqlu = "UPDATE SNPOR set Score = '%s' WHERE  DiseaseID=%d " % (template_score,disease)
try: cursor.execute(sqlu)
except: print "Error: unable to fecth data %s" % sqlu


sql = "SELECT DISTINCT dbsnp FROM SNPOR WHERE DiseaseID = %d " % (disease)
#print sql
try:
        # Execute the SQL command
        cursor.execute(sql)
        # Fetch all the rows in a list of lists.
        dbsnp = cursor.fetchall()
        print "Scoring disease %d, found %d SNPs " %(disease, len(dbsnp))
        for snp in dbsnp: ### loop on SNPs
                
############### Disease associated #########################
                sqld=  "SELECT *  FROM SNPOR WHERE DiseaseID = %d AND dbsnp='%s'" % (disease,snp[0])
                #print sqld
                try:
                        cursor.execute(sqld)
                        row= cursor.fetchall()
                        desc = cursor.description
                        ORv1 = {}
                        for (name, value) in zip(desc, row[0]) :
                                ORv1[name[0]] = value
                        nd=0;
                        if ORv1['NextBioDisease'] != 'NBNA' and ORv1['NextBioDisease'] != "NBNF"  and   ORv1['NextBioDisease'] is not None:
                                for d in NextBio:
                                   nd += len(re.split(d,ORv1['NextBioDisease'])) -1
                                if nd>0:
                                        sqlu = "UPDATE SNPOR set Score = replace(Score,'di= ','di=%d ')\
WHERE DiseaseID=%d  AND dbsnp='%s'" % (nd,disease,snp[0])
                                        try: cursor.execute(sqlu)
                                        except: print "Error: unable to fecth data %s" % sqlu
                except:  print "Error: unable to fecth data %s" % sqld

##################### LD replication ######################
                nld =  ""
                sqll= "SELECT *  FROM SNPLD WHERE (dbsnp1='%s' and dbsnp2 in (select dbsnp from SNPOR where diseaseid=%d)) or (dbsnp2='%s' and dbsnp1 in (select dbsnp from SNPOR where diseaseid=%d))" % (snp[0],disease,snp[0],disease)
                #print sqll
                try:
                        cursor.execute(sqll)
                        row= cursor.fetchall()
                        snpl={}
                        str=""
                        n=0
                        nl=nlsign=0
                        for ld in row:
                                if n==0: str = "%s'%s','%s'" % (str,ld[0],ld[1])
                                else:   str = "%s,'%s','%s'" % (str,ld[0],ld[1])
                                n+=1
                        if n>0:
                                sqls= "SELECT distinct dbsnp  FROM SNPOR WHERE dbsnp in (%s) AND diseaseid=%d" % (str,disease)
                                #print sqls
                                try:
                                        cursor.execute(sqls)
                                        sn= cursor.fetchall()
                                        for s in sn:
                                                sqlor= "SELECT *  FROM SNPOR WHERE dbsnp = '%s' AND diseaseid=%d" % (s[0],disease)
                                        #print sqlor
                                                try:
                                                        cursor.execute(sqlor)
                                                        nl,nlsign = signOR(nl,nlsign,cursor)
                                                        if nl>0:
                                                                nld = "%d/%d" % (nlsign,nl)
                                                except:  print "Error: unable to fecth data %s" % sqlor                                       
                                except:  print "Error: unable to fecth data %s" % sqls
                except:  print "Error: unable to fecth data %s" % sqll

########################## Count ethnicity significant/num ethnicity given & find 
                neth = ""
                sqle=  "SELECT DISTINCT  BroadEthnicity FROM SNPOR WHERE DiseaseID = %d AND dbsnp='%s'" % (disease,snp[0])
                try:
                        cursor.execute(sqle)
                        eth= cursor.fetchall()
                        nesign=0;
                        ne=0
                        for e in eth: ## Loop on ethnicities
                                sqlor= "SELECT * FROM SNPOR WHERE DiseaseID=%d AND BroadEthnicity='%s' AND dbsnp='%s'" % (disease,e[0],snp[0])
                                #print sqlor
                                try:
                                        cursor.execute(sqlor)
                                        ne,nesign = signOR(ne,nesign,cursor)
                                        if ne>0:
                                                neth = "%d/%d" % (nesign,ne)
                                except: print "Error: unable to fecth data %s" % sqlor    
                except:  print "Error: unable to fecth data %s" % sqle


########################## Count PMID significant/num
                npm=""               
                sqlp=  "SELECT DISTINCT  pmid FROM SNPOR WHERE DiseaseID = %d AND dbsnp='%s'" % (disease,snp[0])
               # print sqlp
                try:
                        cursor.execute(sqlp)
                        pmid= cursor.fetchall()
                        desc = cursor.description
                        npsign=0;
                        np=0
                        for p in pmid: ## Loop on pmid
                                sqlor= "SELECT * FROM SNPOR WHERE DiseaseID=%d AND PMID=%d AND dbsnp='%s'" % (disease,p[0],snp[0])
                               # print sqlor, p
                                try:
                                        cursor.execute(sqlor)
                                        np,npsign = signOR(np,npsign,cursor)
                                        if np>0:
                                                npm = "%d/%d" % (npsign,np)
                                except:   print "Error: unable to fecth data %s" % sqlor
                except:
                        print "Error: unable to fecth data %s" % sqlp
                  
########################## Count replications
                nrep = ""
                sqlc=  "SELECT * FROM SNPOR WHERE DiseaseID = %d AND dbsnp='%s' AND comments is not NULL" % (disease,snp[0])
                #print sqlc
                try:
                        cursor.execute(sqlc)
                        row = cursor.fetchall()
                        desc = cursor.description
                        NF=sig=nots=0
                        mt=''
                        mtden=mtnum=0
                        for co in row:
                                ORv2 = {}
                                for (name, value) in zip(desc, co) :
                                        ORv2[name[0]] = value
                                if ORv2['comments'] is not None:
                                        if len(re.split('NF',ORv2['comments']))>1: NF += 1
                                        if len(re.split('nr',ORv2['comments']))>1 :
                                                nots+= len(re.split('nr',ORv2['comments']))-1
                                        if len(re.split('ns',ORv2['comments']))>1 : 
                                                nots+=len(re.split('ns',ORv2['comments']))-1
                                        if len(re.split('sg',ORv2['comments']))>1 : 
                                                sig +=len(re.split('sg',ORv2['comments']))-1
                                        com=re.split(' |,|-|:|\.|\+|\[|\]',ORv2['comments'])
                                                        ### meta numbers in Comments
                                        for c in com:
                                                if len(re.split('/',c))==2 :
                                                        mtden +=  string.atoi(re.split('/',c)[1])
                                                        mtnum +=  string.atoi(re.split('/',c)[0])  
                        if mtden>0 :
                                mt = "%d/%d" % (mtnum,mtden)
                        if sig >0 and NF>0:
                                print "Error: I found NF and significant for %s\n\t%s" % (sqlc,ORv2['comments'])
                        if nots+sig>0: nrep = "%d-%d" % (sig,nots) 
                except:
                        print "Error: unable to fecth data %s" % sqlc
          
                sqlu = "UPDATE SNPOR set Score = replace(Score,'et= ','et=%d/%d '), \
 Score = replace(Score,'pm= ','pm=%d/%d '), \
 Score = replace(Score,'ld= ','ld=%s '), \
 Score = replace(Score,'re= ','re=%s ') \
WHERE DiseaseID=%d  AND dbsnp='%s'" % (nesign,ne,npsign,np,nld,nrep,disease,snp[0])
                try:
                        while cursor.execute(sqlu) >0 :
                                cursor.execute(sqlu)
                except: print "Error: unable to fecth data %s" % sqlu

                                       ########################### Information on PMID and SNP ###############
                sqlp=  "SELECT DISTINCT PMID,BroadEthnicity FROM SNPOR WHERE DiseaseID = %d AND dbsnp='%s'" % (disease,snp[0])
                try:
                        cursor.execute(sqlp)
                        pmideth =  cursor.fetchall()
                        npe = 0
                        for p in pmideth: ## Loop on PMID to recover informations                       
                                sqlor= "SELECT * FROM SNPOR WHERE DiseaseID = %d AND PMID = %d AND BroadEthnicity='%s' AND dbsnp='%s'" % (disease,p[0],p[1],snp[0])
                               # print sqlor
                                caseno=samplesize=0
                                try:
                                        hapmap = {}
                                        ok=0
                                        if p[1] == 'Caucasian': hapmap['CEU']=0
                                        elif p[1] == 'Asian':
                                                hapmap['CHB']=0
                                                hapmap['CHD']=0
                                                hapmap['JPT']=0
                                        elif p[1] == 'AfricanAmerican': hapmap['ASW']=0
                                        elif p[1] == 'Hispanic': hapmap['MEX']=0
                                        else:
                                                ok = 1
                                                hapmap[p[1]] = 0
                                        ############### do I have HAPMAP frequncies #########################
                                        if ok == 0 :
                                                if len(snp[0].split(',')) == 1:
                                                        info = {}
                                                        sqli=  "SELECT * FROM SNPInfo WHERE dbsnp ='%s'" % (snp[0])
                                                        try:
                                                                cursor.execute(sqli)
                                                                ri= cursor.fetchall()
                                                                desc = cursor.description
                                                                for (name, value) in zip(desc, ri[0]) :
                                                                        info[name[0]] = value
                                                                for h in hapmap:
                                                                        freq = 'Freq_%s' % h
                                                                        if info[freq] != None:
                                                                                hapmap[h] += 1
                                                        except:  print "Error: unable to fecth data %s" % sqli
                                        cursor.execute(sqlor)
                                        snpor = cursor.fetchall()
                                        desc = cursor.description
                                        no = 0
                                        freq =0
                                        for OR in snpor: ### check all CI and Pvalues for this pmid
                                                ORv3 = {}
                                                for (name, value) in zip(desc, OR) :
                                                        ORv3[name[0]] = value
                                                 ### record paper specific info caseno, sampleno, gene,gwas? 
                                                if no==0 and ORv3['CaseNo'] is not None and ORv3['SampleSize'] is not None: 
                                                        caseno = string.atof(ORv3['CaseNo'])
                                                        samplesize = string.atof(ORv3['SampleSize'])
                                                elif ORv3['CaseNo'] is None or ORv3['SampleSize'] is None:
                                                        print "Error: Null caseno and samplesize %s" % sqlor
                                                no += 1
                                                ### check if freq control here
                                                if freq == 0 and (ORv3['FreqControl'] is not None or ORv3['FreqPop'] is not None):
                                                        for h in hapmap:
                                                                hapmap[h] += 2
                                                                freq += 1
                                                                break
                                        fr = 0
                                        for h in hapmap:
                                                fr += hapmap[h]
                                except: print "Error: unable to fecth data %s" % sqlor
                                npe += 1
                                ####### Infor from Publications #######
                                sqlj= "SELECT * FROM Publications WHERE PMID = %d" % (p[0])
                                #print sqlj
                                try:
                                        cursor.execute(sqlj)
                                        study = cursor.fetchall()
                                        desc = cursor.description
                                        PMIDinfo={}
                                        for (name, value) in zip(desc, study[0]) :
                                                PMIDinfo[name[0]] = value
                                        st=0
                                        if PMIDinfo['GWAS'] == 'yes': st +=1
                                        if PMIDinfo['Experiment'] is not None and len(re.split('meta',PMIDinfo['Experiment']))>1: st +=2
                                except: print "Error: unable to fecth data %s" % sqlj
                                cano=sano=fq=''
                                if fr > 0:
                                        fq = "%d" % fr
                                if caseno >0  and samplesize>0:
                                        if st>=2: 
                                                caseno =caseno / 1000
                                                samplesize /= 2000
                                        elif st==1:
                                                caseno =caseno /  500
                                                samplesize /= 1000
                                        else: 
                                                caseno =caseno /  200
                                                samplesize /= 400
                                        cano= "%.2f" % caseno
                                        sano= "%.2f" % samplesize
                                sqlu = "UPDATE SNPOR set Score = replace(Score,'mt= ','mt=%s '), \
 Score = replace(Score,'st= ','st=%d '),\
 Score = replace(Score,'ca= ','ca=%s '),\
 Score = replace(Score,'sa= ','sa=%s '),\
 Score = replace(Score,'fr= ','fr=%s ') \
WHERE DiseaseID=%d AND PMID= %d AND BroadEthnicity='%s' AND dbsnp='%s'" % (mt,st,cano,sano,fq,disease,p[0],p[1],snp[0])
                                #print sqlu
                                try:
                                        while cursor.execute(sqlu) >0 :
                                                cursor.execute(sqlu)
                                except: print "Error: unable to fecth data %s" % sqlu
                except: print "Error: unable to fecth data %s" % sqlp
except: print "Error: unable to fecth data %s" % sql
######################### END SNP SPECIFIC ##########################


################ replication GenePaper ###################
sqlgp = "SELECT DISTINCT Genepaper FROM SNPOR WHERE DiseaseID = %d AND GenePaper is not NULL AND GenePaper != 'NA' AND GenePaper != 'intergenic' " % (disease)
try:
        cursor.execute(sqlgp)
        genep = cursor.fetchall()
        snpgp=''
        for g in genep: ## Loop on genepaper
                sqls= "SELECT distinct dbsnp FROM SNPOR WHERE DiseaseID = %d AND GenePaper like '%s%s%s'" % (disease,'%',g[0],'%')
                ngsign= ng=0
                try:
                        cursor.execute(sqls)
                        gsnp = cursor.fetchall()
                        for s in gsnp: ## Loop on ethnicities
                                sqlor= "SELECT * FROM SNPOR WHERE DiseaseID=%d AND GenePaper like '%s%s%s' AND dbsnp='%s'" % (disease,'%',g[0],'%',s[0])
                                try:
                                        cursor.execute(sqlor)
                                        ng,ngsign = signOR(ng,ngsign,cursor)
                                except: print "Error: unable to fecth data %s" % sqlor  
                except: print "Error: unable to fecth data %s" % sqls 
                if(ng>0):
                        snpgp='%d/%d' %(ngsign,ng)
                        sqlu = "UPDATE SNPOR set Score = replace(Score,'gp= ','gp=%s ') \
WHERE DiseaseID=%d  AND GenePaper='%s'" % (snpgp,disease,g[0])
                        try:
                                while cursor.execute(sqlu) >0 :
                                        cursor.execute(sqlu)
                        except: print "Error: unable to fecth data %s" % sqlu
except: print "Error: GenePaper  %s" % sqlgp
######################### END GenePaper ##########################

################ replication GeneNextBio #####################
sqlgn = "SELECT DISTINCT GeneNextBio FROM SNPOR WHERE DiseaseID = %d AND GeneNextBio is not NULL AND GeneNextBio != 'NBNA' AND GeneNextBio != 'NBNF'" % (disease)
try:
        cursor.execute(sqlgn)
        genen = cursor.fetchall()
        snpgn=''
        for g in genen: ## Loop on genepaper
                sqls= "SELECT distinct dbsnp FROM SNPOR WHERE DiseaseID = %d AND  GeneNextBio like '%s%s%s'" % (disease,'%',g[0],'%')
                ngsign= ng=0
                try:
                        cursor.execute(sqls)
                        gsnp = cursor.fetchall()
                        for s in gsnp: ## Loop on ethnicities
                                sqlor= "SELECT * FROM SNPOR WHERE DiseaseID=%d AND GeneNextBio like '%s%s%s' AND dbsnp='%s'" % (disease,'%',g[0],'%',s[0])
                                try:
                                        cursor.execute(sqlor)
                                        ng,ngsign = signOR(ng,ngsign,cursor)
                                except: print "Error: unable to fecth data %s" % sqlor  
                except: print "Error: unable to fecth data %s" % sqls 
                if(ng>0):

                        snpgn='%d/%d' %(ngsign,ng)
                        sqlu = "UPDATE SNPOR set Score = replace(Score,'gn= ','gn=%s ') \
WHERE DiseaseID=%d  AND GeneNextBio='%s'" % (snpgn,disease,g[0])
                        try:
                                while cursor.execute(sqlu) >0 :
                                        cursor.execute(sqlu)
                        except: print "Error: unable to fecth data %s" % sqlu
except: print "Error:GeneNextBio %s" % sqlgn
######################### END GeneNexBio ##########################


 ################ Total score #####################
sqlsc = "SELECT * FROM SNPOR WHERE DiseaseID = %d " % (disease)
#print sqlsc
try:  
        cursor.execute(sqlsc)
        ORrow = cursor.fetchall()
        desc = cursor.description
        for s in ORrow: ## Loop OR raws
                ORcol={}
                for (name, value) in zip(desc, s) :
                        ORcol[name[0]] = value
                if len(re.split("to= ",ORcol['Score']))>1:
                        list=re.split("to=| pm=| et=| st=| ca=| sa=| mt=| ld=| gp=| gn=| di=| re=| fr=",ORcol['Score'])
                        score = 0
                        tot=0
                        for i in list:                                
                                if len(re.split("/",i))>1:
                                        num = string.atof(re.split("/",i)[0])
                                        den = string.atof(re.split("/",i)[1])
                                        score += num/den
                                        tot += den
                                elif  len(re.split("-",i))>1:
                                        num = string.atof(re.split("-",i)[0])
                                        den = string.atof(re.split("-",i)[1])
                                        score += num-den
                                        tot += 1
                                elif i=='' or i == ' ':
                                        tot=tot
                                else:
                                        score += string.atof(i)
                                        tot += 1
                        sqlor= "SELECT * FROM SNPOR WHERE DiseaseID=%d AND PMID=%d AND dbsnp='%s' AND broadethnicity='%s'" % (disease,ORcol['PMID'],ORcol['dbSNP'],ORcol['BroadEthnicity'])
                        ngsign= ng=0
                        try:
                                cursor.execute(sqlor)
                                np,ngsign = signOR(np,ngsign,cursor)
                        except: print "Error: unable to fecth data %s" % sqlor
                        if ngsign <1 : score=0
                        sqlu = "UPDATE SNPOR set Score = replace(Score,'to= ','to=%.2f/%d ') WHERE DiseaseID=%d AND PMID= %d AND BroadEthnicity='%s' AND dbsnp='%s'" % (score,tot ,disease,ORcol['PMID'],ORcol['BroadEthnicity'],ORcol['dbSNP'])
                #print sqlu
                        try:
                                while cursor.execute(sqlu) >0 :
                                        cursor.execute(sqlu)
                        except: print "Error: unable to fecth data %s" % sqlu
except: print "Error:  to= calculation %s" % sqlsc 



# disconnect from server
db.close()


