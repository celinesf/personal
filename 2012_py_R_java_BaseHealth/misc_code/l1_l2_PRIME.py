#!/usr/bin/python
############ ALgo
################### ! for non dbsnp data check make_call for gt0 gt1 swuitch
# 285: look for position / check dbsnp name match too (accept position first
# 292 found dbSNP ID in chromosome file
# 293: found build37 position in DBANP_temp
# 294 multi alleles noted by ctg_hits 
# 295 record coverage
# 296 find snp by build37 position if not -> dbsnp/ + make sure single hit in genotype data
# 297 other snp sequencing + look for that SNP in genotyping data too
# ?? added position from DBSNP_temp in seq_position
# 298 'bcalls_used' in output for celine
# 299 output coverage in output for celine
# 301 Coverage <10 not good  + low statistics. changed ING position and genotype information output
# 302 if hetero: if any allele <5 <5 not good
# 303 if qsnp<40 and gt = ref alleles in make_call
# 305 if qsnp<40 and gt != ref alleles in make_call
# 308 309 qsnp>40 gt!=ref 
# 310 311 gt==ref
# 313 check_ref_sninfo -> ref allele is ok and dbsnp alleles as snpinfo
# 314 check_reference input values + test p of chec-ref_snpinfo
# 316 check_ref_sninfo test negative strand
# 317 test check_ref_sninfo multi alleles
# 318 NA for conflict genotyping and sequencing call + checked problem genotype individual error + genotyping tests
# 319 change check_genotyping: make sure same alleles as SNPInfo unless - strand
# 322 added tab in findit with position + added gt- unref... as stat tag that were missing
# 323 flag new alleles in gt call from SNPInfo + data test for reza
#######
####### FIXED nhit<=1-> merge case find genotype not second  snp (need find positon for it? plus once merged need to add in reseq
# ALL_EQ_REV- add AT CG count
# - 01/10 merged SNPs not counted as not found.
# - 01/10 mnew rsid from merged SNPs record the merged SNP? check both together???~
# - 01/11-0 SNPs with In_del column not empty need to be ignored 
# - 01/11-0 mutihit as missing data/ or take the closest position?
############### tested
# ---01/18/12 change -- to NA
# + no more position or chromosome on genotyping. only forward alleles 
# - accept the second genotype
# - check sequencing - changes:
# --- not found sequencing, in genotyping but reverse strand : rs2187668 chr 1 _ check snpinfo
# --- Genotyping strand- + sequencing strand- :rs646776 chr 1 get genotyping
# --- in dbsnp take the second genotype NOW only always but show if error (make call)
# --- check ref
# - do not give dbsnp too much weight.
# - the sign are only for dbsnp so change function somehow
# - count any conflict with genotyping
#
# speed0:
# loop on file
# -> functions in separate file
# check for time
# output less and in table form only
# 3-> check_ref (should be called concordance?)
# 4 -> make call from different strand + change strand, calculate stat+ edit problem 
# s1 -> reza file with rsid chromosome genotype only
# s2 -> check same alleles/genotype in + strand
# s3 -> check with genotyping
# s4 -> check fill not found SNP with genotyping data
# s5 -> check that XY and M are all homozygous
# s6 check genotyping and genotype call are the same
# 11/16 -> find genotype from position and chr number
#       -> get alleles from ins-del when minor/maor are empty
import re
import sys
import time
import function #@UnresolvedImport
import MySQLdb #@UnresolvedImport
import datetime

def main(DIR,IND):
    now = datetime.datetime.now()
    ############## CUSTOMIZE per Member? ##########
    dirdbsnp= '%s/%s/Variations/dbSNP/' % (DIR,IND)
    dirgeno= '%s/%s/Genotyping/' % (DIR,IND)
    fngt= '%sFinalReport_HumanOmni2.5-8v1_%s.txt' % (dirgeno,IND)
    startf = 'chr'
    endf = '_custom_dbSNP.txt'
    fnout = '%s_dbSNPs_%s'  % (IND,"z")
    #fnout = '%s_dbSNPs_%s' % (IND,now.strftime("%Y-%m-%d_%H-%M"))
    ########################
    sys.stdout =open('%s_OUT' % fnout, 'w')
    ####################### what is the sex of the member
    sex = 0
    chrfile = '%s%s%s' % (startf, 'Y', endf)
    fnin = '%s%s' % (dirdbsnp,chrfile)
    try:
        fileIN =    open(fnin)
        fileIN.close()
        sex = 1
        print "%s is MALE" % IND
    except :
        print "%s is FEMALE" % IND
        sex = 0
    #################### connection to SNPInfo -> get SNPInfo list of dbSNPs ID
    db = MySQLdb.connect(host="192.168.2.104",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    rssnpinfo=[] # rsid in SNPInfo
    sql0 = "SELECT distinct chromosome from SNPInfo "
    #where  (dbsnp= 'rs12743918' or dbsnp= 'rs71206208'or dbsnp= 'rs2229276' or  dbsnp='rs16834521')
    # or dbsnp='g1182C>T' or dbsnp='rs1042164' or dbsnp='rs11084753' or dbsnp='rs12049330' or dbsnp='rs12748299' or dbsnp='rs1353456' or dbsnp='rs150160' or dbsnp='rs2038366' or dbsnp='rs2187668' or dbsnp='rs2228564' or dbsnp='rs2273289' or dbsnp='rs2476601' or dbsnp='rs2524005' or dbsnp='rs2568958' or dbsnp='rs2580520' or dbsnp='rs259940' or dbsnp='rs3008621' or dbsnp='rs3737964' or dbsnp='rs4712524' or dbsnp='rs5916245' or dbsnp='rs5944185' or dbsnp='rs5970164' or dbsnp='rs61886492' or dbsnp='rs646776' or dbsnp='rs673' or dbsnp='rs7533552' or dbsnp='rs8321' or dbsnp='rs898165' or dbsnp='rs9268858')"
    try:
        cursor.execute(sql0)
        tp = cursor.fetchall()
        desc = cursor.description
        chrlist={}
        for s in tp:
            if s[0] != None and s[0] != '':
                chrlist[s[0]] = s[0]
        print 'list of chromosome names:', chrlist
    except: print "Error: unable to fecth data %s" % sql0
    ################
    chrlist={}
    chrlist['6']='6'
    ### Init clock
    t0=time.clock()
    ### Init Output files: data for systen
    fileING = open('%s.txt' % fnout,'w')
    hing=['dbsnp','chr','pos','gt']
    ## output for QC
    fileOUT = open(fnout,'w')
    hout=['dbsnp','chr','pos','ref','strand','alleles','gt','Q(gt)','Q(snp)','bcalls_used','base','sequencing','SNPInfo','genotyping']
    ## statistics summaries
    fileSTAT = open('%s_STAT' % fnout,'w')
    hstat=['File','nSNP','time','low','lhet','high','found','Mhit','mhit','nhit','un','homMYX','hetMYX',\
           'notfound','nfOK','nfm','nfmnf','nfgm','findel',\
           'nopos','nfpos','nfs','copymerged',\
           'rok','rpb',\
           'psf','psc','psd','prf','prd','nsf','nsc','nsd','nrf','nrd','Mpsame','Mprev','Mpdif','Mnsame','Mnrev','Mndif',\
           'genodata','geno','gs+','gs-','gr+','gr-','gnog','gmiss','gind',\
           'gna','gsame','gdif',\
           'OK+','OK-','g+','g-','gq+','gq-','strandg+','strandg-','strandr-','strandr+','rev+','rev-',\
           'ARG_NO','G_NOAR','AR_A+','AR_A-','AR_RG+','AR_RG-','ARGNO',
           'gt-','unall','unref','gfill',\
           'reh0f','reh0d','reh1f','reh1d','rel0f','rel0d','rel1f','rel1d','rnh0f','rnh0d','rnh1f','rnh1d','rnl0f','rnl0d','rnl1f','rnl1d',\
           'snh0f','snh0d','snh1f','snh1d','snl0f','snl0d','snl1f','snl1d','seh0f','seh0d','seh1f','seh1d','sel0f','sel0d','sel1f','sel1d',]
    totstat = {} # total for all sequencing data       
    for h in range(0,max(len(hout),len(hstat),len(hing))):
        if h<len(hing):
            fileING.write('%s\t'%hing[h])
        if h<len(hstat):
            fileSTAT.write('%s\t'%hstat[h])
            totstat[hstat[h]] = 0
        if h<len(hout):
            fileOUT.write('%s\t'%hout[h])
    fileOUT.write("\n")
    fileING.write("\n")
    fileOUT.flush()
    fileING.flush()
    fileSTAT.flush()
    ############### Recover the genotyping data
    print '----------------->\tStarting to read %s' % fngt
    sys.stdout.flush()   
    (numrsgt,rsgt,gtdata) = function.get_genotyping_data(fngt)
    t1=time.clock()
    print '\tI got %d SNPs in \'%s\' (%ss).' % (numrsgt-1,fngt,t1-t0)
    sys.stdout.flush()   
    fileSTAT.write('\n%s\t%d\t%s\n' % (fngt,numrsgt-1,t1-t0)) # stats on genotyping file
    fileSTAT.flush()
    ### file summary for each CHR sequencing data
    ########### LOOP on CHROMOSOME sequencing data
    rsseq=[] # rsid in sequencing data
    for CHR in chrlist:
        t1=time.clock()
        chrfile = '%s%s%s' % (startf, CHR, endf)
        print  '----------------->\tI am in file', chrfile
        ## Init summary on chr
        statchr={}        
        for h in hstat: statchr[h] = 0
        statchr['File']=chrfile
        fnin = '%s%s' % (dirdbsnp,chrfile)
        fileok=0
        try:
            fileIN =  open(fnin)
            fileok = 1
        except :
                print 'NOFILE: I could not find the file \'%s\'' % (fnin)
        if fileok ==1 :
            dbdata = fileIN.read()
            sql1 = "SELECT * from SNPInfo where  chromosome='%s' and dbsnp='rs2187668'" % CHR
            #sql1 = "SELECT * from SNPInfo where  chromosome='%s' and (dbsnp='rs16834521' or dbsnp='rs2229276')" % CHR
            #sql1 = "SELECT * from SNPInfo where chromosome='%s' and (dbsnp= 'rs16948048' or dbsnp= 'rs1042522'or dbsnp= 'rs1671021' or  dbsnp='rs3737964' or dbsnp='g1182C>T' or dbsnp='rs1042164' or dbsnp='rs11084753' or dbsnp='rs12049330' or dbsnp='rs12748299' or dbsnp='rs1353456' or dbsnp='rs150160' or dbsnp='rs2038366' or dbsnp='rs2187668' or dbsnp='rs2228564' or dbsnp='rs2273289' or dbsnp='rs2476601' or dbsnp='rs2524005' or dbsnp='rs2568958' or dbsnp='rs2580520' or dbsnp='rs259940' or dbsnp='rs3008621' or dbsnp='rs3737964' or dbsnp='rs4712524' or dbsnp='rs5916245' or dbsnp='rs5944185' or dbsnp='rs5970164' or dbsnp='rs61886492' or dbsnp='rs646776' or dbsnp='rs673' or dbsnp='rs7533552' or dbsnp='rs8321' or dbsnp='rs898165' or dbsnp='rs9268858')" % CHR
            SNPInfo={}
            try:
                cursor.execute(sql1)
                tab = cursor.fetchall()
                desc = cursor.description
                for s in tab:
                    rssnpinfo.append(s[0].lower())
                    SNPInfo[s[0].lower()]={}
                    for (name, value) in zip(desc, s) :
                        SNPInfo[s[0].lower()][name[0]] = value
            except: print "Error: unable to fecth data %s" % sql1
            print '# SNP in genotyping %d. Total # SNPs: %d' % (len(set(SNPInfo).intersection(set(rsgt))), len(SNPInfo))
            #a=0
            for snp in SNPInfo:
                #if a==1: snp='rs2229276'
                #else: snp= 'rs16834521'
                #a+=1
                print snp
                build37,statchr,foundpos, foundsnp,found37 = function.get_build37(cursor ,CHR,snp,SNPInfo,statchr,dbdata,fnin)
                statchr['nSNP'] += 1
                if  len(re.split('[,|_| ]',snp)) >1 or (SNPInfo[snp]['MinorAllele'] == '' or SNPInfo[snp]['MajorAllele'] == ''  or SNPInfo[snp]['dbSNP'] == ''  or SNPInfo[snp]['Position'] == '' or SNPInfo[snp]['MinorAllele'] == None or SNPInfo[snp]['MajorAllele'] == None  or SNPInfo[snp]['dbSNP'] == None  or SNPInfo[snp]['Position'] == None ) and (SNPInfo[snp]['Ins_Del'] == '' or SNPInfo[snp]['Ins_Del'] ==  None):
                    print 'PROBLEM: I found empty data or wrong ID for dbsnp %s in SNPInfo: chr %s pos %s maj %s min %s, INDEL %s ' % (SNPInfo[snp]['dbSNP'], SNPInfo[snp]['Chromosome'],SNPInfo[snp]['Position'], SNPInfo[snp]['MajorAllele'] ,SNPInfo[snp]['MinorAllele'],SNPInfo[snp]['Ins_Del'])
                elif rsseq.count(snp) == 0:
                    ########## SNP not found in sequencing  ################
                    snp2 = snp 
                    com_indel = 0
                    nfindel='NA'
                    if found37 == 0: nfindel = function.fill_problem(nfindel,"37" ) # record no position 37 in algophen
                    if foundsnp == None and foundpos == None  and rsseq.count(snp) == 0:
                        ########## merged SNPs : new SNP should be in
                        statchr['notfound'] += 1
                        nfindel = function.fill_problem(nfindel,"NF" )
                        if  SNPInfo[snp]['Comments'] != None: 
                            if SNPInfo[snp]['Comments'].strip().split()[0] == '!new' or SNPInfo[snp]['Comments'].strip().split()[0]  == '!original': # merged SNPs
                                com_indel = 1
                                snp2 = SNPInfo[snp]['Comments'].strip().split()[2]
                                if re.search('\t%s\t' % snp2.lower(), dbdata.lower()) == None and rsseq.count(snp2) == 0:
                                    statchr['notfound'] += 1
                                    statchr['nfmnf'] += 2
                                    nfindel = function.fill_problem(nfindel,"MNF" )
                                    print 'NOTFOUNDM: I could not find dbSNP merged %s and new %s in file \'%s\'. Comments %s, Indel %s' % (snp,snp2, fnin,SNPInfo[snp]['Comments'],SNPInfo[snp]['Ins_Del'])
                                elif rsseq.count(snp2) == 0:   ##### Call New RSID - same as old RSID
                                    statchr['nfm'] += 1
                                    nfindel = function.fill_problem(nfindel,"M" )
                                    print 'NF_M: I could not find dbSNP merged %s (new/old= %s)in file \'%s\'. Comments %s, Indel %s' % (snp,snp2, fnin,SNPInfo[snp]['Comments'],SNPInfo[snp]['Ins_Del'])   
                                    build372,statchr,foundpos, foundsnp,found372 = function.get_build37(cursor ,CHR,snp2,SNPInfo,statchr,dbdata,fnin)
                                    for i in range(0,len(build372)):
                                        print snp2, snp
                                        (statchr, rsseq) = function.sequencing_call(IND, sex,CHR,snp2,snp,build372[i],foundpos,found372,dbdata,statchr,rsseq,SNPInfo,gtdata,rsgt,fngt,fnin,fileOUT,fileING, hout, hing,hstat,nfindel)
                            elif SNPInfo[snp]['Ins_Del'] != None and SNPInfo[snp]['Ins_Del'] != '':
                                statchr['nfOK'] += 1
                                com_indel = 2
                                nfindel = function.fill_problem(nfindel,"OK1" )
                                print 'NF_OK1: I could not find dbSNP %s in file \'%s\'. Comments %s, Indel %s' % (snp, fnin,  SNPInfo[snp]['Comments'],SNPInfo[snp]['Ins_Del'] )
                            else:
                                nfindel = function.fill_problem(nfindel,"nf1" )
                                print 'NOTFOUND1 I could not find dbSNP %s in file \'%s\'. Comments %s, Indel %s' % (snp, fnin ,  SNPInfo[snp]['Comments'],  SNPInfo[snp]['Ins_Del'] )   
                        elif SNPInfo[snp]['Ins_Del'] != None and SNPInfo[snp]['Ins_Del'] != '':
                            statchr['nfOK'] += 1
                            com_indel = 2
                            nfindel = function.fill_problem(nfindel,"OK2" )
                            print 'NF_OK2: I could not find dbSNP %s in file \'%s\'. Comments %s, Indel %s' % (snp, fnin, SNPInfo[snp]['Comments'], SNPInfo[snp]['Ins_Del'] )
                        else:
                            nfindel = function.fill_problem(nfindel,"nf2" )
                            print 'NOTFOUND2 I could not find dbSNP %s in file \'%s\'. Comments %s, Indel %s' % (snp, fnin ,  SNPInfo[snp]['Comments'],  SNPInfo[snp]['Ins_Del'] )
                        if rsseq.count(snp) == 0: ## check if can find in genotyping?
                            (seqdata,statchr, rsgt) = function.get_genotyping_call(gtdata,rsgt,snp,IND,CHR,statchr,fngt,hstat, nfindel)
                            seqdata['dbsnp2'] = snp2
                            if snp !=  snp2 and rsseq.count(snp2) == 0 and seqdata['gt'][len(seqdata['gt'])-1] == 'NA': 
                                (seqdata,statchr, rsgt) = function.get_genotyping_call(gtdata,rsgt,snp2,IND,CHR,statchr,fngt,hstat, nfindel)
                                tpsnp=snp2
                                snp2=snp
                                snp=tpsnp
                                seqdata['dbsnp2'] = snp2
                            if com_indel == 2: 
                                seqdata['sequencing'] = function.fill_problem(seqdata['sequencing'],"INDEL")
                            elif com_indel == 1: 
                                seqdata['sequencing'] = function.fill_problem(seqdata['sequencing'],"MERGE")
                            if  seqdata['genotyping'] == 'DATA': # checked
                                rsseq.append(snp)
                                (seqdata,statchr) = function.check_ref_dbsnp_SNPInfo(SNPInfo[snp],seqdata,statchr)
                                (seqdata,statchr, rsgt) = function.check_genotyping(IND,seqdata,gtdata,statchr,rsgt,fngt)
                                (seqdata,statchr) = function.genotype_call(snp,seqdata,statchr,hstat)
                                #seqdata['genotyping'] = function.fill_problem(seqdata['genotyping'],"NF_FGENO")
                                #print('FOUNDGENO\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s'%(seqdata['dbsnp'],seqdata['chr'],seqdata['ref'][0],seqdata['strand'],seqdata['alleles'],seqdata['gt']))#check
                            function.write_output(hout,seqdata,fileOUT,0) # write stats for chr
                            function.write_output(hing,seqdata,fileING,1) # write stats for chr
                            if seqdata['dbsnp'] !=  seqdata['dbsnp2'] and rsseq.count(seqdata['dbsnp2']) == 0: 
                                seqdata2,statchr,rsseq = function.merge(seqdata['dbsnp'],statchr,rsseq,seqdata,hstat)
                                function.write_output(hout,seqdata2,fileOUT,0) # write stats for chr
                                function.write_output(hing,seqdata2,fileING,1) # write stats for chr
                    else: ######### call found in sequencing data
                        if  SNPInfo[snp]['Comments'] != None and SNPInfo[snp]['Comments'] != '': 
                            if SNPInfo[snp]['Comments'].strip().split()[0] == '!original' or SNPInfo[snp]['Comments'].strip().split()[0] == '!new': # merged SNPs
                                snp2 = SNPInfo[snp]['Comments'].strip().split()[2]
                                if re.search('\t%s\t' % snp2.lower(), dbdata.lower()) == None:
                                    statchr['notfound'] += 1
                                    statchr['nfm'] += 1
                                    nfindel = function.fill_problem(nfindel,"NF_M2" )
                                    print 'NF_M2: dbSNP merge %s (old new %s) in file \'%s\'. Comments %s, Indel %s' % (snp2, snp, fnin,  SNPInfo[snp]['Comments'], SNPInfo[snp]['Ins_Del'] )
                        if SNPInfo[snp]['Ins_Del'] != None and SNPInfo[snp]['Ins_Del'] != '':
                            statchr['findel'] += 1
                            nfindel = function.fill_problem(nfindel,"INDEL" )
                            print 'INDEL_FOUND: dbSNP %s in file \'%s\'. Comments %s, Indel %s' % (snp, fnin,  SNPInfo[snp]['Comments'], SNPInfo[snp]['Ins_Del'] )
                        if rsseq.count(snp) == 0:
                            if found37 == 0:
                                (statchr, rsseq) = function.sequencing_call(IND,sex, CHR,snp,snp2,build37,foundpos,found37,dbdata,statchr,rsseq,SNPInfo,gtdata,rsgt,fngt,fnin,fileOUT,fileING, hout, hing,hstat,nfindel)  
                            for i in range(0,len(build37)):
                                (statchr, rsseq) = function.sequencing_call(IND,sex, CHR,snp,snp2,build37[i],foundpos,found37,dbdata,statchr,rsseq,SNPInfo,gtdata,rsgt,fngt,fnin,fileOUT,fileING, hout, hing,hstat,nfindel)  
                sys.stdout.flush()                      
            ## END loop on SNP
            print '----------------->\tnow I have # SNP in genotyping %d' % (len(set(SNPInfo).intersection(set(rsgt))))
            fileIN.close()
            t2=time.clock()
            statchr['time'] =  t2-t1
            function.write_stat(hstat,statchr,fileSTAT) # write stats for chr
            totstat = function.add_stat(hstat,statchr,totstat) # add total stat
            fileOUT.flush()
            sys.stdout.flush()
    print rsseq
    ## END LOOP on CHROMOSOME sequencing data
    fileING.close()
    fileOUT.close()
    ### Final Statistics
    t2=time.clock()
    totstat['time'] = t2-t0
    totstat['File'] = 'Total'
    function.write_stat(hstat,totstat,fileSTAT)
    fileSTAT.close()
    print '# SNP I found in genotyping %d. Total # SNPs in dbSNPs: %d' % (len(set(rssnpinfo).intersection(set(rsgt))), len(rssnpinfo))
    print '# SNPs with no data found: %d' % len(set(rssnpinfo).difference(set(rsseq)))
    sys.stdout.flush()


