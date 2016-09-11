#!/usr/bin/python
############ Function modules
# compare_genotype
# remove_slash
# compare_allele
# reverse_strand
# genotype_call
# get_sequencing_data 
# add_stat
# write_stat
# get_genotyping_data

import re
import string
import linecache

#######################
def sequencing_call(IND,sex, CHR,rsid,rsid2,dbsnpdata,statchr,rsseq,SNPInfo,build37,gtdata,rsgt,fngt,fnin,fOUT,fING, hout, hing,hstat,nfindel):
    statchr['found'] += 1
    nhit = build37[rsid]['total_hits']
    
    lmin = 1e9
    #for m in re.finditer('\t%s\t' % rsid.lower(), dbsnpdata.lower()):
    #    nhit += 1
    if nhit>1:
        print "MHIT"
        statchr['mhit'] += 1
        statchr['nhit'] += nhit
        for m in re.finditer('\t%s\t' % rsid.lower(), dbsnpdata.lower()):
            start = m.start()
            lineno = dbsnpdata.count('\n', 0, start) + 1
            word = m.group(0)
            line = linecache.getline(fnin, lineno)
            word =(line.strip().split())   
            if abs(string.atoi(word[1]) - SNPInfo[rsid]['Position']) < lmin:
                lmin = abs(string.atoi(word[1]) - SNPInfo[rsid]['Position'])
    ########### Chose the genotype for the closest to SNPInfo?
    nrsid = 0 # hit number
    for m in re.finditer('\t%s\t' % build37[rsid]['chr_pos'], dbsnpdata.lower()):
        nrsid += 1
        start = m.start()
        lineno = dbsnpdata.count('\n', 0, start) + 1
        word = m.group(0)
        line = linecache.getline(fnin, lineno)
        word =(line.strip().split())      
        ## read dbSNP resequencing data/ genotype call
        if word[0] == ('chr%s' % CHR) and nhit>=1: 
            #rsid = word[10+4].lower()
            rsseq.append(rsid)
            (seq) = get_sequencing_data(word,rsid,CHR,statchr,nfindel)
            print 'A',seq
            ## check the comments and idel of SNPINFO
            (seq,statchr) = check_SNPInfo(SNPInfo[rsid],seq,statchr)
            seq['pos'].append(build37[rsid]['chr_pos'])
            print 'B',seq
            if nhit > 1:
                seq['sequencing'] = fill_problem(seq['sequencing'],"MHIT_%s_%s" % (seq['gt'][0],seq['gt'][1]) )
                seq['gt'].append('NA')
                print('MULTI_HIT\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s'  %(seq['dbsnp'],seq['chr'],seq['ref'][0],seq['strand'],seq['alleles'],seq['gt']))
            print 'C'
            if  SNPInfo[rsid]['Comments'] != None and SNPInfo[rsid]['Comments'] != '' and SNPInfo[rsid]['Comments'].strip().split()[0] == "MultiHit": # mutiple hit has in SNPIfo
                seq['SNPInfo'] = fill_problem(seq['SNPInfo'],"mhit" )
                if nrsid == 1:
                        statchr['Mhit'] +=  1
                seq['gt'].append('NA')
                nhit += 1
                print('mhit\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s'  %(seq['dbsnp'],seq['chr'],seq['ref'][0],seq['strand'],seq['alleles'],seq['gt']))#check
            print 'D'
            if nhit == 1 :   
                seq['dbsnp2'] = rsid2
                (seq,statchr, rsgt) = check_genotyping(IND,seq,gtdata,statchr,rsgt,fngt)
                (seq,statchr) = genotype_call(rsid,seq,statchr,hstat)
            print 'e'
            ### Check Genotype and genotyping are the same
            if nhit == 1 and len(seq['gtg']) > 1 and  seq['genotyping'] != 'NA':
                revgeno = reverse_strand(seq['gtg'][0], rsid,CHR)
                same = compare_genotype(seq['gt'][len(seq['gt'])-1],seq['gtg'][0])
                samerev = compare_genotype(seq['gt'][len(seq['gt'])-1],revgeno)
                if same == 0 and samerev == 0:
                    if seq['gt'][len(seq['gt'])-1] != 'NA' and seq['gtg'][0] != 'NA' :
                        statchr['gdif'] +=1
                        print seq
                        print "DIF dbSNP %s, the call is %s (score %s) while genotyping is %s (score %s)" %(rsid,seq['gt'][len(seq['gt'])-1],seq['Q(gt)'][len(seq['Q(gt)'])-1],seq['gtg'][0],seq['gtg'][1]  )
                else:  
                    statchr['gsame'] +=1
            print 'f'
            ### Check M YX are homozygous
            if nhit == 1 and  (CHR == 'M' or CHR == 'Y' or  CHR == 'X'):
                if seq['gt'][len(seq['gt'])-1][0] != seq['gt'][len(seq['gt'])-1][1]:
                    statchr['hetMYX'] +=1
                    if sex==1 and  (CHR == 'Y' or  CHR == 'X'):
                        print "HET for a MALE or mtDNA I found an heterozygous error dbSNP %s, genotype %s" %(rsid,seq['gt'][len(seq['gt'])-1]) 
                else:
                    if sex==1 and  (CHR == 'Y' or  CHR == 'X'):
                        seq['gt'][len(seq['gt'])-1]=seq['gt'][len(seq['gt'])-1][0]
                    statchr['homMYX'] +=1
            print 'g'
            seq['sequencing'] = fill_problem(seq['sequencing'],seq['Q(snp)'][0])
            print 'h'
            write_output(hout,seq,fOUT,0) # write stats for chr
            if  nrsid ==1:             
                write_output(hing,seq,fING,1) # write stats for chr
            print 'i'
            if seq['dbsnp'] !=  seq['dbsnp2'] and rsseq.count(seq['dbsnp2']) == 0: # merge old into new
                seq2,statchr,rsseq = merge(rsid,statchr,rsseq,seq,hstat)
                write_output(hout,seq2,fOUT,0) # write res
                write_output(hing,seq2,fING,1) # write ING final file
            print 'j'
        elif word[0] != 'seq_name':
            if nhit >1: print 'MHIT: I found  (Ignored) \'%s\' in file %s' % (word,fnin)
            else:       print 'CHR_SEQ: I found \'%s\' in file %s' % (word,fnin)
    return statchr, rsseq
############### call from sequencing

################
def merge(snp,stat,rslist,seq1,stath):
    stat['copymerged'] += 1
    if rslist.count(seq1['dbsnp2']) == 0: rslist.append(seq1['dbsnp2'])
    if rslist.count(seq1['dbsnp']) == 0: rslist.append(seq1['dbsnp'])
    seq2 = seq1
    seq2['dbsnp'] = seq1['dbsnp2']
    seq2['dbsnp2'] = snp
    for h in stath: # add this merged SNP in the statistics - so that numbers fit with summary table
        if h != 'File':
            stat[h] += seq2['OK'][h]
    seq2['sequencing'] = fill_problem(seq2['sequencing'],"M_%s" %  seq2['dbsnp2']) 
    print('MERGED\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(seq2['dbsnp'],seq2['dbsnp2'],seq2['chr'],seq2['ref'][0],seq2['strand'],seq2['alleles'],seq2['gt']))#check
    return seq2, stat, rslist
################### end merge

##########
def get_build37(curs, chro,rsid,snpinfo,stat,data,fin):
    sql2 = "SELECT * from DBSNP_temp where  chr='%s' and genome='GRCh37.p5' and rs=%s" % (chro,rsid.split('rs')[1])
    build37={}
    print sql2
    try:
        curs.execute(sql2)
        tab = curs.fetchall()
        desc = curs.description
        for s in tab:
            print s
            build37['rs%d' % s[0]]={}
            for (name, value) in zip(desc, s) :
                build37['rs%d' % s[0]][name[0]] = value 
            if build37[rsid]['total_hits'] > 1:
                
                break
    except: print "Error: unable to fecth data %s" % sql2
    foundpos = None
    try: 
        foundpos = re.search('\t%s\t' % build37[rsid]['chr_pos'], data.lower())
    except :
        stat['nfpos'] += 1
        print 'NF_POS: I could not find dbSNP %s in file \'%s\'. Comments %s, Indel %s' % (rsid, fin,  snpinfo[rsid]['Comments'],snpinfo[rsid]['Ins_Del'] )
    return(build37,stat,foundpos)
        
################
def get_genotyping_call(genotyping,genodb,rsid,indv,chrnum,stat,fgeno,hstat,nfindel):
    sequencing = {}
    sequencing['dbsnp'] = rsid
    sequencing['dbsnp2'] = rsid
    sequencing['chr'] =chrnum
    sequencing['pos'] = []
    sequencing['pos'].append(-1) # position in 0 Sequencing/Genotyping, 1 in dbsnp
    sequencing['OK'] = {} # record all stat added for that SNP
    for h in hstat: sequencing['OK'][h] = 0
    sequencing['ref'] = []
    sequencing['ref'].append('-') # major/ minor alleles from SNPInfo
    sequencing['Q(snp)'] = []
    sequencing['Q(snp)'].append(0) # quality SNP given in sequencing file
    sequencing['Q(ref)'] = [] # comp to alleles
    sequencing['Q(ref)'].append(0)  
    sequencing['base'] =[] # base_used
    sequencing['alleles'] =[]
    sequencing['alleles'].append('NA')
    sequencing['strand'] = []
    sequencing['strand'].append('+')
    sequencing['R(gt)']=[] # comp to reference gt1/+ gt1/- gt2/+ gt2/-  same in revgt1...
    sequencing['gt']=[] # all genotype calls : 
    sequencing['gt'].append('NA') # genotype in sequencing 0 max_gt,1 poly+gt 2 rev gt 3 rev gt 4 = genotype CALL 
    sequencing['gt'].append('NA')
    sequencing['Q(gt)']=[]
    # Q in sequencing 0 Q max_gt,1 Q poly+gt, 2 comp allele/gt 3 comp allele/pgt, 4 & 5 comp alleles/rev mgt or pgt,
    # 6 are the two GT the same in sequencing  7 GT call/total score from illumina
    sequencing['Q(gt)'].append(0)
    sequencing['Q(gt)'].append(0)
    sequencing['sequencing']=nfindel
    sequencing['genotyping']='NA'
    sequencing['SNPInfo']='NA'
    sequencing['gti'] = ''
    sequencing['gtg'] = [] # genotype in genotyping
    okkey,gsnp, genodb = find_genotyping(sequencing,genotyping,genodb,fgeno)
    if okkey ==0:
        if genotyping[gsnp][1].find(indv) >-1 :#genotyping[gsnp][1] == indv:
            if chrnum == chrnum: #((genotyping[gsnp][2][0] == chrnum or genotyping[gsnp][2] == chrnum) ):
                print('GENOFOUND: in \'%s\' I found genotyping data for \'%s\' (called \'%s\').' % (fgeno,rsid,gsnp))# check
                gt = '%s%s' %(genotyping[gsnp][2],genotyping[gsnp][3]) #(genotyping[gsnp][6],genotyping[gsnp][7])
                sequencing['ref'][0] = genotyping[gsnp][2]
                sequencing['alleles'][0] = '%s%s' % (genotyping[gsnp][2],genotyping[gsnp][3]) #!!! fake allele known ( re.split('[\[\]]',genotyping[gsnp][4])[1])
                sequencing['gt'][0] = gt
                sequencing['gt'][1] = gt
                sequencing['gtg'].append(gt)
                sequencing['gtg'].append( genotyping[gsnp][6])
                sequencing['genotyping']='DATA'
                stat['genodata'] +=1
                if  sequencing['chr'] !=chrnum:
                    print sequencing
                    print('CHROMOSOME\tSEQ\t%s\t%s\t%s' % (gsnp, chrnum,sequencing['chr'])) #check
    return sequencing, stat, genodb
########### end .get_sequencing_data

##################
def find_genotyping(seqdata, genodata,dbg,fg):
    ok=0
    gtid = seqdata['dbsnp'] # find the rs number in genotype data
    try:
        genodata[gtid]
    except: 
        ok = 1
        if seqdata['dbsnp'] != seqdata['dbsnp2'] : # try to find the merged SNP
            try:
                genodata[seqdata['dbsnp2']]
                ok = 2
            except: ok = 1
        ########### no position or chromosome number in HumanOmni2.5-8v1_A.bpm ###
        #if ok == 1 : # try to find the position in the chromosome
        #       for g in genodata:
        #              if  string.atoi(genodata[g][3]) == seqdata['pos'] and (genodata[g][2][0] == seqdata['chr'] or genodata[g][2] == seqdata['chr']):                                               
        #                     dbg[dbg.index(g)]=gtid
        #                    fPROB.write('GENOPOS_NOTdbSNP: in \'%s\' I found dbSNP %s (chr %s pos %d) named %s on chr \'%s\' (pos %s).' % (fg,gtid, seqdata['chr'], seqdata['pos'], genodata[g][2],genodata[g][3]))# check
        #                   ok = 0
        #                  gtid = genodata[g][0]
        #                 print 'Rescue data with position', seqdata['dbsnp'], gtid,genodata[g]
    return ok,gtid, dbg
################ end  find_genotyping

###################
def check_genotyping(indv,sequencing,genotyping,stat,dbgeno,fgeno):
    okkey,gsnp, dbgeno = find_genotyping(sequencing,genotyping,dbgeno,fgeno)
    if okkey == 2:
        gsnp = sequencing['dbsnp2']
        stat['nfgm']+=1
        okkey = 0
        sequencing['genotyping'] = fill_problem(sequencing['genotyping'],"GM_%s" % gsnp)
        print('GMERGE\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(sequencing['dbsnp'],sequencing['dbsnp2'],sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt']))#check
    if okkey == 0:
        ################### no more chr or position
        if genotyping[gsnp][1].find(indv) >-1 :
            if indv == indv:
                print len(sequencing['alleles'])-1
                gt = sequencing['alleles'][len(sequencing['alleles'])-1]
                gtrev = reverse_strand(gt,sequencing['dbsnp'],sequencing['chr'])
                geno = '%s%s'% (genotyping[gsnp][2],genotyping[gsnp][3])
                genorev = reverse_strand(geno,sequencing['dbsnp'],sequencing['chr'])
                stat['geno'] +=1
                sequencing['gtg'].append(geno) # genotyping
                sequencing['gtg'].append(genotyping[gsnp][6]) # score
                sameall = compare_alleles(geno,gt)
                sameallrev = compare_alleles(geno,gtrev)
                stat['gall'] +=1
                if sameall >=2:
                    stat['ggt+'] += 1 # checked
                    sequencing['OK']['ggt+'] += 1
                    sequencing['genotyping'] = fill_problem(sequencing['genotyping'],"ggt+_%s" % geno)
                    #fPROB.write('G_GT+\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s'%(sequencing['dbsnp'],sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt']))#check
                elif sameallrev >= 2:  # checked
                    stat['ggt-'] += 1
                    sequencing['OK']['ggt-'] += 1
                    sequencing['genotyping'] = fill_problem(sequencing['genotyping'],"ggt-_%s" % genorev)
                    sequencing['gtg'][0] = genorev
                    #fPROB.write('G_GT-\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s'%(sequencing['dbsnp'],sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt']))#check
                else:
                    if geno == 'NA':
                        stat['gmiss'] += 1 ###################### NONE?
                        sequencing['OK']['gmiss'] += 1
                        sequencing['genotyping'] = fill_problem(sequencing['genotyping'],"gmiss_%s" % geno) 
                        #fPROB.write
                        print('G_MISS\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s'%(sequencing['dbsnp'],sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt']))#check
                    else:   
                        stat['gnog'] += 1  ###################### NONE?
                        sequencing['OK']['gnog'] += 1
                        sequencing['genotyping'] = fill_problem(sequencing['genotyping'],'PRBALL_%s' % (geno)) 
                        print sequencing, genotyping[gsnp]
                        print('G_ALL\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s\tgeno: %s'%(sequencing['dbsnp'],sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'],geno))#check        
        else:  ###################### NONE?
            print 'd'
            stat['gind'] += 1
            sequencing['OK']['gind'] += 1
            sequencing['genotyping']=fill_problem(sequencing['genotyping'],'PBLIND_%s_%s' % (genotyping[gsnp][1],indv))
            print('G_IND: in \'%s\' I found data for \'%s\' when I am looking for %s' % (fgeno,genotyping[gsnp][1],indv)) #che
    return sequencing, stat, dbgeno
############### end check_genotyping

####################
def get_alleles_SNPInfo(info):
    if (info['MajorAllele'] == None or info['MajorAllele'] == '' or info['MinorAllele'] == None or info['MinorAllele'] == '') == False:
        gt = '%s/%s' % (info['MajorAllele'], info['MinorAllele'])
    elif info['Ins_Del'] != None and info['Ins_Del'] != '':
        gt = info['Ins_Del']
        ###################### NONE?
    else : print('NOALLELES: SNPInfo has no data at \'%s\': maj %s min %s indel %s' % (info['dbSNP'],info['MajorAllele'],info['MinorAllele'],info['Ins_Del']))
    return gt
################### end get_alleles_SNPInfo

################
### check correlation between SNPInfo and dbSNP alleles (from Sequencing ddata)
def check_SNPInfo(snpinfo,data,stat):
    data['pos'].append(snpinfo['Position'])
    gt = get_alleles_SNPInfo(snpinfo) # alleles M/m as in SNPInfo: + strand
    data['alleles'].append(gt)
    data['gti'] = gt #possible alleles in SNPINFO
    if  data['strand'][0] == '+': #  db alleles+ == SNPInfo alleles
        stat['db+'] +=1
        data['OK']['db+'] += 1
        data['SNPInfo'] = fill_problem(data['SNPInfo'],'DB+')
    else: #  db alleles+ == SNPInfo alleles
        stat['db-'] +=1
        data['OK']['db-'] += 1
        data['SNPInfo'] = fill_problem(data['SNPInfo'],'DB-')
        print('SHOULDNEVERBEHERE_DB-\tSNPINFO\t%s\t%s\t%s\t%s\t%s\t%s\tdbgt=%s'%(data['dbsnp'],data['chr'],data['ref'][0],data['strand'],data['alleles'],data['gt']))   
    return data,stat
################end check_spninfo
        
#################
def make_call(rs,sequencing, stat,hstat):
    tot = sequencing['Q(gt)'][0] + sequencing['Q(gt)'][1]
    ### case Allele or ref
    qt = [sequencing['Q(gt)'][i] for i in range(2,4)] # quality for both genotype]
    ref = sequencing['OK']['AR_RG+']+sequencing['OK']['AR_RG-']
    if ref >0:   qt = [sequencing['R(gt)'][i] for i in range(0,2)] # fit with reference + strand
    ### dbsnp reference strand
    strand=  sequencing['strand'][0]
    if len(sequencing['strand']) > 1:
            strand = sequencing['strand'][len(sequencing['strand'])-1]
    ### call + or -
    gt = [sequencing['gt'][i] for i in range(0,2)] # +
    if strand == '-':
            sequencing['strand'].append('+')
            stat['gt-'] +=1
            sequencing['OK']['gt-'] += 1
            gt = [sequencing['gt'][i] for i in range(2,4)] #-
            sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'gt-' )
            if ref == 0:  qt = [sequencing['Q(gt)'][i] for i in range(4,6)] # re
            print ('PROBLEM -STRAND\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rs,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
    ### report unknow allele/error
    if ref > 0 and sum(qt) <4 : ##################### NONE?
        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'UNALL' )
        stat['unall'] +=1
        sequencing['OK']['unall'] += 1
        print('UNKNWON_ALL\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rs,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
    elif  ref == 0 and sum(qt) <2 :  ##################### NONE?
        stat['unref'] +=1
        sequencing['OK']['unref'] += 1
        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'UNREF' )
        print('UNKNWON_REF\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rs,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
    ### Call genotype
    if  sequencing['Q(gt)'][6] == 1: # same call
        sequencing['gt'].append(gt[1])
        stat['2gt'] +=1
        sequencing['OK']['2gt'] += 1
        sequencing['Q(gt)'].append('%s/%s'% (sequencing['Q(gt)'][0]+sequencing['Q(gt)'][1], tot) )
        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'2gt' )
    else: # different
        sequencing['gt'].append(gt[1]) ########## for dbSNP calls only
        if qt[0] >= qt[1] and sequencing['Q(gt)'][0] >= sequencing['Q(gt)'][1]: ## in dbsnp  should take the other one instead
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'G1>G2' )
                sequencing['Q(gt)'].append('%s/%s'% (sequencing['Q(gt)'][1], tot) )
                stat['g1'] +=1
                sequencing['OK']['g1'] += 1
                print('G1>G2\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rs,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
        elif  qt[0] >= qt[1] and sequencing['Q(gt)'][0] < sequencing['Q(gt)'][1]:  # checked
                sequencing['Q(gt)'].append('%s/%s'% (sequencing['Q(gt)'][1], tot) )
                stat['g1s2'] +=1
                sequencing['OK']['g1s2'] += 1
                print('FIT>\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rs,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'FIT>%s' % sequencing['Q(gt)'][1])
        elif  qt[0] < qt[1] and sequencing['Q(gt)'][0] < sequencing['Q(gt)'][1]: # checked
                sequencing['Q(gt)'].append('%s/%s'% (sequencing['Q(gt)'][1], tot) )
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'G1<G2' )
                stat['g2'] +=1
                sequencing['OK']['g2'] += 1
        else: #checked 
                sequencing['Q(gt)'].append('%s/%s'% (sequencing['Q(gt)'][1], tot) )
                stat['g2s1'] +=1
                sequencing['OK']['g2s1'] += 1
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'FIT<%s' % sequencing['gt'][0])
                print('FIT<\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rs,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
    return sequencing
############## end make call


############ def check_reference()
# 0 1 OK
# 15 16 g g+/g-
# 17 18 gq-+ (multpile alleles)
                    
# 2 3 AG_STR - strandg+-
# 4 5 strand r+- 
# 6 7 rev +-
# 10 11 AR_A +- ->
# 12 AR_RG + -> ref
                    
# 9 G_NOAR [unresolved]

# 13 AR_RG - [NEVER]
# 8 ARG_NO [NEVER]
# 14 ARGNO [NEVER]                    
#def check_reference(rsnum, sequencing,stat,fprob):
def check_reference(rsnum, sequencing,stat):
    allgt =  sum([sequencing['Q(gt)'][i] for i in range(2,6)]) >0
    ref = sum(sequencing['R(gt)']) >0
    refpm =  sum([sequencing['R(gt)'][i] for i in range(0,2)]) >= sum([sequencing['R(gt)'][i] for i in range(2,4)])
    plus_minus = sum([sequencing['Q(gt)'][i] for i in range(2,4)]) >= sum([sequencing['Q(gt)'][i] for i in range(4,6)])
    strand =  sequencing['strand'][0]
    if len( sequencing['strand']) > 1:
            strand =  sequencing['strand'][len(sequencing['strand'])-1]
    if sum(sequencing['Q(ref)']) >=1: 
            ### ref all same strand
        if allgt:
            if strand =='+' and sequencing['Q(ref)'][0]==1 or strand =='-' and sequencing['Q(ref)'][1]==1:
                if  strand =='+' and plus_minus:
                    sequencing['strand'].append(strand) # strand OK
                    sequencing['OK']['OK+'] += 1 #ok+
                    stat['OK+'] +=1
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'OK+')
                    #fprob.write('OK+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                elif  strand =='-' and plus_minus == False:
                    sequencing['strand'].append(strand) # strand OK
                    sequencing['OK']['OK-'] += 1 #ok+
                    stat['OK-'] +=1
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'OK-')
                    print ('OK-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                elif plus_minus:
                    sequencing['strand'].append('+')
                    if ref: # -ref=call fit allele +/-  
                        stat['g+'] +=1
                        sequencing['OK']['g+'] += 1
                        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'g+')
                        #fprob.write('g+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                    elif  sum([sequencing['Q(gt)'][i] for i in range(2,4)]) ==  sum([sequencing['Q(gt)'][i] for i in range(4,6)]):
                        stat['gq+'] +=1 # multiple allele: ref=all, call=all, ref!=call -> trust the sign +
                        sequencing['OK']['gq+'] += 1
                        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'gq+')
                        print ('gq+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))     
                    else:
                        stat['strandg+'] +=1
                        sequencing['OK']['strandg+'] += 1
                        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AG_STR-+')
                        print('AG_STR-+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                else:
                    if ref:  # +ref=call fit allele +/-  
                        stat['g-'] +=1
                        sequencing['OK']['g-'] += 1
                        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'g-')
                        print('g-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                    elif  sum([sequencing['Q(gt)'][i] for i in range(2,4)]) ==  sum([sequencing['Q(gt)'][i] for i in range(4,6)]):
                        stat['gq-'] +=1 # multiple allele: ref=all, call=all, ref!=call -> trust the sign -
                        sequencing['OK']['gq-'] += 1
                        print sequencing
                        print('SHOULNOTBEHERE_MULTI-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                    else:
                        stat['strandg-'] +=1
                        sequencing['OK']['strandg-'] += 1
                        sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AG_STR+-')
                        print ('AG_STR+-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
            ### ref all change strand
            elif  strand =='+' and sequencing['Q(ref)'][1]==1 or strand =='-' and sequencing['Q(ref)'][0]==1:
                sequencing['strand'].append('+') # trust + strand
                if  strand=='+' and plus_minus == False:
                    stat['strandr-'] +=1
                    sequencing['OK']['strandr-'] += 1
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_STR+-') # Check should stay
                    print('AR_STR+-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))  
                elif strand=='-' and plus_minus:
                    stat['strandr+'] +=1
                    sequencing['OK']['strandr+'] += 1
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_STR-+') # Check should stay
                    print('AR_STR-+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))  
                elif plus_minus: ### unresolved
                    stat['rev+'] +=1
                    sequencing['OK']['rev+'] += 1
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_REV+')
                    print('AR_REV+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                else:  ### unresolved conflict 
                    stat['rev-'] +=1
                    sequencing['OK']['rev-'] += 1
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_REV-')
                    print('AR_REV-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
            else:
                stat['ARG_NO'] +=1
                sequencing['OK']['ARG_NO'] += 1
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'ARG')
                if  sequencing['genotyping'] != 'DATA':
                    print sequencing
                    print('SHOULNOTBEHERE STR_ALL_NO_GT\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                else:
                    print('ARG_NO-GENO\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))   
        else:# ref all dif strand different?
            stat['G_NOAR'] +=1
            sequencing['OK']['G_NOAR'] += 1
            sequencing['strand'].append('+') # keep/change? + strand change
            sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'G_NOAR')
            #fprob.write
            print('GT_NO_ALLR\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1])), sequencing
    else: # ref dont fit all
        if allgt: # gt fit all
            if plus_minus:
                stat['AR_A+'] +=1
                sequencing['OK']['AR_A+'] += 1
                sequencing['strand'].append('+') # keep/change? + strand change
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_A+')
                print ('ALLREF_ALL+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
            else:
                stat['AR_A-'] +=1
                sequencing['OK']['AR_A-'] += 1
                sequencing['strand'].append('+') # keep/change? + strand change
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_A-')
                print sequencing
                print('SHOULNOTBEHERE ALLREF_ALL-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
        else:
            if ref: # gt fit ref
                if refpm: # if comparison + > -
                    stat['AR_RG+'] +=1
                    sequencing['OK']['AR_RG+'] += 1
                    sequencing['strand'].append('+') # keep/change? + strand change
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_RG+')
                    print('ALLREF_NOALL_REFGT+\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
                else :
                    stat['AR_RG-'] +=1
                    sequencing['OK']['AR_RG-'] += 1
                    sequencing['strand'].append('+') # keep/change? + strand change
                    sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'AR_RG-')
                    print('ALLREF_NOALL_REFGT-\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
            else: # gt no ref no all
                stat['ARGNO'] +=1
                sequencing['OK']['ARGNO'] += 1
                sequencing['strand'].append('+') # keep/change? + strand change
                sequencing['sequencing'] = fill_problem(sequencing['sequencing'],'ARGNO')
                print sequencing
                print('SHOULNOTBEHERE ALLREF_NOALL_NOREF\tSEQ\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(rsnum,sequencing['chr'],sequencing['ref'][0],sequencing['strand'],sequencing['alleles'],sequencing['gt'][0],sequencing['gt'][1]))
    return sequencing,stat
############### end check_reference


###############
#def genotype_call(rs,seqdata,stattable, fprob):
def genotype_call(rs,seqdata,stattable,stat):
        refgt =  seqdata['alleles'][0]
        if len( seqdata['alleles']) > 1:
                refgt =  seqdata['alleles'][ len( seqdata['alleles'])-1]
        for i in range(0,4):
                if i<2:
                    seqdata['gt'].append(reverse_strand( seqdata['gt'][i], rs, seqdata['chr']))
                    if i<1:
                        seqdata['ref'].append(reverse_strand( seqdata['ref'][i], rs, seqdata['chr']))
                    seqdata['Q(ref)'].append(compare_alleles (seqdata['ref'][i],refgt))
                seqdata['Q(gt)'].append(compare_alleles (seqdata['gt'][i],refgt))
                seqdata['R(gt)'].append(compare_alleles (seqdata['gt'][i], seqdata['ref'][0]))
        seqdata['Q(gt)'].append(compare_genotype(seqdata['gt'][0], seqdata['gt'][1]))
        seqdata, stattable = check_reference(rs,seqdata,stattable)
        seqdata = make_call(rs,seqdata,stattable,stat)
        return seqdata, stattable
############ end genotype_call



################
def get_sequencing_data(data,rsid,chrnum,stat,nfindel):
    print data
    sequencing = {}
    sequencing['dbsnp'] = rsid
    sequencing['dbsnp2'] = rsid
    sequencing['chr'] = data[0].split('chr')[1]
    sequencing['pos'] = []
    sequencing['pos'].append(string.atoi(data[1]))
    sequencing['OK'] = {} # record all stat added for that SNP
    for h in stat: sequencing['OK'][h] = 0
    sequencing['ref'] = []
    sequencing['ref'].append(data[4])
    sequencing['Q(snp)'] = []
    sequencing['Q(snp)'].append(string.atoi(data[5]))
    sequencing['Q(ref)'] = [] # comp to alleles
    sequencing['base'] =[] # base_used
    for i in range(0,4): sequencing['base'].append(data[i+10])
    sequencing['alleles'] =[]
    #sequencing['alleles'].append(data[11+4])
    sequencing['strand'] = []
    sequencing['strand'].append('+')
    #sequencing['strand'].append(data[12+4])
    sequencing['R(gt)']=[] # comp to reference gt1/+ gt1/- gt2/+ gt2/-  same in revgt1...
    sequencing['gt']=[] # all genotype calls
    sequencing['gt'].append(data[6]) 
    sequencing['gt'].append(data[8])
    sequencing['Q(gt)']=[] # comp to Alleles and each other
    sequencing['Q(gt)'].append(string.atoi(data[7]))
    sequencing['Q(gt)'].append(string.atoi(data[9]))
    sequencing['sequencing']=nfindel
    sequencing['genotyping']='NA'
    sequencing['SNPInfo']='NA'
    sequencing['gti'] = ''
    sequencing['gtg'] = []        
    if  sequencing['chr'] !=chrnum:
        print data
        print('CHROMOSOME\tSEQ\t%s\t%s\t%s' % (rsid, chrnum,sequencing['chr'])) #check
    return sequencing
########### end .get_sequencing_data


#############
def compare_genotype(genotype1,genotype2):
    print genotype1,genotype2
    gt1 = remove_slash(genotype1)
    gt2 = remove_slash(genotype2)
    print gt1,gt2
    gt='%s%s'% (gt1[1],gt1[0])
    print gt
    if gt1 == gt2 or gt==gt2: return 1
    else: return 0
############# end compare_genotype
        
#########
def fill_problem(pbl,str0):
    if(pbl == 'NA'): return '%s' % str0
    else: return '%s_%s' % (pbl,str0)
######## end fill_problem

###############
def add_stat(header,stat,tot):
    for h in header:
        if h != 'File'  and h != 'time(s)':
            tot[h] +=  stat[h]
    return tot
############## end write_stat
###############
def write_stat(header,stat,fstat):
        for h in header:
                fstat.write("%s\t" % stat[h])
        fstat.write("\n")
        fstat.flush()
############## end write_stat
############### 
def get_genotyping_data(filename):
        genotyping={} # data from genotyping 
        rsgeno =[] # set of SNP in genotype data
        fgt = open(filename) 
        numid = 0 # number of SNP in GT data
        for line in fgt.readlines():
                word=(line.strip().split()) # word in genotyping data
                if numid >0 :           
                        idg = word[0].lower()  # rs number in GT data
                        rsgeno.append(idg) 
                        genotyping[idg] =  word   
                        numid += 1
                elif word[0] == 'SNP':  numid += 1
        fgt.close()
        return (numid,rsgeno,genotyping)
################ end  get_genotyping_data

#############
#def reverse_strand(genotype,rsnum,nchr,fileprob):
def reverse_strand(genotype,rsnum,nchr):
        genotype= remove_slash(genotype)
        gt=''
        for allele in  genotype:
                if allele == 'A':
                        allele = 'T'
                elif allele == 'T' :
                        allele = 'A'
                elif allele == 'C':
                        allele = 'G'
                elif allele == 'G' :
                        allele = 'C'
                else:
                        #fileprob.write
                        if allele != '-':
                            print('UNKNOW_ALLELES\tSEQ\t%s\t%s\t%s\t%s' % (rsnum,nchr,genotype, allele)) #check
                gt = '%s%s' % (gt,allele)
        return gt
############### wnd rev_strand
#############
def compare_alleles (genotype1,genotype2):
        gt1 = remove_slash(genotype1)
        gt2 = remove_slash(genotype2)
        sameall = 0
        for a1 in  gt1:
                for a2 in gt2:
                        if a1 == a2:
                                sameall += 1
        return sameall
############ end compare_allele
############ remove_slash
def remove_slash(genotype):
        gt = re.split('[/]',genotype)
        if len(gt ) >1 :
                g=''
                for a in gt:
                        g= '%s%s' %(g,a)
                return g
        else: return gt[0]
############ end remove slash

###############
def write_output(header,data,fout,t):
    for h in header:
        if h != '#':
            #if h == 'pos':
            #      fout.write("%s\t" % data[h])
            if len(data[h]) == 1  :
                fout.write("%s\t" % data[h])
            elif h == 'ref' :
                fout.write("%s\t" % data[h][0])
            elif h == 'alleles' :
                fout.write("%s\t" % data[h])
            elif h == 'Q(gt)' or (h == 'gt' and t== 1): # output
                fout.write("%s\t" % data[h][len(data[h])-1])
            #elif h == 'gt' and t== 1: # ing
            #       fout.write("%s\t" % data[h][len(data[h])-1])
            elif h == 'strand' :
                fout.write("%s\t" % (data[h]))
            else:
                fout.write("%s\t" % data[h])
    fout.write("\n")
    fout.flush()
############## end write_stat
