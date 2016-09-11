#!/usr/bin/python

# select distinct dbSNP from IDHap 
#  Read L2_prim data and extract the genotype for each dbsnp
# write in chrname files each dbsnp genotype order by positions
# in input file format for beagle
# 360 : get_argv, get_genotype_file_names
# output the marker files
# run beagle for all chromosome
# select distinct dbSNP from IDHap 
#  Read L2_prim data and extract the genotype for each dbsnp
# write in chrname files each dbsnp genotype order by positions
# in input file format for beagle
# 360 : get_argv, get_genotype_file_names
# output the marker files


import os
import sys
import gzip
import shlex, subprocess
import pymongo
import re
import linecache
import copy
import l1l2


def get_argv(argv):
    file_illu = 'output.txt'
    dir_illu = os.getcwd()
    if len(argv) >1:
        dir_illu =sys.argv[1]
        file_illu =sys.argv[2]
        
    global DISEASE,IDHAP,HAPMAP_ID,HAPMAP_DATA
    cw = pymongo.Connection(l1l2.MONGO, l1l2.PORT);
    DISEASE = cw[l1l2.GENOPHENDB][l1l2.DISEASE_COLLECTION];
    IDHAP = cw[l1l2.GENETICSBD][l1l2.IDHAP_COLLECTION];
    HAPMAP_ID = cw[l1l2.GENETICSBD][l1l2.HAPMAP_ID_COLLECTION];
    HAPMAP_DATA = cw[l1l2.GENETICSBD][l1l2.HAPMAP_DATA_COLLECTION];

    return  dir_illu, file_illu

def add_hapmap_files(dir_illu,chr_num,cmdline, dirdata,type_file,type_beagle):
    file_name = [i for i in os.listdir('%s/%s' % (dir_illu,dirdata)) if i.startswith(type_file) and i.endswith('_chr%s' % chr_num.split('_chr')[1])] 

    for t in file_name:
        cmdline= '%s %s=%s/%s/%s ' % (cmdline, type_beagle,dir_illu,dirdata,t)
    return cmdline

def get_snp_list():
    #genophen30 = CW["genophen30"];
    #finddoccursor = genophen30.genetics.rep.find({});     
    finddoccursor = DISEASE.find({});
    snp_list = []    
    hap_list = []
    for c in finddoccursor:     
        for snp in c['snps']:
            if len(snp['snp'].split('rs'))<=1:
                if snp['snp'] not in hap_list:
                    hap_list.append(snp['snp'])
                    snp_list = get_hap_snp(snp['snp'], snp_list)

    from operator import itemgetter
    sorted_snp_list = sorted(snp_list, key=itemgetter('chr_pos')) 
    return sorted_snp_list
'''
update get_hapID
'''
def get_hapID(hap) :
    snplist=[]
#    geneticsdb = CW["genetics"];
#    finddoccursor = geneticsdb.idhap.find({'dbhap':hap});     
    finddoccursor = IDHAP.find({'dbhap':hap});              
    for c in finddoccursor:
        snplist.append(c)
    return snplist 
     
'''
update get_hap_snp
'''
def get_hap_snp(hap,snplist) :
#    geneticsdb = CW["genetics"];
#    finddoccursor = geneticsdb.idhap.find({'dbhap':hap});
    finddoccursor = IDHAP.find({'dbhap':hap}); 
    for c in finddoccursor:
        newsnp = {}
        newsnp['rs'] = c['dbsnp']
        newsnp['chr'] = c['chromosome']
        newsnp['chr_pos'] = c['position']
        newsnp['major_allele'] = c['major_allele']
        newsnp['minor_allele'] = c['minor_allele']
        if newsnp not in snplist:
            snplist.append(newsnp)
    return snplist
        
def get_genotype_l2(snp_list,directory,f_name):
    genotype={}   
    genotype[f_name]=[]
    file_name='%s/%s'%(directory,f_name)
    data_file = open(file_name,'r')
    data = data_file.read()
    for snp in snp_list:
        tpsnp=copy.deepcopy(snp)
        #print tpsnp
        tpsnp['chr_pos'] = snp['chr_pos']
        tpsnp['genotype']= '??'
        i=0
        ############ TO CHANGE ############
        #for m in re.finditer('\t%s\t' % snp['chr_pos'], data.lower())  :
        for m in re.finditer('%s\t' % snp['rs'], data.lower())  :
            if i==0:
                start = m.start()
                lineno = data.count('\n', 0, start) + 1
                line = linecache.getline(file_name, lineno)
                word =(line.strip().split()) 
                tpsnp['genotype']=word[3]
            else:
                print 'Issue I got more than one snp named %s' % snp['rs']
            i+=1 
        if i == 0:
            print 'I could not find %s at chr %s position %s in file %s.' % (snp['rs'], snp['chr'],snp['chr_pos'],f_name)
        genotype[f_name].append(tpsnp)
    data_file.close() 
    return genotype

def hapmap_file_ID(hapmap_id,data,file_type, nfile,nindv,dir_unphased, chromosome, dir_illu):
    hapmap_id.append({})
    hapmap_id[nfile]['nfile'] ='%s/%s/%s%d_chr%s' % (dir_illu,dir_unphased,file_type,nfile, chromosome)
    hapmap_id[nfile]['indv']=[]
               
    hamap_file = open(hapmap_id[nfile]['nfile'],'w')
    hamap_file.write("I\tid\t")
    for i in range(nindv,nindv+60):
        try:
            hapmap_id[nfile]['indv'].append(copy.deepcopy(data[i]))
            hamap_file.write("%s\t%s\t" % (data[i]['indv1'],data[i]['indv1']))
        except:
            break;

    if i+1 == nindv+60:
        hamap_file.write('\n')
        hamap_file.close()
        nfile += 1
        hapmap_id = hapmap_file_ID(hapmap_id,data,file_type, nfile,i+1,dir_unphased,  chromosome, dir_illu)
    else:
        hamap_file.write('\n')
        hamap_file.close()
    return hapmap_id

def hapmap_beagle(snp, data):
    for nfile in range(0,len(data)):
        beagle_file = open(data[nfile]['nfile'],'a')
        beagle_file.write("M\t%s\t" % snp['rs'])
        for g in data[nfile]['indv']:
            try:
                if  g['genotype'][snp['rs']] != 'NN': # write both known alleles
                    beagle_file.write("%s\t%s\t" % (g['genotype'][snp['rs']][0] , g['genotype'][snp['rs']][1]))
                else:
                    beagle_file.write("?\t?\t")
            except:
                beagle_file.write("?\t?\t") 
        beagle_file.write("\n" ) 
        beagle_file.close()       

def Beagle_input(chromosome,snp_list, genotype,trio_data, unrelated_data, dir_unphased,output_name, dir_illu):
    # 128 alleles only per file
# get ceu urelated data
# 1 male, 2 female, indv (2 father 3 mother)
    marker_file = open('%s/%s/markers_chr%s' % (dir_illu,dir_unphased, chromosome),'w')
    
    beagle_file = open('%s/%s/%s_chr%s' % (dir_illu,dir_unphased,output_name, chromosome),'w')
    beagle_file.write("I\tid\t")   
    for g in genotype:
        beagle_file.write("%s\t%s\t" % (g,g)) 
    beagle_file.write("\n" ) 

    trio_file=[] 
    trio_file = hapmap_file_ID(trio_file,trio_data,'trios',0,0,dir_unphased,  chromosome, dir_illu)

    unrelated_file=[] 
    unrelated_file = hapmap_file_ID(unrelated_file,unrelated_data,'unrelated',0,0,dir_unphased,chromosome, dir_illu)
    
    pos1 = pos2 = 0
    snp_no=0
    for snp in snp_list:
        marker_file.write('%s\t%s\t%s\t%s\t%s\n' % (snp['rs'],snp['chr'],snp['chr_pos'],snp['major_allele'],snp['minor_allele']))
        pos2 = snp['chr_pos']
        if pos1 >= pos2:
            print "PROBLEM POSITION NOT IN ORDER", pos1 , pos2,snp
        else:
            hapmap_beagle(snp,trio_file)    
            hapmap_beagle(snp,unrelated_file) 
            beagle_file.write("M\t%s\t" % snp['rs'])
            for g in genotype:
                if genotype[g][snp_no]['genotype'] != 'NA' and len(genotype[g][snp_no]['genotype']) == 2: # write both known alleles
                    beagle_file.write("%s\t%s\t" % (genotype[g][snp_no]['genotype'][0] , genotype[g][snp_no]['genotype'][1]))
                elif genotype[g][snp_no]['genotype'] != 'NA' and len(genotype[g][snp_no]['genotype']) == 1: # one allele and one missing data
                    beagle_file.write("%s\t?\t" % (genotype[g][snp_no]['genotype'][0] ))
                else:# missing data
                    beagle_file.write("?\t?\t")
            beagle_file.write("\n" ) 
        snp_no += 1
        pos1 = pos2
    
    beagle_file.close()
    marker_file.close()

def get_hapmap_indv_data():
#    geneticsdb = CW["genetics"];
#    finddoccursor = geneticsdb.HapMapFamily.find({})
    finddoccursor =  HAPMAP_ID.find({})
    hapmap_data = {}
    for c in finddoccursor:
        hapmap_data[c['indv1']] = c
        
#    finddoccursor = geneticsdb.FamilyAllele.find({})
    finddoccursor =  HAPMAP_DATA.find({})
    for c in finddoccursor:
        try:
            hapmap_data[c['Indv']]['genotype'][c['dbSNP']] = c['A1']+ c['A2']
        except:
            hapmap_data[c['Indv']]['genotype'] = {}
    return hapmap_data 

def have_data(data,data_type):
    try:
        data[data_type]
        return True
    except:
        return False

def get_hapmap_trio(hapmap_data):
    trio_data = []
    for i in hapmap_data:
        if have_data(hapmap_data[i],'genotype')is True :
            mother = hapmap_data[i]['indv3'] 
            father = hapmap_data[i]['indv2'] 
            if father != '0' and  mother != '0':
                if hapmap_data[father]['Gender'] == 1 \
                and hapmap_data[mother]['Gender'] == 2 \
                and hapmap_data[mother]['Family'] == hapmap_data[i]['Family'] \
                and hapmap_data[father]['Family'] == hapmap_data[i]['Family'] :            
                    if have_data(hapmap_data[mother],'genotype')is True  and have_data(hapmap_data[father],'genotype')is True :
                        hapmap_data[mother]['in']=hapmap_data[father]['in']=hapmap_data[i]['in']=0
                        trio_data.append(copy.deepcopy(hapmap_data[mother]))
                        trio_data.append(copy.deepcopy(hapmap_data[father]))
                        trio_data.append(copy.deepcopy(hapmap_data[i]))  
                    elif have_data(hapmap_data[mother],'genotype')is True:
                        hapmap_data[i]['in']=0
                    elif have_data(hapmap_data[father],'genotype')is True :  
                        hapmap_data[i]['in']=0          
                else:
                    print "ISSUE 00 with Family %s of indv %s" % (hapmap_data[i]['Family'] ,i)
            elif (father != '0' ) :
                if hapmap_data[father]['Gender'] == 1 \
                and hapmap_data[father]['Family'] == hapmap_data[i]['Family'] :              
                    if have_data(hapmap_data[father],'genotype') is True :
                        hapmap_data[i]['in']=0 # only use the father unrelated 
                else:
                    print "ISSUE dad with Family %s of indv %s" % (hapmap_data[i]['Family'] ,i)
            elif (mother != '0' ): 
                if hapmap_data[mother]['Gender'] == 2 \
                and hapmap_data[mother]['Family'] == hapmap_data[i]['Family'] :  
                    if have_data(hapmap_data[mother],'genotype') is True :
                        hapmap_data[i]['in']=0 # only use the mother unrelated 
                else:
                    print "ISSUE mum with Family %s of indv %s" % (hapmap_data[i]['Family'] ,i)
   
    return trio_data,hapmap_data

def get_hapmap_unrelated(hapmap_data):
    unrelated_data = []
    for i in hapmap_data:
        if have_data(hapmap_data[i],'genotype') is True and have_data(hapmap_data[i],'in') is False :
            hapmap_data[i]['in']=0
            unrelated_data.append(copy.deepcopy(hapmap_data[i]))
    return unrelated_data
             
def get_hapmap_data():
    hapmap_data = get_hapmap_indv_data()
    trio_data,hapmap_data = get_hapmap_trio(hapmap_data)
    unrelated_data = get_hapmap_unrelated(hapmap_data)
    return  trio_data,unrelated_data

def get_phased_data(dir_illu,dirdata, start,end ):
    hap_geno={} # name of member, rsid[allele 1, allele 2, chr, position]
    phased_name = [i for i in os.listdir('%s/%s' % (dir_illu,dirdata)) if i.startswith(start) and i.endswith(end) ]
    nf = 0
    for fname in phased_name:
        phased_file = gzip.open('%s/%s/%s' %(dir_illu,dirdata,fname), 'r:gz')
        nl = 0
        for line in phased_file:
            word=line.strip().split()
            if nl != 0:
                hap_geno[word[1]]={}
                for w in range(2,len(word),2):
                    hap_geno[word[1]]['A1']= word[w]
                    w += 1
                    hap_geno[word[1]]['A2']= word[w]  
            nl += 1
        nf += 1
    return hap_geno
    
def get_IDHap_info():
    
#    genophen30 = CW["genophen30"];
#    finddoccursor = genophen30.genetics.rep.find({});     
    finddoccursor = DISEASE.find({});       
    hap_list={}         
    for c in finddoccursor:      
        for snp in c['snps']:
            if len(snp['snp'].split('rs'))<=1:
                if snp['snp'] not in hap_list:
                    hap_list[snp['snp']] = get_hapID(snp['snp'])
    
    sorted_hap_list ={}
    for hap in hap_list:
        from operator import itemgetter
        sorted_hap_list[hap]= sorted(hap_list[hap], key=itemgetter('position')) 
    return  sorted_hap_list
      
def get_genotype(ind_geno,haplotype):
    genotype = {}
    genotype['A1'] = genotype['A2'] = ""
    chromosome = 0
    position = 0
    for place in range(0,len(haplotype)):
        if chromosome != 0 and chromosome !=  haplotype[place]['chromosome']:
            print "PROBLEM CHROMOSOME NUMBER HAPLO %s SNP %s" %(haplotype[place]['dbhap'], haplotype[place]['dbsnp'])
        chromosome =  haplotype[place]['chromosome']
        if position == 0:
            position = haplotype[place]['position']
        elif haplotype[place]['position'] < position:
            position = haplotype[place]['position']
            
        try:
            genotype['A1'] = "%s%s" % (genotype['A1'],ind_geno[haplotype[place]['dbsnp']]['A1'])
            genotype['A2'] = "%s%s" % (genotype['A2'],ind_geno[haplotype[place]['dbsnp']]['A2'])           
        except:  
            genotype['A1'] = genotype['A2'] = "NA"
            break
    return chromosome, position, genotype

def prep_beable_input(dir_l2,file_l2,dir_unphased,funphased):
    print "---> Phasing file: %s" % file_l2
    trio_data, unrelated_data = get_hapmap_data()
    snp_list= get_snp_list() #order by position 
    os.system('mkdir %s/%s' % (dir_illu,dir_unphased))   
    for chromosome in range(1,23):
        chr_snp_list = [i for i in snp_list if i['chr'] == '%d' % chromosome]   # list of snp order for chr 
        if len(chr_snp_list)>1:
            genotype = get_genotype_l2(chr_snp_list,dir_l2,file_l2)
            Beagle_input(chromosome,chr_snp_list,genotype,trio_data, unrelated_data ,dir_unphased,funphased,dir_l2) 
    return genotype

def run_Beagle(dir_illu,trios,unrelated,dir_unphased,dir_phased, fphased):
    unphased_name = [i for i in os.listdir('%s/%s' % (dir_illu,dir_unphased)) if i.startswith('unphased_chr') ]  
    os.system('mkdir %s/%s' % (dir_illu,dir_phased))
    for chr_num in unphased_name:  
        cmdline = 'java -Xmx1000m -jar %sbeagle.jar missing=? markers=%s/%s/%s%s out=%s/%s/%s unphased=%s/%s/%s' % (l1l2.BEAGLE_DIR,dir_illu,dir_unphased,marker, chr_num.split('_chr')[1], dir_illu,dir_phased, fphased,dir_illu,dir_unphased,chr_num)
        cmdline= add_hapmap_files(dir_illu,chr_num,cmdline, dir_unphased,trios,'trios')
        cmdline= add_hapmap_files(dir_illu,chr_num,cmdline, dir_unphased,unrelated,'unphased')
        p= subprocess.Popen(shlex.split(cmdline))
        p.communicate()      
    cmdline = 'rm -r %s/%s/*trio* %s/%s/*unrelated* %s/%s/*.log %s/%s' %(dir_illu,dir_phased,dir_illu,dir_phased,dir_illu,dir_phased,dir_illu,dir_unphased)
    os.system(cmdline)  

def postprocessing_beagle_output(dir_l2,file_l2,haplo_file,dir_phased):
    print "\n ---> Generating Haplotype genotypes for %s" % (file_l2)
    snp_geno = get_phased_data(dir_illu,dir_phased, 'phased.unphased_chr','phased.gz')
    IDHap = get_IDHap_info()
    output_file = open('%s/%s_%s' % (dir_l2,haplo_file,file_l2),'w')
    num_line = 0
    for hap in IDHap:
        chromosome, position, genotype = get_genotype(snp_geno,IDHap[hap])
        numsnp = 0
        if num_line > 0:
            output_file.write('\n') 
        for g in genotype['A1']:
            if g == '?':
                break
            numsnp += 1
        if numsnp == len(genotype['A1']) and genotype['A1'] != 'NA':
            output_file.write('%s\t%s\t%s\t%s/%s' %(chromosome, position, hap,genotype['A1'],genotype['A2']))
        else :
            output_file.write('%s\t%s\t%s\tNA' %(chromosome, position, hap))
        num_line += 1
    output_file.close()  
    cmdline = 'rm -r %s/%s' %(dir_illu,dir_phased)
    os.system(cmdline)

def append_L2_haplotype_data(dir_l2,haplo_file,file_l2):
    output_file=open('%s/%s_%s' % (dir_l2,haplo_file,file_l2),'a')
    l2file=open('%s/%s' %(dir_l2,file_l2),'r')
    snpdata=l2file.read()
    output_file.write(snpdata) 
    output_file.close()

if __name__ == '__main__':
    
    dir_illu = os.getcwd( )
    funphased='unphased'
    dir_unphased='UnphasedData'
    fphased='phased'
    dir_phased='PhasedData'
    phased_start='phased.unphased_chr'
    phased_end='phased.gz'
    haplo_file = 'haplotype' 
    trios='trios'
    marker = 'markers_chr'
    unrelated='unrelated'

    dir_illu, file_illu =get_argv(sys.argv) #recover arguments

    prep_beable_input(dir_illu,file_illu,dir_unphased,funphased) # preprocessing

    run_Beagle(dir_illu,trios,unrelated,dir_unphased,dir_phased, fphased) #phasing

    postprocessing_beagle_output(dir_illu,file_illu,haplo_file,dir_phased) #post pocessing  

    append_L2_haplotype_data(dir_illu,haplo_file,file_illu)

    pass
