#!/usr/bin/python

# select distinct dbSNP from IDHap 
#  Read L2_prim data and extract the genotype for each dbsnp
# write in chrname files each dbsnp genotype order by positions
# in input file format for beagle

# 360 : get_argv, get_genotype_file_names
# output the marker files

import os
import sys
import MySQLdb #@UnresolvedImport
import re
import linecache
import copy

def get_genotype_file_names(directory,start,end):
    return [i for i in os.listdir(directory) if i.startswith(start) and  i.endswith(end)]

def get_argv(argv):
    directory=os.getcwd( )
    start='QCed_'
    end=''
    fout='unphased'
    dir_out='UnphasedData'
    if len(argv) >1:
        i=0
        for arg in argv:
            if arg == '-i': directory = sys.argv[i+1]
            elif arg == '-s':  start = sys.argv[i+1]  
            elif arg == '-e': end = sys.argv[i+1]  
            elif arg == '-o':  fout = sys.argv[i+1]
            elif arg == '-d': dir_out = sys.argv[i+1]
            i+=1
    return directory, start, end, dir_out,fout

def get_snp_list():
    db = MySQLdb.connect(host="192.168.2.104",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    sql= "select distinct chr,chr_pos,rs, Q.MajorAllele, Q.MinorAllele from DBSNP_temp S, SNPInfo Q where S.rs=Q.dbsnp and rs in (select distinct dbsnp from IDHap) and genome='GRCh37.p5' "
    snp_list=[]
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab:
            snp_dic={}
            for (name, value) in zip(desc, s) :
                snp_dic[name[0]] = value 
            snp_list.append(snp_dic)
    except: 
        print "Error: unable to fecth data %s" % sql    
    cursor.close ()
    from operator import itemgetter
    sorted_snp_list = sorted(snp_list, key=itemgetter('chr_pos')) 
    return sorted_snp_list

def get_genotype(snp_list,directory,file_list):
    genotype={}
    for f_name in file_list: 
        genotype[f_name]=[]
        file_name='%s/%s'%(directory,f_name)
        data_file = open(file_name,'r')
        data = data_file.read()
        for snp in snp_list:
            tpsnp=copy.deepcopy(snp)
            tpsnp['chr_pos']=snp['chr_pos']
            tpsnp['genotype']= '??'
            i=0
            for m in re.finditer('\t%s\t' % snp['chr_pos'], data.lower())  :
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

def hapmap_file_ID(hapmap_id,data,file_type, nfile,nindv,dir_unphased, chromosome):
    hapmap_id.append({})
    hapmap_id[nfile]['nfile'] ='%s/%s%d_chr%s' % (dir_unphased,file_type,nfile, chromosome)
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
        hapmap_id = hapmap_file_ID(hapmap_id,data,file_type, nfile,i+1,dir_unphased,  chromosome)
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

def Beagle_input(chromosome,snp_list, genotype,trio_data, unrelated_data, dir_unphased,output_name):
    # 128 alleles only per file
# get ceu urelated data
# 1 male, 2 female, indv (2 father 3 mother)
    marker_file = open('%s/markers_chr%s' % (dir_unphased, chromosome),'w')
    
    beagle_file = open('%s/%s_chr%s' % (dir_unphased,output_name, chromosome),'w')
    beagle_file.write("I\tid\t")   
    for g in genotype:
        beagle_file.write("%s\t%s\t" % (g,g)) 
    beagle_file.write("\n" ) 

    trio_file=[] 
    trio_file = hapmap_file_ID(trio_file,trio_data,'trios',0,0,dir_unphased,  chromosome)

    unrelated_file=[] 
    unrelated_file = hapmap_file_ID(unrelated_file,unrelated_data,'unrelated',0,0,dir_unphased,chromosome)
    
    pos1 = pos2 = 0
    snp_no=0
    for snp in snp_list:
        marker_file.write('%s\t%s\t%s\t%s\t%s\n' % (snp['rs'],snp['chr'],snp['chr_pos'],snp['MajorAllele'],snp['MinorAllele']))
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
        pos1=pos2
    beagle_file.close()
    marker_file.close()

def get_hapmap_indv_data():
    db = MySQLdb.connect(host="192.168.2.104",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    hapmap_data = {}
    sql = "select * from HapMapFamily"
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab:
            hapmap_dic={}
            for (name, value) in zip(desc, s) :
                hapmap_dic[name[0]] = value 
            hapmap_data[hapmap_dic['indv1']] = hapmap_dic
    except: 
        print "Error: unable to fecth data %s" % sql 
        
    sql = "select * from FamilyAllele"
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab:
            hapmap_dic={}
            for (name, value) in zip(desc, s) :
                hapmap_dic[name[0]] = value 
            try:
                hapmap_data[hapmap_dic['Indv']]['genotype'][hapmap_dic['dbSNP']] = hapmap_dic['A1']+ hapmap_dic['A2']
            except:
                hapmap_data[hapmap_dic['Indv']]['genotype'] = {}
    except: 
        print "Error: unable to fecth data %s" % sql  
    cursor.close ()
    return hapmap_data 

def have_data(data,data_type):
    try:
        tp=data[data_type]
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

if __name__ == '__main__':
    dirdata, start, end, dir_unphased,output_name =get_argv(sys.argv) #recover arguments
    trio_data, unrelated_data = get_hapmap_data()
    listmembers = get_genotype_file_names(dirdata,start,end ) # list o datafile with genotype to extract
    snp_list= get_snp_list() #order by position 
    os.system('mkdir %s' % dir_unphased)
    for chromosome in range(1,22):
        chr_snp_list = [i for i in snp_list if i['chr'] == '%d' % chromosome]   # list of snp order for chr 
        if len(chr_snp_list)>1:
            genotype = get_genotype(chr_snp_list,dirdata,listmembers)
            Beagle_input(chromosome,chr_snp_list,genotype,trio_data, unrelated_data ,dir_unphased,output_name)     
    pass