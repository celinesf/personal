#!/usr/bin/python

# get data from phased data file
# output in L2'-format of haplotype for each individuals
# use the 'place' in IDHap s order for each hapname
# need to output hapname, chr, position1, pos2, genotype
 
import os
import sys
import gzip
import MySQLdb #@UnresolvedImport

def get_argv(argv):
    directory='PhasedData'
    start='phased.unphased_chr'
    end='phased.gz'
    dir_out = os.getcwd( )
    fout='haplotype'
    if len(argv) >1:
        i=0
        for arg in argv:
            if arg == '-i': directory = sys.argv[i+1]
            elif arg == '-s':  start = sys.argv[i+1]  
            elif arg == '-e': end = sys.argv[i+1]  
            elif arg == '-o':  fout = sys.argv[i+1]
            elif arg == '-d': dir_out = sys.argv[i+1]
            i+=1
    return directory, start, end,dir_out, fout
 
def get_phased_data(dirdata, start,end ):
    hap_geno={} # name of member, rsid[allele 1, allele 2, chr, position]
    os.system('mkdir %s' % dirdata)
    phased_name = [i for i in os.listdir(dirdata) if i.startswith(start) and i.endswith(end) ]
    nf = 0
    for fname in phased_name:
        phased_file = gzip.open('%s/%s' %(dirdata,fname), 'r:gz')
        nl = 0
        for line in phased_file:
            word=line.strip().split()
            if nf == 0 and nl ==0:
                ind=0
                for w in range(2,len(word),2):
                    hap_geno[ind]={}
                    hap_geno[ind]['ID'] = word[w]
                    ind +=1
                    w += 2
            elif nl != 0:
                for ind in hap_geno:
                    hap_geno[ind][word[1]]={}
                ind = 0
                for w in range(2,len(word),2):
                    hap_geno[ind][word[1]]['A1']= word[w]
                    w += 1
                    hap_geno[ind][word[1]]['A2']= word[w]
                    ind +=1
            nl += 1
        nf += 1
    return hap_geno
    
def get_IDHap_info():
    db = MySQLdb.connect(host="192.168.2.104",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    hap_list={}
    sql= "select distinct dbHAP from IDHap "
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab: 
            hap_list[s[0]]=[]
    except: 
        print "Error: unable to fecth data %s" % sql  
    sql= "select dbHAP,dbsnp, place from IDHap "
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab:
            snp_dic={}
            for (name, value) in zip(desc, s) :
                snp_dic[name[0]] = value 
            
            hap_list[snp_dic['dbHAP']].append(snp_dic)
    except: 
        print "Error: unable to fecth data %s" % sql 
    cursor.close ()
    sorted_hap_list ={}
    for hap in hap_list:
        from operator import itemgetter
        sorted_hap_list[hap]= sorted(hap_list[hap], key=itemgetter('place')) 
    return  sorted_hap_list
    
def get_genotype(ind_geno,haplotype):
    genotype = {}
    genotype['A1'] = genotype['A2'] = ""
    for place in range(0,len(haplotype)):
        if place+1 == haplotype[place]['place']:
            try:
                genotype['A1'] = "%s%s" % (genotype['A1'],ind_geno[haplotype[place]['dbsnp']]['A1'])
                genotype['A2'] = "%s%s" % (genotype['A2'],ind_geno[haplotype[place]['dbsnp']]['A2'])           
            except:  
                genotype['A1'] = genotype['A2'] = "NA"
                break
        else:
            print "ISSUE Place at %s" % haplotype
    return genotype

if __name__ == '__main__':
    
    dirdata, start,end, dir_out, output_name =get_argv(sys.argv) #recover arguments    
    snp_geno = get_phased_data(dirdata, start,end)
    IDHap = get_IDHap_info()
    
    for ind in snp_geno:
        output_file = open('%s/%s_%s' % (dir_out,output_name,snp_geno[ind]['ID']),'w')
        output_file.write('Haplotype_name\tgenotype\n')
        for hap in IDHap:
            genotype = get_genotype(snp_geno[ind],IDHap[hap])
            output_file.write('%s\t%s,%s\n' %(hap,genotype['A1'],genotype['A2']))
        output_file.close()
        
    cmdline = 'rm -r %s' %(dirdata)
    os.system(cmdline)
    
    pass