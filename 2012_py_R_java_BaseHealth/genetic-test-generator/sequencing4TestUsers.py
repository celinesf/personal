#!/usr/bin/env python
'''
Created on October, 22 2012
@author: celine

Abstract classes:
1) added method
2) sql to mongo
3) random based on frequencies to be more realistic
4) positions are relative now
5) added gene names
6) all disease risk are calulated
'''

#import MySQLdb
from numpy import random
import math, csv
import pymongo, json, copy
import re
from datetime import * 
import collections
DATE = str(datetime.now())

DIRECTORY = "/Users/celine/Box Documents/Celine/ScriptInOut/sequencing4TestUsers"
# tiger_dev
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"

#mars_qa
# mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.244.24.253/admin"

CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genophen30"]
GENETICS = CW["genetics"]

class UserGenome():
    ''' to generate randome Genotypes for an imaginary individual '''
     
    def __init__(self,filename , ethn, gender, type,gid):
        self.directory = DIRECTORY
        self.genetic     = type
        self.gender     = gender
        self.ethn       = ethn.lower()
        self.study ={}
        self.id         = filename
        self.gid         = gid
        self.risk       = {}
        self.snp_info   = {}
        cursor = GENOPHEN30['general'].find({'name':'diseases'})
        self.bible_mongo ={ 'Caffeine':200, 'Lactose':201, 'Alcohol':202, 'Carbohydrates':203, 'Calcium':204,'Vitamin D':205,
                            'Metformin':300,'Citalopram':301,'Fluoxetine':302,'Clopidogrel':303,'Desipramine':304,'Irbesartan':305 }

        for doc in cursor:
            for d in doc['list']:
                if 'v' in d:
                    self.bible_mongo[d['v']] = 0
    
    ### END of generate_user_data
    
    ### generate_user_data
    ###
    def generate_user_data(self):     
#         print(' Function: generate_user_data')   

        ### all the snps in our db
        snp_info = self.get_snp_info()
#         print snp_info
#         print self.study
        snp_info = self.get_snp_info_blue(snp_info)
        snp_info = self.get_snp_genomics(snp_info)

        gid = '%s_%s' %(self.id, self.ethn+self.gender)
        outfile = open("%s/%s.txt"%(self.directory,gid), 'w')
        outfile.write("rsid\tchr_num\tpos\tgene\tGeno\t")

        for key, value in sorted(self.bible_mongo.iteritems(), key=lambda (k,v): (v,k)):
            outfile.write('%s\t' % (key))
        outfile.write("\n")
        
        ### generate_genotype_data
        new_doc, data = self.generate_genotype_data(outfile,snp_info)
        new_doc['gid'] = self.gid
        new_doc['name'] = gid   
        new_doc['file_id'] = gid  
        new_doc["active"] = True
        new_doc["ready"] = True
        new_doc["genome_type"] = "sequencing"
        new_doc["vendor"] = "genophen"

        ''' UNCHECK for SEQUENCING '''
        new_doc = self.update_genomsequence(gid, self.gid, new_doc)
        ''' END '''

        sequencing_file = open("%s/%s.json"%(self.directory,gid), 'w')
        sequencing_file.write(json.dumps(new_doc,indent=4))
        
        ### output test file
        self.output_summary(outfile,"\t")
        outfile.close()    
        
        ''' UNCHECK for 23 & me data file '''
        self.write_23andme_file(gid,snp_info,data, 36)
        self.write_23andme_file(gid,snp_info,data, 37)
        ''' END '''     
    ### END of generate_user_data
    
    ### write_23andme_file
    ###
    def write_23andme_file(self,gid,snp_info,newdoc, build):     
#         print(' Function: write_23andme_file %s' % build)   
        ### init
        file23 = open("%s/23&me_%s.txt" % (self.directory ,build),"r")
        out23 = open("%s/build_%s_%s.txt" % (self.directory ,build,gid),"w")

        ### reinit risk
        for disease in self.risk:
            self.risk[disease] = {}
            self.risk[disease]['count'] = 0
            self.risk[disease]['OR'] = 0

        for line in file23:
            if "#" in line:
                if "rsid\t" in line:
                    out23.write("%s\t" % line.strip())
                    for key, value in sorted(self.bible_mongo.iteritems(), key=lambda (k,v): (v,k)):
                        out23.write('%s\t' % (key))
                    out23.write("\n")
                else:
                    out23.write(line)
            else:
                word = line.strip().split('\t')
                snp = word[0]
                if snp in newdoc:
                    if newdoc[snp] == "NA":
                        newdoc[snp] = "--"
                    newline = "%s\t%s\t%s\t%s\t" %(word[0],word[1],word[2],newdoc[snp])
                    out23.write(newline)
                    self.write_genotype_risk(out23,snp,newdoc[snp],snp_info[snp]['disease'],"\t")  
              
#                else :
#                    out23.write(line)
        # add haplotype data
        out23.write("\n")
        for snp in snp_info:
            if 'rs' not in snp :
                out23.write('#%s\t%s\t%s\t%s\t' % (snp,str(snp_info[snp]['chr']) ,snp_info[snp]['pos'],newdoc[snp]))
                self.write_genotype_risk(out23,snp,newdoc[snp],snp_info[snp]['disease'],"\t")    

        self.output_summary(out23,"\t")
        file23.close()
        out23.close()
    ### END of write_23andme_file
    
    
    ### get_snp_genomics
    ###
    def get_snp_genomics(self, unique_snp):
#         print(' Function: get_snp_genomics')
        finddoccursor = GENOPHEN30.genomics.rep.aggregate([
            {'$match':{'$or':[
                   {'gid':'drugs'}, {'gid':'nutrition'}
                   ]}} ,
            {'$project':{'drugs':1,'nutrition':1,'_id':0}},
            {'$group':{'_id':{'drugs':'$drugs','nutrition':'$nutrition'}}}])  
        for doc1 in finddoccursor['result']:
            for  doc2 in doc1['_id']:
                pharma = doc1['_id'][doc2]
                for chem_doc in pharma:
                    chem = chem_doc['display_name']
                    for snpdoc in chem_doc['snp']:
                        snp = snpdoc['name']
                        try:
                            unique_snp[snp]
                        except:
                            unique_snp[snp] = {}
                            unique_snp[snp]['chr'] = snpdoc['chromosome_number']
                            unique_snp[snp]['gene'] = snpdoc['gene']
                            unique_snp[snp]['pos'] = snpdoc['relative_position']
                            unique_snp[snp]['disease'] = {}
                            
                            ########## TO DO ADD DESCRIPTION?
                        try:
                            unique_snp[snp]['disease'][chem]
                        except:
                            unique_snp[snp]['disease'][chem] = {'toto':{}}  
                            unique_snp[snp]['disease'][chem]['toto']={}
                            unique_snp[snp]['disease'][chem]['toto']['gt'] ={}
                        for geno in snpdoc['genotype']:
                            unique_snp[snp]['disease'][chem]['toto']['gt'][geno['genotype']] = '%s - %s' % (geno['effect'], geno['ethnicity_studied']  )
                        #unique_snp[snp]['disease'] = sorted(unique_snp[snp]['disease'], key=unique_snp[snp]['disease'].get) 
        return unique_snp
    #END OF get_snp_genomics


    ### get_snp_info
    ###
    def get_snp_info(self):
#         print(' Function: get_snp_info')
        finddoccursor = GENOPHEN30.genetics.rep.find({'gt':{'$ne':'NA'}})
        unique_snp = {}
        for doc in finddoccursor  :
            disease = doc['disease']
            for snpdoc in doc['snps']:
                snp = snpdoc['snp']
                eth =  snpdoc['ethnicity'].lower()
                study = snpdoc['study_in'].lower()
                gt = snpdoc['gt']
                if eth == self.ethn:
                    # populate self.study
                    if disease not in self.study:
                        self.study[disease] = []
                    if study not in self.study[disease] :
                        self.study[disease].append(study)
                    if 'rs' not in snp and snp not in unique_snp:
                        unique_snp, list_snp = self.get_snp_haplotypes(snp, unique_snp, disease,study)
                    if snp not in unique_snp:
                        unique_snp[snp] = {}
                        unique_snp[snp]['chr'] = snpdoc['chromosome']
                        unique_snp[snp]['gene'] = snpdoc['gene']
                        unique_snp[snp]['pos'] = snpdoc['relative_position']
                        unique_snp[snp]['disease'] = {}
                        if 'rs' not in snp:
                            unique_snp[snp]['list'] = list_snp
                    if disease not in unique_snp[snp]['disease']:
                        unique_snp[snp]['disease'][disease] = {}  
                    if study not in unique_snp[snp]['disease'][disease]  :
                        unique_snp[snp]['disease'][disease][study] ={}
                        unique_snp[snp]['disease'][disease][study]['gt'] ={}
                    unique_snp[snp]['disease'][disease][study]['gt'][gt] = snpdoc['or']  
        return unique_snp
    #END OF get_snp_info


    ### get_snp_haplotypes
    ###
    def get_snp_haplotypes(self, hap , unique_snp, disease,study):
#         print(' Function: get_snp_haplotypes %s' % hap)
        list_snp = []
        finddoccursor = GENETICS.idhap.find({'dbhap':hap}).sort([('position',1)])
        for hap_info in finddoccursor:
            snp = hap_info['dbsnp']
            list_snp.append(snp)
            if snp not in unique_snp:
                unique_snp[snp] = {}
                unique_snp[snp]['chr'] = hap_info['chromosome']
                unique_snp[snp]['gene'] = hap_info['gene_next_bio']
                unique_snp[snp]['pos'] = hap_info['relative_position']
                unique_snp[snp]['disease'] = {}
            if disease not in unique_snp[snp]['disease']:
                unique_snp[snp]['disease'][disease] = {}  
            if study not in unique_snp[snp]['disease'][disease]  :
                unique_snp[snp]['disease'][disease][study] ={}
                unique_snp[snp]['disease'][disease][study]['gt'] ={}
               
            minor = hap_info['minor_allele'] 
            major = hap_info['major_allele'] 
            geno = ['%s%s'%(major,minor),'%s%s'%(minor,minor),'%s%s'%(major,major)]
            unique_snp[snp]['haplotype'] = True
            for gt in geno:
                gtrev = self.reverse(snp,gt,hap_info['chromosome'] )
                gt_list =unique_snp[snp]['disease'][disease][study]['gt']
                if gt not in gt_list and gtrev not in gt_list:
                    unique_snp[snp]['disease'][disease][study]['gt'][gt] = "haplotype" 
                    if gt[0] != gt[1]:# add NA
                        unique_snp[snp]['disease'][disease][study]['gt']["NA"] = "haplotype" 
        return unique_snp, list_snp
    #END OF get_snp_haplotypes


    ### get_snp_info_blue
    ###
    def get_snp_info_blue(self, unique_snp):
#         print(' Function: get_snp_info_blue')
        list_snp = []
        finddoccursor = GENOPHEN30.genetics.rep.find({'ethnicity':{'$ne':self.ethn}})
        for doc in finddoccursor  :
            disease = doc['disease']
            for snpdoc in doc['snps']:
                snp = snpdoc['snp']
                eth =  snpdoc['ethnicity'].lower()
                study =  snpdoc['study_in'].lower()
                gt = snpdoc['gt']
                gtrev = self.reverse(snp,gt,snpdoc['chromosome'] )
#                 print disease, study, self.study
                if eth != self.ethn and study != self.ethn and study not in self.study[disease]:
                    if 'rs' not in snp and snp not in unique_snp: # haplotype
                        unique_snp, list_snp = self.get_snp_haplotypes(snp, unique_snp, disease, study)
                    if snp not in unique_snp:
                        unique_snp[snp] = {}
                        unique_snp[snp]['chr'] = snpdoc['chromosome']
                        unique_snp[snp]['gene'] = snpdoc['gene']
                        unique_snp[snp]['pos'] = snpdoc['relative_position']
                        unique_snp[snp]['disease'] = {}
                        if 'rs' not in snp:
                            unique_snp[snp]['list'] = list_snp
                    if disease not in unique_snp[snp]['disease']:
                        unique_snp[snp]['disease'][disease] = {}  
                    add = True
                    if study not in unique_snp[snp]['disease'][disease]:
                        for eth_not_blue in self.study[disease]:
                            if eth_not_blue in unique_snp[snp]['disease'][disease]:
                                add = False
                    if add == True and study not in unique_snp[snp]['disease'][disease] and self.ethn not in unique_snp[snp]['disease'][disease] and study not in self.study[disease]  :
                        unique_snp[snp]['disease'][disease][study] ={}
                        unique_snp[snp]['disease'][disease][study]['gt'] ={}
                    if study in unique_snp[snp]['disease'][disease] : 
                        gt_list = unique_snp[snp]['disease'][disease][study]['gt'] 
                        if gt not in gt_list and gtrev not in gt_list or ((gt in gt_list or gtrev in gt_list) and 'haplotype' in unique_snp[snp]):  
                            unique_snp[snp]['disease'][disease][study]['gt'][gt] = '%s_%s' % (math.exp(snpdoc['or']),study ) 
        return unique_snp
    #END OF get_snp_info_blue
    
    ### generate_genotype_data
    def generate_genotype_data(self,outfile,snp_info):
#         print(' Function: generate_genotype_data ' ) 
        random.seed() 

        hap_doc = {}
        newdoc = {}
        data = {}
        newdoc['date'] = DATE
        newdoc['snps'] = []
        for snp in collections.OrderedDict(sorted(snp_info.items())):
            if 'haplotype' not in snp_info[snp]:
                numd = 0
                for disease in snp_info[snp]['disease']:
                    
                    if disease not in self.risk:
                        self.risk[disease] = {}
                        self.risk[disease]['count'] = 0
                        self.risk[disease]['OR'] = 0
                    if numd == 0:
                        numd += 1
                        tmp = copy.deepcopy(snp_info[snp])
#                         print snp, disease , tmp['disease']
                        tmp['disease'][disease] = self.get_mean_or(disease,tmp['disease'])
                        if  'rs' in snp:
                            genotype = self.get_genotype(disease,tmp)
                        else:
                            genotype = self.get_genotype(disease,tmp)
                            hap_doc[snp]=snp_info[snp]
                if 'rs' in snp:
                    self.output_marker_line( outfile, snp, genotype, snp_info[snp], newdoc, data)  
                else:                  
                    data[snp] = genotype         
        ### get genotype for snp in hapltypes
        hap_snp = {}
        for hap in hap_doc:
            newdoc, data,hap_snp = self.get_haplotype_snp_genotype(hap, snp_info, newdoc, data,outfile,hap_snp)

        return newdoc, data
    ### END of generate_genotype_data

    def get_mean_or (self,disease, snp_info):
#         print(' Function: get_mean_or') 
        tmp = {"caucasian":{'gt':{}}}
        i = 0
        for d in snp_info:
#             print d
            for e in snp_info[d]:
#                 print e
                for gt in snp_info[d][e]["gt"]:
                    OR = self.get_unicod_OR(snp_info[d][e]["gt"][gt])
                    if i ==0:
                        tmp["caucasian"]['gt'][gt] = 0
                    elif gt not in tmp["caucasian"]['gt']:
                        gt = gt[1]+gt[0]
                    if type(OR) is unicode and len(snp_info) ==1 and len(snp_info[d]) == 1:
                        tmp["caucasian"]['gt'][gt] = OR
                    elif type(OR) is not unicode:
                        tmp["caucasian"]['gt'][gt] += OR
#                     else: 
#                         print "ERROR ?",d, len(snp_info) , len(snp_info[d]) 
#                     print gt,d, OR, tmp["caucasian"]['gt']
                i += 1
        for gt in tmp["caucasian"]['gt']:
            if i >0 and type(tmp["caucasian"]['gt'][gt]) is not unicode:
                tmp["caucasian"]['gt'][gt]/= i
        return tmp

    ### output_marker_line
    def output_marker_line(self, outfile, marker, genotype, snp_info, newdoc, data):
#        print(' Function: output_marker_line', marker, genotype, snp_info) 
        ### randomly reverse the genotype for testing purposes
        genotype = self.reverse_genotype(marker,genotype,snp_info['chr'])
        newsnp = {}
        newsnp['gt'] = genotype
        newsnp['id'] = marker
        newdoc['snps'].append(newsnp)
        data[marker] = genotype
        outfile.write('%s\t%s\t%s\t%s\t%s\t' % (marker,str(snp_info['chr']) ,snp_info['pos'],snp_info['gene'],genotype))
        
        ### Write/calculate genetic risk per disease for each snp
        self.write_genotype_risk(outfile,marker,genotype,snp_info['disease'], '\t')
        return newdoc, data
    ### END of output_marker_line

    ### get_haplotype_snp_genotype
    # get SNP genotype from their haplotypes
    def get_haplotype_snp_genotype(self, hap, snp_info, newdoc, data,outfile,hap_snp):
#        print(' Function: get_haplotype_snp_genotype', hap,snp_info[hap]['list']) 
        num = 0
        for disease in snp_info[hap]['disease']:
            if num == 0:
                #### populate the allele for this haplotype
                if data[hap] != 'NA':
                    alleles = data[hap].split('/')
                    new_alleles = ['','']
                    pos = 0
                    for snp in snp_info[hap]['list']:
                        if snp in data:
                            a1 = hap_snp[snp][0]
                            a2 = hap_snp[snp][1]
                            new_alleles[0] = '%s%s' % (new_alleles[0], a1)
                            new_alleles[1] = '%s%s' % (new_alleles[1], a2)
                        else :
                            a1 = alleles[0][pos]
                            a2 = alleles[1][pos]
                            hap_snp[snp] = '%s%s' % (a1, a2)
                            
                            new_alleles[0] = '%s%s' % (new_alleles[0], a1)
                            new_alleles[1] = '%s%s' % (new_alleles[1], a2)
                            genotype = '%s%s' % (a1,a2)
                            self.output_marker_line( outfile, snp, genotype, snp_info[snp], newdoc, data)
                        pos += 1
                    hap_geno = '%s/%s' % (new_alleles[0],new_alleles[1])
                    self.output_marker_line( outfile, hap, hap_geno, snp_info[hap], newdoc, data, )      
                else:
                    hap_geno = 'NA'
                    for snp in snp_info[hap]['list']:
                        if snp not in data:
                            genotype = "NA"
                            self.output_marker_line( outfile, snp, genotype, snp_info[snp], newdoc, data)
                    self.output_marker_line( outfile, hap, hap_geno, snp_info[hap], newdoc, data, )               
            else: break
        return newdoc, data,hap_snp
    ### END of get_haplotype_snp_genotype


    ### write_genotype_risk
    def write_genotype_risk(self,outfile,snp,genotype,snp_info, type_file):   
#            print(' Function: write_genotype_risk', snp,genotype,snp_info) 
        for disease, value in sorted(self.bible_mongo.iteritems(), key=lambda (k,v): (v,k)):
            ### get OR and add to disease OR
            gtype_risk =''
            if disease in snp_info:
                for study in snp_info[disease]:
                    g = genotype
                    if re.match(re.compile('rs*'),snp) is not None and len(genotype) > 1:
                        g = genotype[1]+genotype[0]
                    elif len(genotype.split('/')) >1:
                        g = genotype.split('/')
                        g = '%s/%s' % (g[1],g[0])
                    if len(genotype) > 1 and (genotype == '--' or g == '--'):
                        g = 'NA'
                    if genotype in snp_info[disease][study]['gt']:
                        tmp =  (self.GenotypeRisk(disease,snp_info[disease][study]['gt'][genotype]))
                    elif g in snp_info[disease][study]['gt']:
                        tmp = (self.GenotypeRisk(disease,snp_info[disease][study]['gt'][g]))
                    else:
                        tmp =''
                    if gtype_risk !='' and tmp != '':
                        gtype_risk = '%s*%s' %(gtype_risk , tmp)
                    elif gtype_risk == '':
                         gtype_risk = tmp
            else:
                gtype_risk =''
            outfile.write('%s%s' % (gtype_risk,type_file))
        outfile.write("\n")
    ### END of write_genotype_risk

    ### get_genotype
    def get_genotype(self,disease,snp_info):    
#         print(' Function: get_genotype', disease, snp_info) 
        for study in snp_info['disease'][disease]:
            break    
        genotypes = {}
        OR_list = []
        for gt in snp_info['disease'][disease][study]['gt']:
            OR = snp_info['disease'][disease][study]['gt'][gt]
            self.get_unicod_OR(OR)  
            if (len(gt) == 1 or gt == 'NA') and  ((snp_info['chr'] == 'X' and self.gender == 'M') or snp_info['chr'] == 'Y'):
                genotypes[OR]=gt
                OR_list.append(OR)
                break
            elif (len(gt) == 2 ) and  ((snp_info['chr'] == 'X' and self.gender == 'F')):
                OR_list.append(OR)
                genotypes[OR]=gt
            elif (snp_info['chr'] == 'Y' or (len(gt) ==1  and snp_info['chr'] == 'X'))  and self.gender == 'F':
                ''' do nothin '''
            else :
                OR_list.append(OR)
                genotypes[OR]=gt
        ### choose genotype for the member
        OR = 0
        gt = "NA"
        if self.genetic == 'NA':
            gt = "NA"
        elif self.genetic[0] == 'R':
            if self.genetic[1] == 'H':
                new = self.find_higher(sorted(OR_list), 'R')
                gt = genotypes[new[random.randint(0,len(new))]]
#                 print 'RH', gt, new,genotypes
            else:
                new = self.find_lower(sorted(OR_list), 'R')
                gt = genotypes[new[random.randint(0,len(new))]]
#                 print 'RL', gt, new,genotypes
        elif self.genetic == 'H':
            new = self.find_higher(sorted(OR_list), 'H')
            gt = genotypes[new[random.randint(0,len(new))]]
#             print 'H', gt, new,genotypes
        else :
            new = self.find_lower(sorted(OR_list), 'L')
            gt = genotypes[new[random.randint(0,len(new))]]
#             print 'L', gt, new,genotypes
        return gt
    ### END of get_genotype
    
    def get_unicod_OR(self,OR):
        if type(OR) is unicode:
            if len(OR.split("_"))>1:
                OR = float(OR.split("_")[0])
                OR = math.log(OR)
            elif len(OR.split(": "))>1:
                OR = OR.split(": ")[1]
                OR = OR.split("'")[0]
            elif  len(OR.split("'"))>1:
                OR = OR.split("'")[1]
        return OR 
 
    ''' find higher '''
    def find_higher(self, OR_list, tag):
        if 'Similar' in OR_list:
            if 'Faster' in OR_list:
                return  ['Faster']
            elif 'Higher' in OR_list:
                return  ['Higher']
            else:
                return  ['Similar']
        elif 1 in OR_list:
            print 'PROBLEM find_higher?'
            OR_list = self.get_higher(1,OR_list)
        elif 0 in OR_list:
            if tag == 'R':
                OR_list = self.get_higher(0,OR_list)  
            if tag =="H":
                OR_list = [max(OR_list)]
        return OR_list
            
#      Slower     Similar     Faster
#  Lower     Similar     Higher
    ''' find lower '''
    def find_lower(self, OR_list, tag):
        if 'Similar' in OR_list:
            if 'Slower' in OR_list:
                return  ['Slower']
            elif 'Lower' in OR_list:
                return  ['Lower']
            else:
                return  ['Similar']
        elif 1 in OR_list:
            print 'PROBLEM find_lower?'
            OR_list = self.get_lower(1,OR_list)
        elif 0 in OR_list:
            if tag == 'R':
                OR_list = self.get_lower(0,OR_list)  
            if tag =="L":
                OR_list = [min(OR_list)]
        return OR_list
        
    def get_higher(self,lim, OR_list):
        new_OR = []
        for OR in OR_list:
            if OR>lim:
                new_OR.append(OR)
        if len(new_OR) == 0:
            new_OR.append(lim)
        return new_OR

    def get_lower(self,lim, OR_list):
        new_OR = []
        for OR in OR_list:
            if OR<lim:
                new_OR.append(OR)
        if len(new_OR) == 0:
            new_OR.append(lim)
        return new_OR
    
    ### reverse_genotype
    def reverse_genotype(self,snp,genotype,chr_num):    
        genotype = self.reverse(snp,genotype,chr_num)

        return genotype
    ### END of reverse_genotype

    ### reverse
    def reverse(self,snp,genotype,chr_num):    
#        print(' Function: reverse - snp: %s, genotype: %s' %(snp, genotype)) 
        rev = 1
        #rev = random.randint(0,2)
        if genotype != 'NA' and chr_num != 'X' and chr_num != 'Y':
            if re.match(re.compile('rs*'),snp) is not None:
                if genotype[1] != genotype[0] and rev > 0:
                    genotype = genotype[1]+genotype[0]
            else:
                g = genotype.split('/')
                if g[1] != g[0] and rev > 0:
                    genotype = '%s/%s' % (g[1],g[0])
#        print genotype
        return genotype
    ### END of reverse
    
    ### GenotypeRisk
    ###
    def GenotypeRisk(self,disease, odd_ratio):
#        print(' Function: GenotypeRisk')   
        ''' find the risk of the specific genotype for the speficied ethnicity and diffeent disease ids'''
        # snp is found for that disease:
        self.risk[disease]['count'] += 1
        # calculate the risk based on the Genotype
        try: 
            gtype_risk = round(math.exp(odd_ratio),4)
            self.risk[disease]['OR'] += odd_ratio  
        except: 
            gtype_risk = odd_ratio
        return gtype_risk 
    ### END of GenotypeRisk
     
    ### update_genomsequence
    ###
    def update_genomsequence(self,gid, realgid, newdoc):
#         print(' Function: update_genomsequence', gid, realgid)   

        finddoccursor = GENOPHEN30.sasha.authentication.find({'static.username':realgid},{"static.gid":1}) 
        i =0
        for doc in finddoccursor  :
            
            true_gid = doc['static']['gid']
            i +=1
        if i == 1:
            GENOPHEN30.genomsequence.remove({"gid":true_gid}) 
            newdoc['gid'] = true_gid
            
            GENOPHEN30.genomsequence.insert(newdoc) 
            newdoc.pop('_id')
            
        else:
            print 'PROBLEM update_genomsequence', gid, realgid
        return newdoc
    ### END of update_genomsequence
     
    ### output_summary
    ###
    def output_summary(self,outfile,t_f):
#         print(' Function: output_summary')
#         outfile.write('\n\n#<---------------------------------------------------------------->\n')
        outfile.write('#ETHNICITY : ' + self.ethn + "\n")
        outfile.write("##\tDiseaseID\tnum_SNPs\tlogOR\tOR\tstrenght\n"  )

        for disease, value in sorted(self.bible_mongo.iteritems(), key=lambda (k,v): (v,k)):
            if disease in self.risk:
                log_or = self.risk[disease]['OR']
                if log_or < 0:
#                     msg = round(1/math.exp(log_or),2)
                    msg = 'protective'
#                     outfile.write("###### DiseaseID%s%s%snum SNPs%s%d%srisk log%s%s%s%s%s%s %s##########\n"  
#                                   % (t_f,disease,t_f,t_f, self.risk[disease]['count'],t_f,t_f,log_or,t_f,msg,t_f,round(1/math.exp(log_or),2),t_f))
                else :
#                     msg = round(math.exp(log_or),2)
                    msg = 'risk'
#                     outfile.write("###### DiseaseID,%s,num SNPs,%d,risk log,%s,%s,%s ,##########\n"  % (disease, self.risk[disease]['count'],log_or,msg,round(math.exp(log_or),2)))
#                 outfile.write("###### DiseaseID%s%s%snum SNPs%s%d%srisk log%s%s%s%s%s%s %s##########\n"  
#                                   % (t_f,disease,t_f,t_f, self.risk[disease]['count'],t_f,t_f,log_or,t_f,msg,t_f,round(1/math.exp(log_or),2),t_f))
                outfile.write("##\t%s\t%s\t%s\t%s\t%s\n"  
                                % (disease, self.risk[disease]['count'],log_or,round(math.exp(log_or),2),msg))
    
    
#         outfile.write('#<---------------------------------------------------------------->\n')
    ### END of output_summary
    
def getInput(): 
#     print(' Function: getInput')
    input = open("%s/%s" % (DIRECTORY,"memberGenetic.csv"), 'rU')
    reader = csv.reader(input)
    members = {}
    header = []
    i =0

    for line in reader:
        if i == 0:
            header = line
        else:
            j =0
            username = ""
            for h in header:
                if h == 'username':
                    username= line[j]
                    members[username] = {}
                elif username != "":
                    members[username][h] = line[j]
                j+=1
        i += 1
    return members
### main
###     
if __name__ == "__main__":

##############  female caucasian ################
    ethn = 'african_american'
    g = 'F'
    gid='aut_601'

##############  FeMale Asian ################

#    ethns = ['Caucasian', 'Hispanic', 'native_american', 'Asian', 'african_american', 'pacific_islander']

    members = getInput()
    for m in members:
        if len(members[m]["ethnicity"].split(" ")) >1:
            words = members[m]["ethnicity"].split(" ")
            members[m]["ethnicity"] = words[0]+"_"+words[1]
        print m, members[m]
        member_genetics = UserGenome("%s" % m,members[m]["ethnicity"], 
                                     members[m]["gender"], members[m]["genetic"],m)
        if member_genetics.genetic != "NULL":
            member_genetics.generate_user_data()
    print('DONE')
