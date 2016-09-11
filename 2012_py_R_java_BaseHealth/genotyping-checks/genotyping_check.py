#!/usr/bin/env python
'''
Created on dec, 12 2012
@author: celine

check the presence/ absence of 2.1 genophen yes  and 3.0 SNPs in Affi and illumina chip
'''


import pymongo
import re


######## dev
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.142.17/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genophen30"]

class UserGenome():
    ''' to generate randome Genotypes for an imaginary individual '''
     
    def __init__(self,filename , ethn, gen, gid):
        self.gender     = gen
        self.ethn       = ethn.lower()
        self.id         = filename
        self.gid         = gid
        self.risk       = {}
        self.snp_info   = {}
        self.bible_mongo      = {'type_2_diabetes':1, 'colorectal_cancer':4, 'alzheimer':5,'heart_disease':6, 'stroke':7,'depression':8,
                                 'hypertension':9,'anxiety':101,'backpain':102,'breast_cancer': 103,'crohns_disease':104, 
                           'hair_loss':105,'lung_cancer': 106,'migraines':107,'obesity':108, 'osteoporosis':109,'prostate_cancer':110,
                            'sleep_apnea':111,'ulcerative_colitis':112,'erectile_dysfunction':113, 
                            'Caffeine':200, 'Lactose':201, 'Alcohol':202, 'Carbohydrates':203, 'Calcium':204,'Vitamin D':205,
                            'Metformin':300,'Citalopram':301,'Fluoxetine':302,'Clopidogrel':303,'Desipramine':304,'Irbesartan':305 }

    
    ### END of generate_user_data
    
    ### generate_user_data
    ###
    def generate_user_data(self):     
        print(' Function: generate_user_data')   
        ### all the snps in our db
        snp_info = self.get_snp_info()
        snp_info = self.get_snp_info_blue(snp_info)
        snp_info = self.get_snp_genomics(snp_info)
        gid = '%s_%s' %(self.id, self.ethn+self.gender)
        outfile = open("%s.xls"%(gid), 'w')
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
        new_doc["genome_type"] = "sequencing"

        ### UNCHECK for SEQUENCING
#        sortder_new_doc = sorted(new_doc, key=itemgetter('chromosome', 'position'))
        self.update_genomsequence(gid, self.gid, new_doc)
        ### output test file
        self.output_summary(outfile)
        outfile.close()    
        
        ### write 23 & me data file
        self.write_23andme_file(gid,snp_info,data, 36)
        self.write_23andme_file(gid,snp_info,data, 37)
            
    ### END of generate_user_data
    
    ### write_23andme_file
    ###
    def write_23andme_file(self,gid,snp_info,newdoc, build):     
        print(' Function: write_23andme_file %s' % build)   
        ### init
        file23 = open("23&me_%s.txt" % build,"r")
        out23 = open("build_%s_%s.txt" % (build,gid),"w")

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
                    self.write_genotype_risk(out23,snp,newdoc[snp],snp_info[snp]['disease'])                    
#                else :
#                    out23.write(line)
        # add haplotype data
        out23.write("\n")
        for snp in snp_info:
            if len(snp.split("rs")) ==1 :
                out23.write('#%s\t%s\t%s\t%s\t' % (snp,str(snp_info[snp]['chr']) ,snp_info[snp]['pos'],newdoc[snp]))
                self.write_genotype_risk(out23,snp,newdoc[snp],snp_info[snp]['disease'])    

        self.output_summary(out23)
        file23.close()
        out23.close()
    ### END of write_23andme_file
    
    
    ### get_snp_genomics
    ###
    def get_snp_genomics(self, unique_snp):
        print(' Function: get_snp_genomics')
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
                            unique_snp[snp]['disease'][chem] = {}  
                            
                            unique_snp[snp]['disease'][chem]['gt'] ={}
                        for geno in snpdoc['genotype']:
                            unique_snp[snp]['disease'][chem]['gt'][geno['genotype']] = '%s - %s' % (geno['effect'], geno['ethnicity_studied']  )
                        #unique_snp[snp]['disease'] = sorted(unique_snp[snp]['disease'], key=unique_snp[snp]['disease'].get) 
        return unique_snp
    #END OF get_snp_genomics


    ### get_snp_info
    ###
    def get_snp_info(self):
        print(' Function: get_snp_info')
        finddoccursor = GENOPHEN30.genetics.rep.find({'gt':{'$ne':'NA'}})
        unique_snp = {}
        for doc in finddoccursor  :
            disease = doc['disease']
            for snpdoc in doc['snps']:
                snp = snpdoc['snp']
                eth =  snpdoc['ethnicity'].lower()
                gt = snpdoc['gt']
                if eth == self.ethn:
                    try:
                        unique_snp[snp]
                    except:
                        unique_snp[snp] = {}
                        unique_snp[snp]['chr'] = snpdoc['chromosome']
                        unique_snp[snp]['gene'] = snpdoc['gene']
                        unique_snp[snp]['pos'] = snpdoc['relative_position']
                        unique_snp[snp]['disease'] = {}
                        ########## TO DO ADD DESCRIPTION?
                    try:
                        unique_snp[snp]['disease'][disease]
                    except:
                        unique_snp[snp]['disease'][disease] = {}  
                        unique_snp[snp]['disease'][disease]['gt'] ={}
                    unique_snp[snp]['disease'][disease]['gt'][gt] = snpdoc['or']   
#                    unique_snp[snp]['disease'] = sorted(unique_snp[snp]['disease'], key=unique_snp[snp]['disease'].get) 
       
        return unique_snp
    #END OF get_snp_info


    ### get_snp_info_blue
    ###
    def get_snp_info_blue(self, unique_snp):
        print(' Function: get_snp_info_blue')
        #logging.debug(' Function: get_snp_info_blue: %s')
        finddoccursor = GENOPHEN30.genetics.rep.find({'ethnicity':{'$ne':self.ethn}})
        for doc in finddoccursor  :
            disease = doc['disease']
            for snpdoc in doc['snps']:
                snp = snpdoc['snp']
                eth =  snpdoc['ethnicity'].lower()
                gt = snpdoc['gt']
                if eth != self.ethn:
                    if snp not in unique_snp:
                        unique_snp[snp] = {}
                        unique_snp[snp]['chr'] = snpdoc['chromosome']
                        unique_snp[snp]['gene'] = snpdoc['gene']
                        unique_snp[snp]['pos'] = snpdoc['relative_position']
                        unique_snp[snp]['disease'] = {}
                        ########## TO DO ADD DESCRIPTION?
                    if disease not in unique_snp[snp]['disease']:
                        unique_snp[snp]['disease'][disease] = {}  
                        unique_snp[snp]['disease'][disease]['gt'] ={}
                    if gt not in unique_snp[snp]['disease'][disease]['gt']:
                        unique_snp[snp]['disease'][disease]['gt'][gt] = '%s_%s' % (math.exp(snpdoc['or']),eth )
                #unique_snp[snp]['disease'] = sorted(unique_snp[snp]['disease'], key=unique_snp[snp]['disease'].get)           
        return unique_snp
    #END OF get_snp_info_blue
    
    ### generate_genotype_data
    def generate_genotype_data(self,outfile,snp_info):
        print(' Function: generate_genotype_data ' ) 
        random.seed(12345)
        newdoc = {}
        data = {}
        newdoc['date'] = DATE
        newdoc['snps'] = []
        for snp in collections.OrderedDict(sorted(snp_info.items())):

#        for snp in snp_info:
            numd = 0
            for disease in snp_info[snp]['disease']:
                try :
                    self.risk[disease]
                except:
                    self.risk[disease] = {}
                    self.risk[disease]['count'] = 0
                    self.risk[disease]['OR'] = 0
                if numd == 0:
                    genotype = self.get_genotype(disease,snp_info[snp])
                    numd += 1

            ### randomly reverse the genotype for testing purposes
            genotype = self.reverse_genotype(snp,genotype,snp_info[snp]['chr'])

            newsnp = {}
            newsnp['gt'] = genotype
            newsnp['id'] = snp
            newdoc['snps'].append(newsnp)
            data[snp] = genotype
            outfile.write('%s\t%s\t%s\t%s\t%s\t' % (snp,str(snp_info[snp]['chr']) ,snp_info[snp]['pos'],snp_info[snp]['gene'],genotype))
            
            ### Write/calculate genetic risk per disease for each snp
            self.write_genotype_risk(outfile,snp,genotype,snp_info[snp]['disease'])

        #from operator import itemgetter
        #newdoc['snps'] = sorted(newdoc['snps'], key=itemgetter('chromosome','position')) 
        return newdoc, data
    ### END of generate_genotype_data


    ### write_genotype_risk
    def write_genotype_risk(self,outfile,snp,genotype,snp_info):   
#        print(' Function: write_genotype_risk', snp,genotype,snp_info) 
        for disease, value in sorted(self.bible_mongo.iteritems(), key=lambda (k,v): (v,k)):
            ### get OR and add to disease OR
            if disease in snp_info:
                if re.match(re.compile('rs*'),snp) is not None and len(genotype) > 1:
                    g = genotype[1]+genotype[0]
                elif len(genotype.split('/')) >1:
                    g = genotype.split('/')
                    g = '%s/%s' % (g[1],g[0])
                if genotype in snp_info[disease]['gt']:
                    gtype_risk = self.GenotypeRisk(disease,snp_info[disease]['gt'][genotype])
                elif g in snp_info[disease]['gt']:
                    gtype_risk = self.GenotypeRisk(disease,snp_info[disease]['gt'][g])
                else:
                    gtype_risk =''
            else:
                    gtype_risk =''
            outfile.write('%s\t' % (gtype_risk))
        outfile.write("\n")
    ### END of write_genotype_risk

    ### get_genotype
    def get_genotype(self,disease,snp_info):    
#        print(' Function: get_genotype') 
        genotypes = []
        for gt in snp_info['disease'][disease]['gt']:
            if (len(gt) == 1 or gt == 'NA') and  ((snp_info['chr'] == 'X' and self.gender == 'M') or snp_info['chr'] == 'Y'):
                genotypes.append(gt)
            elif (len(gt) == 2 ) and  ((snp_info['chr'] == 'X' and self.gender == 'F')):
                genotypes.append(gt)
            elif snp_info['chr'] == 'Y' and self.gender == 'F':
                ''' do nothin '''
            else :
                genotypes.append(gt)
        ### choose genotype for the member
        return genotypes[random.randint(0,len(genotypes)-1)]
    ### END of get_genotype
    
    ### reverse_genotype
    def reverse_genotype(self,snp,genotype,chr_num):    
#        print(' Function: reverse_genotype') 
        rev = random.randint(0,2)
        if rev == 1 and genotype != 'NA' and chr_num != 'X' and chr_num != 'Y':
            if re.match(re.compile('rs*'),snp) is not None:
                genotype = genotype[1]+genotype[0]
            else:
                g = genotype.split('/')
                genotype = '%s/%s' % (g[1],g[0])
        return genotype
    ### END of reverse_genotype

    
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
        print(' Function: update_genomsequence')   
        ### add/update this data
        try:
            GENOPHEN30.genomsequence.remove({"gid":realgid}) 
#            finddoccursor = GENOPHEN30.genomsequence.find({"gid":realgid})
#            for doc in finddoccursor  :
#                 doc['gid'] = ""
#                 GENOPHEN30.genomsequence.update({"gid":realgid},doc)
#            print 'removed', gid
        except:
            print 'creating', gid
        GENOPHEN30.genomsequence.insert(newdoc)  
    ### END of update_genomsequence
     
    ### output_summary
    ###
    def output_summary(self,outfile):
        print(' Function: output_summary')
        outfile.write('\n\n#<---------------------------------------------------------------->\n')
        outfile.write('#ETHNICITY : ' + self.ethn + "\n")

        for disease, value in sorted(self.bible_mongo.iteritems(), key=lambda (k,v): (v,k)):
            if disease in self.risk:
                log_or = self.risk[disease]['OR']
                if log_or < 0:
                    msg = 'protective'
                    outfile.write("###### DiseaseID\t%s\tnum SNPs\t%d\trisk log\t%s\t%s\t%s \t##########\n"  % (disease, self.risk[disease]['count'],log_or,msg,round(1/math.exp(log_or),2)))
    
                else :
                    msg = 'risk'
                    outfile.write("###### DiseaseID\t%s\tnum SNPs\t%d\trisk log\t%s\t%s\t%s \t##########\n"  % (disease, self.risk[disease]['count'],log_or,msg,round(math.exp(log_or),2)))

        outfile.write('#<---------------------------------------------------------------->\n')
    ### END of output_summary
     
### main
###     
if __name__ == "__main__":
##############  female caucasian ################
    ethn = 'Caucasian'
    g = 'F'
    gid = "40b93e00-3dab-11e2-a25f-0800200c9a66"
    name = "celine"
    
#    gid = "e02c7c90-3dab-11e2-a25f-0800200c9a66"
#    name = "Soma"
#    
#    gid = "cab72b80-3dab-11e2-a25f-0800200c9a66"
#    name = "Archana"
#    
#    gid = "203e2fd0-3e48-11e2-a25f-0800200c9a66"
#    name = "Malekeh"
##############  Male caucasian ################
#    ethn = 'Caucasian'
#    g = 'M'
#    gid = "04977440-3dac-11e2-a25f-0800200c9a66"
#    name = "Pouria"
#    
#    gid = "0f826f30-3e48-11e2-a25f-0800200c9a66"
#    name = "Harnarayan"
#    
#    gid = "fd5432d0-3e47-11e2-a25f-0800200c9a66"
#    name = "Javid"
#    
#    gid = "f5b75170-3dab-11e2-a25f-0800200c9a66"
#    name = "Boris"
#    
#    gid = "dcf3930a-f7e9-4687-8858-f5edf81cacad"
#    name = "BorisDev1"
#
#    gid = "ec630330-3dab-11e2-a25f-0800200c9a66"
#    name = "Krishna"
#    
#    gid = "e663d360-3dab-11e2-a25f-0800200c9a66"
#    name = "Hossein"
#    
#    gid = "dab3bad0-3dab-11e2-a25f-0800200c9a66"
#    name = "Karim"
#    
#    gid = "9fa4c0a0-3dac-11e2-a25f-0800200c9a66"
#    name = "Anush"
#    
#    gid = "0b045a00-3dac-11e2-a25f-0800200c9a66"
#    name = "sunil"
##############  FeMale Asian ################ 
#    ethn = 'Asian'
#    g = 'F'
#    gid = "fc08bff0-3dab-11e2-a25f-0800200c9a66"
#    name = "Chantra"
#    
#    gid = "18c81050-3dac-11e2-a25f-0800200c9a66"
#    name = "Portia"
#  
#    gid = "139cf960-3dac-11e2-a25f-0800200c9a66"
#    name = "Thao"
#    
#    gid = "1f0b4370-3dab-11e2-a25f-0800200c9a66"
#    name = "Karen"
##############  FeMale Asian ################

#    ethns = ['Caucasian', 'Hispanic', 'native_american', 'Asian', 'african_american', 'pacific_islander']
#    ethns = ['Caucasian']
#    for ethn in ethns:
#    gender = ['M', 'F']
#     gender = ['F']
#        for g in gender:
    member_genetics = UserGenome("2dev_test_%s" % name,ethn, g, gid)
    member_genetics.generate_user_data()
    print('I am done with %s %s' % (ethn, g))
