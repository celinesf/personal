#!/usr/bin/python
'''
Created on October 9, 2012
@author: celine

update 2.1 genomics documents (in json format)
- add real and relative position
- add gene name
- update/clean effect
'''
import json
import pymongo
import MySQLdb #@UnresolvedImport
import copy

import RelativePosition as RP
RP_class = RP.RelativePosition()

# tiger_dev
mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.239.40/admin"

#tiger_qa
#mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@54.245.229.194/admin"

#mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.140.237/admin"
CW = pymongo.Connection(mongo_cnxn_string, 27017);
GENOPHEN30 = CW["genophen30"]
DIR = '/Users/celine/Box Documents/Celine/ScriptInOut/mongodb-update'

'''
obtain real position from mongo dbsnp
'''
def get_snp_info(snp) :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    geneticsdb = cw["genetics"];
    chromosome = None 
    position = None
    finddoccursor= geneticsdb.dbsnp.find({"genome":"GRCh37.p5","rs":snp});            
    for c in finddoccursor:
        if c["chr_pos"] != '':
            chromosome = c["chr"] 
            position =  c["chr_pos"]
    return chromosome, int(position)

'''
obtain gene of a SNP from snpor (sql)
'''
def get_gene(snp, dbname):
    ### SQL INIT
    db = MySQLdb.connect(host="127.0.0.1",user="root",passwd="",db="algobck")
    cursor = db.cursor()
    sql= " SELECT distinct  genenextbio FROM %s S where dbsnp ='%s' " % (dbname, snp)
    #print sql
    gene = None
    cursor.execute(sql)
    tab = cursor.fetchall()
    for s in tab:
        gene = s[0]
    cursor.close ()
    return  gene

'''
obtain position from human genome 36 of a SNP from SNPInfo (sql)
'''
def get_position36(snp):
    ### SQL INIT
    db = MySQLdb.connect(host="127.0.0.1",user="root",passwd="",db="algobck")
    cursor = db.cursor()
    sql= " SELECT position FROM SNPInfo S where dbsnp ='%s' " % (snp)
    #print sql
    position = None
    cursor.execute(sql)
    tab = cursor.fetchall()
    for s in tab:
        position = s[0]
    cursor.close ()
    return  int(position)

'''
obtain position from human genome 36 of a SNP from SNPInfo (sql)
'''
def get_position36_allele(snp):
    ### SQL INIT
    db = MySQLdb.connect(host="127.0.0.1",user="root",passwd="",db="algobck")
    cursor = db.cursor()
    sql= " SELECT position, MajorAllele, MinorAllele FROM SNPInfo S where dbsnp ='%s' " % (snp)
    #print sql
    position = None
    cursor.execute(sql)
    tab = cursor.fetchall()
    for s in tab:
        position = s[0]
        MajorAllele = s[1]
        MinorAllele = s[2]
    cursor.close ()
    return  int(position), MajorAllele, MinorAllele

'''
update the position and gene from snpor(sql) and dbsnp (mongo)
'''
def update_info(genomics, dbname, update):
    genomics2 = []
    for chemical in genomics:
        if update == None or chemical['display_name'] in update:
            chemical2 = {}
            chemical2['display_name'] = chemical['display_name'] ### drug or food name
            chemical2['snp'] = []
            for snp in chemical['snp']:
                eth = snp['genotype'][0]['ethnicity_studied']
                if update == None or (snp['name'] in update[chemical2['display_name']] and eth in update[chemical2['display_name']][snp['name']]):
                    ### obtain real positions and gene name
                    snp2 = copy.deepcopy(snp)
                    snp2['chromosome_number'], snp2['position']  = get_snp_info(snp['name']) 
                    snp2['position36']  = get_position36(snp['name']) 
                    snp2['gene'] =  get_gene(snp['name'], dbname)
                    
                    ### obtain relative data
                    snp2['relative_position'] = RP_class.get_relative_position(snp2['chromosome_number'], snp2['position'])
                    #print chromosome, int(position), snp2['relative_position']
                    for genotype in snp2['genotype']:
                        ### clean drug effect
                        if len(genotype['effect'].split(' to most people'))>1:
                            genotype['effect'] = genotype['effect'].split(' to most people')[0]
                        if len(genotype['effect'].split(' than most people'))>1:
                            genotype['effect'] = genotype['effect'].split(' than most people')[0]
                            
                        ### update nutri effect
                        if dbname == 'NutrioG':
                            genotype['effect'] = update[chemical['display_name']][snp['name']][eth][genotype['genotype']]
                chemical2['snp'].append(snp2)
            genomics2.append(chemical2)
    return genomics2

''' 
update effects from file that Malekeh updated 
'''
def get_updated_nutrio():
    nutri = {}
    inputfile = open('nutri_malekeh.csv','r')
    for line in  inputfile:
        snp = line.split(',')[3]
        if len(snp.split('rs')) >1 :
            chem = line.split(',')[2]
            eth = line.split(',')[4]
            geno = line.split(',')[5]
            effect = line.split(',')[6]
            if chem not in nutri:
                nutri[chem] = {}
            if snp not in nutri[chem]:
                nutri[chem][snp] = {}
            if eth not in nutri[chem][snp]:
                nutri[chem][snp][eth] = {}
            nutri[chem][snp][eth][geno] = effect
    inputfile.close()
    return nutri

'''
update hapid
'''
def update_hapid() :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    geneticsdb = cw["genetics"];
    geneticsdb.IDHap.rep.remove({});  
    finddoccursor = geneticsdb.IDHap.find({});    
            
    for c in finddoccursor:
        print c['dbHAP'], c['_id']
        newdoc = c
        newdoc['chromosome'] , newdoc['position'] = get_snp_info(c['dbSNP'])
        newdoc['position36'] , newdoc['MajorAllele'] ,newdoc['MinorAllele'] = get_position36_allele(c['dbSNP'])
        newdoc['relative_position'] = RP_class.get_relative_position(newdoc['chromosome'] , newdoc['position'])
  
        geneticsdb.IDHap.rep.insert(newdoc)
        
'''
update update_hapid_place
'''
def update_hapid_place() :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    geneticsdb = cw["genetics"];
    finddoccursor = geneticsdb.IDHap.rep.find({});    
    hapinfo = {}
    happlace = {}
    for c in finddoccursor:
        hapname = c['dbHAP']
        snp = c['dbSNP']
        position = c['position']
        if hapname not in hapinfo:
            happlace[hapname] = {}
            print 'first',snp
            c['place2'] = 1
            hapinfo[hapname] = []
            hapinfo[hapname].append(c)
            happlace[hapname][snp] = {}
            happlace[hapname][snp]['place'] = 1
        else:
            newhap = []
            counter = 1
            done = False
            print 'starting',snp
            for snpinfo  in hapinfo[hapname]:
                if position < snpinfo['position'] and done == False:
                    print 'found place', hapname, snp,position ,snpinfo['dbSNP'],snpinfo['position']
                    c['place2'] = counter
                    newhap.append(c)
                    happlace[hapname][snp] = {} ## new
                    happlace[hapname][snp]['place'] = counter 
                    counter += 1
                    done = True
                print 'adding', hapname, snp,position ,snpinfo['dbSNP'],snpinfo['position'], counter
                snpinfo['place2'] = counter
                newhap.append(snpinfo)
                if snpinfo['dbSNP'] not in happlace[hapname]:
                    happlace[hapname][snpinfo['dbSNP']] = {}
                happlace[hapname][snpinfo['dbSNP']]['place'] = counter  ## update
                counter += 1
            if position > hapinfo[hapname][len(hapinfo[hapname])-1]['position']:
                print 'last place', hapname, snp,position ,snpinfo['dbSNP'],snpinfo['position']
                c['place2'] = counter
                newhap.append(c)   
                happlace[hapname][snp] = {} ## new
                happlace[hapname][snp]['place'] = counter 
                counter += 1
            print counter        
            hapinfo[hapname] = copy.deepcopy(newhap)
    #geneticsdb.IDHap.rep2.remove({});  
    finddoccursor = geneticsdb.IDHap.rep.find({});    
    for c in finddoccursor:
        newdoc = c
        hapname = c['dbHAP']
        snp = c['dbSNP']
        if happlace[hapname][snp]['place'] != c['place']:
                print 'lIssue', hapname, snp,c['place'] ,happlace[hapname][snp]['place'],c['position'],c['position36']
        newdoc['place'] = happlace[hapname][snp]['place']
        #geneticsdb.IDHap.rep2.insert(newdoc)
        geneticsdb.IDHap.rep.update({'_id':c['_id'] },newdoc)
   
#    for hapname in hapinfo:
#        for snp in hapinfo[hapname]:
#            if snp['place'] != snp['place2']:
#                print 'lIssue', hapname, snp['dbSNP'],snp['place'] ,snp['place2'],snp['position'],snp['position36']
#    

'''
update get_hap_info
'''
def get_hap_info(hap) :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    geneticsdb = cw["genetics"];
     
    finddoccursor = geneticsdb.IDHap.rep.find({'dbHAP':hap,'place':1});                
    for c in finddoccursor:
        return c['position'],c['position36'], c['relative_position']
        
    return None,None,None
        

'''
update update_haplotype_doc
'''
def update_haplotype_doc() :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    genophen30 = cw["genophen30"];
     
    finddoccursor = genophen30.genetics.rep.find({});     
              
    for c in finddoccursor:
        print c['disease']
        resultsnplist = [];        
        for snp in c['snps']:
            if len(snp['snp'].split('rs'))<=1:
                #print snp['snp'],
                newhap = snp
                newhap['position'],newhap['position36'], newhap['relative_position'] = get_hap_info(snp['snp'])
                #print newhap['position'],newhap['position36'], newhap['relative_position'] 
                resultsnplist.append(newhap);
            else:
                resultsnplist.append(snp);
        
        genophen30.genetics.rep.update({'disease':c['disease']},{"$set":{"snps":resultsnplist}});  

'''
update change_gene_genetics
'''
def change_gene_genetics() :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    genophen30 = cw["genophen30"]; 
    finddoccursor = genophen30.genetics.rep.find({});               
    for c in finddoccursor:
        print c['disease']
        resultsnplist = [];        
        for snp in c['snps']:
            newhap = snp            
            if len(snp['gene'].split('_')) > 1:
                counter = 0
                gene1 = ''
                for gene in snp['gene'].split('_'):                    
                    if counter == 0:
                        newhap['gene'] = gene
                    elif counter > 0 and gene1.lower() == 'hcg':
                        newhap['gene']  = '%s_%s' %(newhap['gene'] , gene)
                    elif counter > 0 and gene1.lower() != 'hcg':
                        newhap['gene']  = '%s,%s' %(newhap['gene'] , gene) 
                    else: 
                        print "WHERE AM I?"    
                    gene1 = gene               
                    counter += 1
            resultsnplist.append(newhap);  
        genophen30.genetics.rep.update({'disease':c['disease']},{"$set":{"snps":resultsnplist}});  

'''
update change_gene_genomics
'''
def change_gene_genomics() :
    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@ec2-50-112-118-154.us-west-2.compute.amazonaws.com:25565/admin"
    cw = pymongo.Connection(mongo_cnxn_string, 25565);
    genophen30 = cw["genophen30"]; 
    finddoccursor = genophen30.genomics.find({});  
    
    print 'here'             
    for c in finddoccursor:
        if c['gid'] == 'drugs' or c['gid'] == 'nutrition' :
            #genophen30.genomics.rep.insert(c);  
            print c['gid'] 
            resultlist = [];  
            for chem in c[c['gid']]:
                newchem = {}
                newchem['display_name'] = chem['display_name']
                newchem['snp'] = []
                for snp in chem['snp']:
                    newsnp = snp 
                               
                    if len(snp['gene'].split('_')) > 1:
                        counter = 0
                        gene1 = ''
                        for gene in snp['gene'].split('_'):                    
                            if counter == 0:
                                newsnp['gene'] = gene
                            elif counter > 0 and gene1.lower() == 'hcg':
                                newsnp['gene']  = '%s_%s' %(newsnp['gene'] , gene)
                            elif counter > 0 and gene1.lower() != 'hcg':
                                newsnp['gene']  = '%s,%s' %(newsnp['gene'] , gene) 
                            else: 
                                print "WHERE AM I?"    
                            gene1 = gene               
                            counter += 1
                        newchem['snp'].append(newsnp);
                    else:
                        newchem['snp'].append(snp);  
                resultlist.append(newchem)
           
            genophen30.genomics.rep.update({'gid':c['gid']},{"$set":{c['gid']:resultlist}});  
        
def get_lab_ranges():
    f1 = open('old_lab_ranges_out.json','r')
    f2 = open('medical_profile.json','r')
    data1 = json.loads(f1.read())
    data2 = json.loads(f2.read())

    dic = []
    for info in data1['dic']:
        dic.append(info['v'])
        if info['v'] not in data2:
            print (info['v'])
            
    for info in data2:
        if info not in dic:
            print info        
        
#    mongo_cnxn_string = "mongodb://adm1n:adm1npw0rd@50.112.142.17/admin"
#    CW = pymongo.Connection(mongo_cnxn_string, 27017);
#    GENOPHEN30 = CW["genophen30"]
#    
#    labform = GENOPHEN30.general.find({'name':'lab_form'})
#    for doc in labform:
#        print json.dumps(doc['list'])
#        for lab in doc['list']:
#            print lab['c'],lab['d'], lab['lower'], lab['upper']


def get_genetic_rep(disease):
    finddoccursor = GENOPHEN30.genetics.rep.find({'disease':disease})  
    
    for document in finddoccursor:
        newdoc = {'disease': document['disease'],'date':"02-20-2013","version":"3.0", 'snps': []}
#        print document['disease']
        pop = {}
        study = {}
        for snps in document['snps']:
            pop[snps['ethnicity']] = 0
            study[snps['study_in']] = 0
            
        print json.dumps(document['snps'], indent=4)
#        print pop
#        print study
#        outfile = open("%s/%s" %(DIR,document['disease']),'w')
#        outfile.write(json.dumps(document['snps'], indent=4))
       
def change_genetic_rep(filename):
    infile = open('%s/%s' %(DIR, filename),'r')
    data = json.loads(infile.read())
    newdoc = {'disease': filename,'date':"02-22-2013","version":"3.0", 'snps': data}
    print json.dumps(newdoc, indent=4)
    GENOPHEN30.genetics.rep.update({'disease': filename},newdoc)

################ MAIN ###############
if __name__ == "__main__":
    ### drugs 
#    print 'updating drugs'
#    inputfile = open('drugs.json','r')
#    drugs = json.loads(inputfile.read())
#    inputfile.close() 
#    drugs['drugs'] = update_info(drugs['drugs'], 'Pharmaco', None)
#    output = open('drugs-udated.json','w')
#    output.write(json.dumps(( drugs), indent=4))
#    output.close()
    ### Nutrition 
#    print 'updating nutrition'
#    updated_nutrio = get_updated_nutrio()
#    inputfile = open('nutrition.json','r')
#    nutrition = json.loads(inputfile.read())
#    inputfile.close()
#    nutrition['nutrition'] = update_info(nutrition['nutrition'], 'NutrioG',updated_nutrio)
#    output = open('nutrition-udated.json','w')
#    output.write(json.dumps(( nutrition), indent=4))
#    output.close()
    
    ### update_hapid  
#    print 'updating hapid'
#    update_hapid()
#    
#    update_hapid_place()

    #update_haplotype_doc()


    ### change gene names
    #change_gene_genetics()
    #print 'I am done changing gene in genetics'
    
    ### get min and max for all lab work
#    get_lab_ranges()
    
#    change_gene_genomics()
    
    change_genetic_rep('ulcerative_colitis')
    get_genetic_rep('ulcerative_colitis')
    
    print 'I am done changing gene in genomics'

