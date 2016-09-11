#!/usr/bin/env python
"""
   Rank Risks
   08/11/12 - 0.0.2 : second version after the first was unceremoniously deleted
"""
__author__ = "Pouria Mojabi"
__copyright__ = "Copyright 2012, Genophen.com"
__maintainer__ = "Pouria Mojabi"
__email__ = "mojabi@genophen.com"
__status__ = "Test" #

"""
   Algorithm & Flow
   
   To collect all the samples and ors for snps per disease
   with the format that Celine's QC routine expects
"""


import MySQLdb as mdb
import json, codecs, csv
import pymongo

class Cauto():
    def __init__(self, d):
        self.disease = d
        self.ethmap = {"Caucasian" : "Freq_CEU", "NativeAmerican":None, "Hispanic":"Freq_MEX", "Asian":["Freq_CHB", "Freq_CHD", "Freq_JPT"] , "AfricanAmerican":"Freq_ASW" }
        
        self.snps    = None #unique SNPS of the disease
        self.celine  = {}
        
        # function calls
        self.unique_snps()
    
  
  
    def unique_snps(self):
        con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
        cur = con.cursor()
        cur.execute("select distinct dbSNP from SNPOR where DiseaseID=%s" % self.disease)
        rows = cur.fetchall()
        
        self.snps = [snp[0] for snp in rows]
        
        con.close()
        
    
    def unique_ethn_persnp(self, snp):
        con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
        cur = con.cursor()
        cur.execute("select distinct BroadEthnicity from SNPOR where DiseaseID=%s and dbSNP='%s'" % (self.disease, snp) )
        rows = cur.fetchall()
        
        be = [snp[0] for snp in rows]
        
        con.close()
        return be
  
    
    
    def unique_pmid_ethn_persnp(self, snp, eth):
        con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
        cur = con.cursor()
        cur.execute("select distinct PMID from SNPOR where DiseaseID=%s and dbSNP='%s' and BroadEthnicity='%s'" % (self.disease, snp, eth) )
        rows = cur.fetchall()
        
        be = [snp[0] for snp in rows]
        
        con.close()
        return be

        
        
        
    
    def papers_snp(self, snp):
        ''' to see how many papers have data for the snp'''
        con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
        cur = con.cursor()
        cur.execute("select distinct PMID from SNPOR where DiseaseID=%s and dbSNP='%s' " % (self.disease, snp) )
        rows = cur.fetchall()
        
        pmids = [snp[0] for snp in rows]
        
        con.close()
        return pmids
    
        
        
    def snp_hapmap(self, snp, ethn):
        ''' to see whether snp exists in hapmap and has frequency for the ethnicity '''
        con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
        cur = con.cursor()
        if self.ethmap[ethn]:
            if ethn == "Asian":
                cur.execute("select %s, %s, %s from SNPInfo where  dbSNP='%s' " % (self.ethmap[ethn][0], self.ethmap[ethn][1], self.ethmap[ethn][2],  snp) )
                rows= cur.fetchall()
                
            else:
                cur.execute("select %s from SNPInfo where  dbSNP='%s' " % (self.ethmap[ethn],  snp) )
                rows = cur.fetchall()
                if rows[0][0]:
                    return True
                else:
                    return False
        else:
            return False
        


        
        
    def get_snp_paper_sample(self, snp):
        ''' for a given snp, pubmedid find the sample'''
        con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
        cur = con.cursor()
        
        self.celine [snp] =  {'variant_id':snp, "nextbio_score":-1, "nextbio_pvalue":-1,
                       "broad_ethnicity" : { }
                        }
        
        
        sample = 0
        eths = self.unique_ethn_persnp(snp)
        for eth in eths:
            pmids = self.unique_pmid_ethn_persnp(snp, eth)
            eth_dict = {"ethnicity_id" : eth, "hapmap_freq" : self.snp_hapmap(snp, eth), "sample" : {} }
            
            for pmid in pmids:
                cur.execute("select * from SNPOR where DiseaseID=%s and dbSNP='%s' and PMID=%s and BroadEthnicity='%s'  " % (self.disease, snp, pmid, eth)  )
                rows = cur.fetchall()
                rows = [row for row in rows]
                
                row = rows[0]
                
                if len(rows) == 2:
                    model_id = "Allelic " + row[9]
                elif len(rows) == 3:
                    model_id = "Genotypic "+ str(row[9] )
                else:
                    model_id = "None"
                      
                if row[15]:
                    control_freq = True
                else:
                    control_freq = False  
            
                s = {"sample_id" : str(pmid), "sample_size" : row[18], "study_type":int(row[23].split()[3].split("=")[1]), "meta_sample":"None",
                     "association": {model_id:{"model_id":model_id, "control_freq": control_freq, "OR_data":{} } } }
                
                
                eth_dict["sample"][pmid] = s
                tmp = eth_dict
                
                or_data = {}
                
                for row in rows:
                    sample += 1
                    if row[11]:
                        ci = [ row[11].split("-")[0] ,   row[11].split("-")[1] ]
                    else:
                        ci = [-1, -1]                                        
                 
                    tmp['sample'][pmid]["association"][model_id]["OR_data"][ row[5] ] = {"OR_id" : str(row[5]),
                                        "odd_ratio":row[10], "CI" : ci,
                                        "pvalue" : row[12]}
                    
                    
                    
                
                
                self.finalize(snp, eth, tmp)
                    
            
        
 
    def finalize(self, snp, be, tmp):
        
        self.celine[snp]["broad_ethnicity"][be] = tmp
        
        #self.celine[snp] = variant
        
        
  
        
        
    def main(self):
        for snp in self.snps:
            sample = 0
            self.get_snp_paper_sample(snp)
            sample += 1
    
    
            
        
            
def compare_celine(file, did):
    json_data=codecs.open(file, encoding='utf-8')
    data = json.load(json_data, encoding="utf-8")
    
    ethns = data.keys()
    
    for ethn in ethns:
        snps={}
        for k, val in data[ethn].items():
            snps[val["LD_best_variant"] ]= val["variant_best_sample"]
        
        gyes_snps = get_ethn_gyes(ethn, did)
        
        print 
        print "="*50
        print "="*50
        print " "*20, ethn, " "*20
        for snp in gyes_snps:
            if snp not in snps.keys():
                
                print ">>>>>>>>>  present in Genophen2.1  but not here : ", snp, "<<<<<<<<<<"
            else:
                print "For SNP : ", snp, " Genophen2.1 PMID : ", get_pmid_snp(ethn, did, snp), " Celine PMID: ", snps[snp]
        
        #for snp in snps:
        #    if snp not in gyes_snps:
        #        print "present in Here but not Genophen yes : ", snp
        
        
            
    

def get_ethn_gyes(ethn, did):
    con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
    cur = con.cursor()
    cur.execute("select distinct dbSNP from SNPOR where Genophen like 'yes' and DiseaseID=%s and BroadEthnicity='%s'" % (did, ethn) )
    rows = cur.fetchall()
        
    snps = [snp[0] for snp in rows]
    
    con.close()
    return snps
    
    
def get_pmid_snp(ethn, did, snp):
    con = mdb.connect('192.168.2.108', 'journal', 'journal', 'AlgoPhen');
    cur = con.cursor()
    cur.execute("select distinct PMID from SNPOR where Genophen = 'yes' and DiseaseID=%s and BroadEthnicity='%s' and dbSNP='%s'" % (did, ethn, snp) )
    rows = cur.fetchall()
        
    pmid = rows[0][0]
    
    con.close()
    return pmid
    

def get_snpor_mongo():
    connection = pymongo.Connection('ec2-50-112-118-154.us-west-2.compute.amazonaws.com',25565)
    dbh = connection["admin"]
    dbh.authenticate("adm1n","adm1npw0rd")
    dbh = connection["genetics"] # admin"]
    
    colxn     = dbh['snpor']
    snps      = colxn.distinct('dbSNP')
    return snps

def push_nbsnpld_mongo():
    snps = get_snpor_mongo()
    connection = pymongo.Connection('ec2-50-112-118-154.us-west-2.compute.amazonaws.com',25565)
    dbh = connection["admin"]
    dbh.authenticate("adm1n","adm1npw0rd")
    dbh = connection["genetics"] # admin"]
    
    colxn     = dbh['snpld']
   
    ldfile = csv.reader(open("/Users/pouria3/Dropbox/Genophen3X/NextBio/snpLD.txt"), delimiter='\t')
    hdrs   = ldfile.next()
    for row in ldfile:
        tmp = {}
        for i in range(len(row)):
            tmp[hdrs[i]] = row[i]
            
        if tmp['snpA'] in snps and tmp['snpB'] in snps:
            print tmp
            colxn.insert(tmp)
    

if __name__ == "__main__":
    #c = Cauto(1)
    #c.main()
    #print json.dumps(c.celine, indent=4)
    
    #compare_celine("/Users/pouria3/Downloads/disease_1_selection.json", 1)
    push_nbsnpld_mongo()
    