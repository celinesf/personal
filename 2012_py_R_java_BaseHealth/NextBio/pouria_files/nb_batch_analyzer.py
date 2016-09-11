#!/usr/bin/env python

"""
   NextBio batch analyzer
   09/22/12 - 0.0.2 : second version after the first was unceremoniously deleted
"""
__author__ = "Pouria Mojabi"
__copyright__ = "Copyright 2012, Genophen.com"
__maintainer__ = "Pouria Mojabi"
__email__ = "mojabi@genophen.com"
__status__ = "dev" #

"""
   measure the independence and overlap of batches
"""
import collections
import pymongo
import pandas as pd
from graph_utils import nextbio_graph as nbg

def gmajors():
    majs = {
"Age-related macular degeneration" : ["Age-related macular degeneration", "Typical age-related macular degeneration"],
"Alzheimer disease" : ["Alzheimer disease", "Late-onset Alzheimer's disease"],
"Androgenic alopecia" : ["Androgenic alopecia", "Male pattern alopecia, male pattern baldness"],
"Basal cell carcinoma of skin" : ["Basal cell carcinoma of skin", "Cutaneous basal cell carcinoma"],
"Chronic viral hepatitis C" : ["Chronic viral hepatitis C", "Chronic Hepatitis C Infection" ],
"Colon cancer" : ["Colon cancer", "cancer of descending colon", "Colorectal Cancer", "caecum cancer", "Rectal cancer"],
"coronary artery disease" : ["coronary artery disease","Myocardial Infarction", "Acute Myocardial Infarction","Heart failure"],
"Diabetes mellitus type 2" : ["Diabetes mellitus type 2", "Diabetes type 2"],
"Epilepsy" : ["Epilepsy", "Partial Epilepsy"],
"Gastric carcinoma" : ["Gastric carcinoma" , "gastric cancer"],
"Glioma of brain"  : ["Glioma of brain", "Intracranial glioma"],
"Hepatitis B" : ["Hepatitis B", "Chronic Hepatitis B"],
"Hodgkin's lymphoma" : ["Hodgkin's lymphoma", "Hodgkin's lymphoma (Epstein- Barr virus negative) (Epstein- Barr virus positive)"],
"Hypertension" : ["Hypertension", "Hypertensive disease", "Hypertention"],
"Inflammatory bowel disease" : ["Inflammatory bowel disease", "IBD - Inflammatory bowel disease"],
"intracranial aneurysm" : ["intracranial aneurysm", "Intracranial arterial aneurysm"],
"Kidney Stone" : ["Kidney Stone", "Nephrolithiasis"],
"Lung cancer" : ["Lung cancer",  "Lung cancer (other subtype)"],
"Melanoma"  : ["Melanoma", "Cutaneous malignant melanoma"],
"myopia" : ["myopia", "high myopia"],
"Paget's disease" : ["Paget's disease", "Paget's disease OS"],
"Systemic sclerosis; Systemic sclerosis (ATA+); Systemic sclerosis" : ["Systemic sclerosis; Systemic sclerosis (ATA+); Systemic sclerosis", "limitied; Systemic sclerosis diffuse"],
"Testicular germ cell cancer" : ["Testicular germ cell cancer", "Testicular germ cancer", "Testicular germ cell tumor"]
}
    for k, val in majs.items():
        val2 = []
        for v in val:
            val2.append(v.lower() )
        majs[k] = val2
        
    return majs


def db2frame():
    connection = pymongo.Connection('50.112.142.17',27017)
    dbh = connection["admin"]   
    dbh.authenticate("adm1n","adm1npw0rd")
    dbh = connection["genetics"] # admin"]
    colxn     = dbh['nextbio_raw']
    curs      = colxn.find()
    
    gmajs     = gmajors() # to have one name for different naming of one disease
    
    all_batches = []
    samples = collections.defaultdict(int)
    #counter = 10
    for sample in curs:
    #    if counter == 0:
    #        break
        if "tag" in sample.keys():
            disease = sample["tag"].split("|")[0].strip()
        elif "tags" in sample.keys():
            disease = sample["tags"].split("|")[0].strip()
        else:
            print "No tags in this sample", sample
            disease = "NA"
        
        for k, val in gmajs.items():
            if disease.lower() in val:
                print disease, "was found!", " its now : ", k
                disease = k
                
        samples[disease] += 1
        batch   = sample["batch_name"]
        try:
            samplesize = [int(s) for s in sample["sample number"].replace(",", "").split() if s.isdigit()]
        except:
            samplesize = [-1, -1]

        #print "reading info for: ", disease.encode('ascii', 'ignore'), " in batch : ", batch
        
        for snp in sample['snps']:
            snp["disease"]    = disease
            snp["batch"]      = batch
            snp["samplesize"] = samplesize
            all_batches.append(snp)
        
    #    counter -=1
        
    d_frame = pd.DataFrame(all_batches)
    #print d_frame[:30][['batch', 'disease', 'snp', "snp_population", "pmid", "samplesize"]]
    
    return d_frame, samples
            


def frame2graph():
    # read everything
    df, samples = db2frame()
    
    #initialze graph
    graph = nbg()
    
    disease   = df.groupby('disease').groups.keys()
    for d in disease:
        # filter based on disease
        df_d = df[df["disease"]==d] 
        
        # find total number of snps per population for the disease
        total_snp = len(df_d.groupby('snp').groups.keys())
        all_pops  = df_d.groupby('snp_population').groups.keys()
        pop_snps  = {}
        if total_snp > 1000:
            print
            # disease is the head node
            graph.add_headnode(d)
            for pop in all_pops:
                if pop.startswith("Global"):
                    continue
                df_d_p    = df_d[df_d['snp_population']==pop]
                pop_unique_snps = len(df_d_p.groupby('snp').groups.keys())
                pop_snps[pop]   = list(df_d_p['snp'])
                pop_total_snps  = len(df_d_p['snp'])
                graph.add_triple(d, pop, d+pop,  "Tsnps: "+str(pop_total_snps), "Usnps: " + str(pop_unique_snps))
                print d.encode('ascii', 'ignore'), pop, " Tsnps: "+str(pop_total_snps), " Usnps: " + str(pop_unique_snps)
                
                # samples
                pop_unique_samples = len(df_d_p.groupby('pmid').groups.keys())
                pop_total_samples  = samples[d]
                #graph.add_triple(d, pop, d+pop, "Tsamp: " + str(pop_total_samples), "Usamp: "+str(pop_unique_samples) )
        
                # sample sizes
                ss     = df_d_p["samplesize"]
                min_ss =  min(ss)
                max_ss = max(ss)
                #try:
                #    graph.add_triple(d, pop, d+pop, "Case: " +str(min_ss[0]), "Cnt: " + str(max_ss[1]) )
                #except Exception, e:
                #    print str(e)
                #    graph.add_triple(d, pop, d+pop, "Case: " +str(min_ss[0]), "Cnt: " + str(max_ss[0]) )
                    
            
            # inter-population snp coverage
            
            pops_ = pop_snps.keys()
            for i in range(len(pop_snps.keys())-1) :
                for j in range(i+1,len(pop_snps.keys())) :
                    l1 = pop_snps[pops_[i]]
                    l2 = pop_snps[pops_[j]]
                    
                    overlap = [snp for snp in l1 if snp in l2]
                    
                    #graph.add_node_middle(d+pops_[i], d+pops_[j], "ovlp: " + str(len(overlap)) )
                    print pops_[i], pops_[j], "ovlp: " + str(len(overlap))
 
                    
                        
                    
            
        # getting batch specifi now
        #batches   = df_d.groupby('batch').groups.keys()
        #for batch in batches:
        #    ba_df = df_d[df_d["batch"] == batch]
        #    pops   = ba_df.groupby('snp_population').groups.keys() 
        #    for pop in pops:
        #        pop_df   = ba_df[ba_df["snp_population"] == pop]
        #        pop_snps = len(pop_df.groupby('snp').groups.keys())
        #        graph.add_quad(d, batch, pop, pop_snps)
            
            
            
            
    graph.write_file()
    graph.draw_graph()

    







if __name__ == "__main__":
    frame2graph()
 