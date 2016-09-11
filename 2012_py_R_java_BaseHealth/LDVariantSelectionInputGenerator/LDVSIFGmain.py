#!/usr/bin/python

'''
Created on Jul 13, 2012
@author: celine
'''

import json
import sys
import logging
import MySQLdb #@UnresolvedImport

### log file name
#from datetime import * 
LOG_FILENAME = "logfile" #_%s" % str(datetime.now())

### get_disease_ID ###
# recover the id of the disease of interest
###
def get_disease_ID(disease_name):
    logging.debug(' Function (get_disease_ID): disease %s' % (disease_name))
    
    ### SQL INIT
    db = MySQLdb.connect(host="",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    disease_ID = 0
    sql= "select IDNo from IDmap where IDName = '%s'" % disease_name
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        if len(tab) == 1:
            disease_ID = tab[0][0]
        else:
            logging.critical(' disease_name: %s (disease_name) -- CANNOT FIND DISEASE ID (sql results: %s)' % (disease_name, tab))
            sys.exit('(disease_name) -- CANNOT FIND DISEASE ID')
    except: 
        sys.exit("Error: unable to fetch data '%s'" % sql ) 
    cursor.close ()
  
    return  disease_ID
### END OF get_disease_ID

### get_distinct_variant ###
# recover SNPOR data
###
def get_distinct_variant(disease_id,eth):
    logging.debug(' Function (get_distinct_variant): disease: %s, ethnicity: %s' % (disease_id, eth))
    
    ### SQL INIT
    db = MySQLdb.connect(host="",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    sql= " SELECT distinct dbsnp, genenextbio FROM SNPOR S where diseaseid = %d and broadethnicity ='%s'" % (disease_id,eth)
    #print sql
    snp_list = {}
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab:
            snpor = {}
            for (name, value) in zip(desc, s) :
                snpor[name[0]] = value 
            if len(snpor['dbsnp'].split('rs'))>1:
                snp_list[snpor['dbsnp']] = snpor['genenextbio'].split('_')  
    except: 
        sys.exit("Error: unable to fetch data '%s'" % sql ) 
    cursor.close ()
    return snp_list
### END OF get_distinct_variant

### define_gene_blocks ###
# define LD block for each ethnicity
# Called from main
###
def define_gene_blocks(disease):
    logging.debug(' Function (define_gene_block): disease: %s' % (disease))

    ### get unic broad ethnicity
    ethnicity_list = get_distinct_ethnicities(disease)
#    ethnicity_list = {}
#    ethnicity_list['NativeAmerican'] = {}
    ### Loop on ethnicities
    for eth in ethnicity_list:
        logging.debug(' ethnicity: %s' % (eth))
        
        ### get unique Genes list in this ethnicity
        gene_list = get_distinct_genes(disease, eth)
    
        ### get unique SNP list in this Gene
        variant_list = get_distinct_variant(disease, eth)
        
        ### get variants in respective gene blocks
        ethnicity_list[eth] = get_variant_in_gene_block(gene_list, variant_list)

    return ethnicity_list
### END OF define_gene_blocks

### get_variant_in_gene_block ###
# lists variants within each gene block
# called by define_gene_blocks
###
def get_variant_in_gene_block(g_list, v_list):
    logging.debug(' Function (get_variant_in_gene_block)' )
    gene_block = {}
    ### Loop on variants
    for variant in v_list:
        blockname = ""
        for gene in v_list[variant]:
            
            ## init block name
            if blockname == "":
                blockname =  '_'.join(g_list[gene])
            
            ## check block name is complete
            else:
                if gene not in blockname:
                    logging.critical(' variant: %s, gene: %s (get_variant_in_gene_block) -- CANNOT FIND IN BLOCK NAME (%s)' % (variant,gene,blockname))
                    sys.exit(' variant: %s, gene: %s (get_variant_in_gene_block) -- CANNOT FIND IN BLOCK NAME (%s)' % (variant,gene,blockname)) 
                else:
                    logging.critical(' variant: %s, gene: %s (get_variant_in_gene_block) -- gene in block name %s' % (variant,gene,blockname))
            
            ## create new gene/LD block if needed
            if blockname not in gene_block:
                gene_block[blockname] = {}
                gene_block[blockname]['LD_block_id'] = blockname
                gene_block[blockname]['variant_list'] = []
                
            ## add variant in this gene block
            gene_block[blockname]['variant_list'].append(variant)
            
    return gene_block
        
    
### END OF get_variant_in_gene_block

### get_distinct_ethnicities ###
# recover list of unique ethnicity studied for this disease
# called by define_gene_blocks
###
def get_distinct_ethnicities(disease_id):
    logging.debug(' Function (get_distinct_ethnicities): disease: %s' % (disease_id))
    
    ### SQL INIT
    db = MySQLdb.connect(host="",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    sql= " SELECT distinct broadethnicity FROM SNPOR S where diseaseid = %s" % (disease_id)
    #print sql
    eth_list = {}
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab: 
            for (name, value) in zip(desc, s) :
                eth_list[value] ={} 
    except: 
        sys.exit("Error: unable to fetch data '%s'" % sql ) 
    cursor.close ()
    return eth_list
### END OF get_distinct_ethnicities

### get_distinct_genes ###
# recover list of unique genes studied for this ethnicity
# called by define_gene_blocks
###
def get_distinct_genes(disease_id, eth):
    logging.debug(' Function (get_distinct_genes): disease: %s, ethnicity: %s' % (disease_id, eth))
    
    ### SQL INIT
    db = MySQLdb.connect(host="",user="celin",passwd="celin",db="AlgoPhen")
    cursor = db.cursor()
    sql= " SELECT distinct genenextbio FROM SNPOR S where diseaseid = %s and broadethnicity = '%s'" % (disease_id, eth)
    #print sql
    gene_list = {}
    try:
        cursor.execute(sql)
        tab = cursor.fetchall()
        desc = cursor.description
        for s in tab: 
            for (name, value) in zip(desc, s) : # value is genenextbio
                gene_list = get_gene_block(gene_list,value.split('_'))
    except: 
        sys.exit("Error: unable to fetch data '%s'" % sql ) 
    cursor.close ()
    return gene_list
### END OF get_distinct_genes


### get_gene_block ###
# recover list of unique ethnicity studied for this disease
# called by get_distinct_genes
###
def get_gene_block(glist, gname):
    logging.debug(' Function (get_gene_block):  gene_name: %s' % ( gname))
    
    ### Loop on genes to add new genes in the list
    g = None # # gene counter
    for gene in gname:
        
        ### new gene in the list
        if gene not in glist:
            logging.info(' gene_name: %s (get_gene_block): I am new ' % ( gene))
            if g is None: 
                logging.info(' gene_name: %s (get_gene_block): g none' % ( gene))
                glist[gene] = gname
                g = gene
            else:
                logging.info(' gene_name: %s (get_gene_block): g exist %s' % ( gene,g))
                if gene not in glist[g]:
                    glist[g].append(gene)
                glist[gene] = glist[g]
            
        ### already existing gene: update block of genes
        else:   
            logging.info(' gene_name: %s (get_gene_block): I already exist ' % ( gene))
            
            if g is None: 
                logging.info(' gene_name: %s (get_gene_block): g none' % ( gene))
                g = gene
            else:
                logging.info(' gene_name: %s (get_gene_block): g exist %s' % ( gene,g))
                for g2 in glist[gene] :
                    if g2 not in glist[g]:
                        glist[g].append(g2)
                glist[gene] = glist[g]

    return glist
### END OF get_gene_block

################ MAIN ###############
if __name__ == '__main__':

    disease_id = 1
    
    ### argument disease id
    if len(sys.argv) >1:
        disease_id = sys.argv[1]
    
    ### configure logs
    LOG_FILENAME=LOG_FILENAME+'disease%s_LD' % disease_id
    logging.basicConfig(filename=LOG_FILENAME, filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')

    ### specify disease ID
#    disease_id = get_disease_ID(disease_name)
#    print disease_name, disease_id 

    ### get * and GeneNextBio for DiseaseID order by ethnicity, chromosome position ethnicity...
    gene_block = define_gene_blocks(disease_id)

    ### output LD/gene block data
    output = open('disease_%s_LD.json' % disease_id, 'w')
    output.write(json.dumps(( gene_block),sort_keys=True, indent=4))
    output.close()

    print "I AM DONE WITH" , disease_id

    pass
