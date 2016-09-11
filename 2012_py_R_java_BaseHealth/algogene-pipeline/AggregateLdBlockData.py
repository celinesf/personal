#!/usr/bin/python

'''
Created on Jul 13, 2012
@author: celine
former LDVariantSelectionInputGenerator -> aggerate_ld_block_data
    - sql to mongo query
    - record ld_blok information data has dictionary (no more output files)
'''


import copy
import sys
import logging
import AlgoGeneConfig as config

###### AggregateLdBlockData class
class AggregateLdBlockData(object):
    ### __init__ ###
    ### Constructor
    def __init__(self, disease_name):
        logging.debug(' Function: __init__ of AggregateLdBlockData class -- disease: %s' % disease_name)
        self.version = '2.0'
        self.disease = disease_name
        self.unique_list = {}
        self.gene_block = {}
    ### END __init__


    ### getVariantGeneData ###
    # recover SNPOR data
    ###
    def getVariantGeneData(self, eth, gene_list):
        logging.debug(' Function (getVariantGeneData): ethnicity: %s' % (eth))
        finddoccursor = config.NEXTBIO_DB.nextbio.cooked.find({'disease_id':self.disease, 'broad_ethnicity':eth},
                                                              {'chromosome':1,'position_37':1,'dbsnp':1,'gene_name':1,'_id':0})
        snp_list = {}
        for doc in finddoccursor:
            if  doc['dbsnp'] not in snp_list:
                snp_list[doc['dbsnp']] = {'dbsnp' : doc['dbsnp'], 'chromosome':doc['chromosome'],'position':float(doc['position_37'])} 
            ### fine gene names
            if doc['gene_name'] is not None:
                snp_list[doc['dbsnp']]['gene'] = doc['gene_name'].split('|')
                logging.info(' found gene (getVariantGeneData) snp: %s, gene: %s' % (doc['dbsnp'], snp_list[doc['dbsnp']]['gene']))
            else:
                gene_list[doc['dbsnp']] = [doc['dbsnp']]
                snp_list[doc['dbsnp']]['gene'] = [doc['dbsnp']]
                logging.info(' NO GENE (getVariantGeneData) snp: %s, gene: %s ' % (doc['dbsnp'], snp_list[doc['dbsnp']]['gene']))
        return snp_list


    ### getEthnicityGeneList ###
    # recover list of unique ethnicity studied for this disease
    # for each ethnicity record the unique genes
    # for each genes record the unique snp and they chr and position
    # called by defineGeneBlocks
    ###
    def getEthnicityGeneList(self):
        logging.debug(' Function (getEthnicityGeneList): disease: %s' % (self.disease))
        finddoccursor = config.NEXTBIO_DB.nextbio.cooked.find({'disease_id':self.disease},
                                                              {'chromosome':1,'position_37':1,'dbsnp':1,'gene_name':1,'broad_ethnicity':1,'_id':0})
        for doc in finddoccursor:
            if doc['broad_ethnicity'] not in self.unique_list:
                self.unique_list[doc['broad_ethnicity']] = {}
            ### fine gene names
            if doc['gene_name'] is not None:
                gene_names = doc['gene_name'].split('|')
                self.unique_list[doc['broad_ethnicity']] = self.getGeneBlockFromNames(copy.deepcopy(self.unique_list[doc['broad_ethnicity']]),gene_names)    
            else:
                logging.warning( ' NO GENE (getEthnicityGeneList) %s' % doc)
            

    ### aggregateGeneBlocks ###
    # if two gene blocks are close to each other, they become one block
    # called by defineGeneBlocks
    ###
    def aggregateGeneBlocks(self,blocks, chromosome):
        logging.debug(' Function (aggregateGeneBlocks)' )
        
        for c in chromosome:
            chr_copy = copy.deepcopy(chromosome[c])
            for b1 in chromosome[c]:
                chr_copy.remove(b1)
                for b2 in chr_copy:
                    if min(abs(blocks[b2]['min']-blocks[b1]['min']), abs(blocks[b2]['max']-blocks[b1]['max']),abs(blocks[b2]['min']-blocks[b1]['max']), abs(blocks[b2]['max']-blocks[b1]['min'])) < config.LD_LIMIT:
                        logging.info('chromosome: %s, block1: %s, block2: %s (aggregateGeneBlocks) -- aggregated (distance: %d)' %( c, b1, b2,min(abs(blocks[b2]['min']-blocks[b1]['min']), abs(blocks[b2]['max']-blocks[b1]['max']),abs(blocks[b2]['min']-blocks[b1]['max']), abs(blocks[b2]['max']-blocks[b1]['min']))))
                        blockname =  '%s_%s' %(b1,b2)
                        blocks[blockname] = {}
                        blocks[blockname]['LD_block_id'] = blockname
                        blocks[blockname]['variant_list'] = blocks[b2]['variant_list']
                        blocks[blockname]['variant_list'].extend( blocks[b1]['variant_list'])
                        blocks[blockname]['chromosome'] = c
                        blocks[blockname]['min'] = min(blocks[b2]['min'], blocks[b1]['min'])
                        blocks[blockname]['max'] = max(blocks[b2]['max'], blocks[b1]['max'])
                        ## delet obsoleat block
                        del blocks[b1]
                        del blocks[b2]
                        chromosome[c].insert(chromosome[c].index(b1),blockname)
                        chromosome[c].remove(b1)
                        chromosome[c].remove(b2)
                        ## resume the loop with new lists
                        b1 = blockname
                        chr_copy = copy.deepcopy(chromosome[c])
                        chr_copy.remove(b1)
        return blocks,chromosome
    ### END OF aggregateGeneBlocks
    
    ### assignVariantToGeneBlock ###
    # lists variants within each gene block
    # called by defineGeneBlocks
    ###
    def assignVariantToGeneBlock(self,g_list, v_list):
        logging.debug(' Function (assignVariantToGeneBlock)' )
        gene_block = {}
        chromosome = {} # list of genes on each chromosome
        ### Loop on variants
        for variant in v_list:
            blockname = ""
            for gene in v_list[variant]['gene']:
                ## init block name
                if blockname == "":
                    blockname =  '_'.join(g_list[gene])
                ## check block name is complete
                else:
                    if gene not in blockname:
                        logging.critical(' variant: %s, gene: %s (assignVariantToGeneBlock) -- CANNOT FIND IN BLOCK NAME (%s)' % (variant,gene,blockname))
                        sys.exit(' variant: %s, gene: %s (assignVariantToGeneBlock) -- CANNOT FIND IN BLOCK NAME (%s)' % (variant,gene,blockname)) 
                    else:
                        logging.info(' variant: %s, gene: %s (assignVariantToGeneBlock) -- gene in block name %s' % (variant,gene,blockname))
                ## create new gene/LD block if needed
                if blockname not in gene_block:
                    gene_block[blockname] = {}
                    gene_block[blockname]['LD_block_id'] = blockname
                    gene_block[blockname]['variant_list'] = []
                    gene_block[blockname]['chromosome'] = ''
                    gene_block[blockname]['min'] = 1e10
                    gene_block[blockname]['max'] = 0
                    
                ## add variant in this gene block
                if variant not in gene_block[blockname]['variant_list']:
                    logging.info('variant: %s, gene: %s (assignVariantToGeneBlock) -- adding new variant' %( blockname, variant))
                    
                    gene_block[blockname], chromosome = self.setBlockPosition(v_list[variant],gene_block[blockname], chromosome)
                    gene_block[blockname]['variant_list'].append(variant)
                else:
                    logging.info('variant: %s, gene: %s (assignVariantToGeneBlock) -- variant already exists in this block ', blockname, variant)         
        return gene_block, chromosome
    ### END OF assignVariantToGeneBlock
    
    ### setBlockPosition ###
    # get genes on a chromosome
    # check that variant is on the right chromosome
    # find min and max positition of a gene block
    # called by assignVariantToGeneBlock
    ###
    def setBlockPosition(self,variant, block, chromosome):
        logging.debug(' Function (setBlockPosition) -- block: %s, variant: %s' % (block['LD_block_id'],variant['dbsnp']))
    
        ## set chromosome or check same chromosome
        if len(block['variant_list']) == 0 and block['chromosome'] == '' :
            block['chromosome'] = variant['chromosome']
            try:
                chromosome[variant['chromosome']].append(block['LD_block_id'])
            except:
                chromosome[variant['chromosome']] = []
                chromosome[variant['chromosome']].append(block['LD_block_id'])
        if block['chromosome'] != '' and block['chromosome'] == variant['chromosome']:
            ## find min and max position
            block['min'] = min(variant['position'], block['min'])
            block['max'] = max(variant['position'], block['max'])
        else:
            #print block['LD_block_id'], variant,str(block['chromosome']) ,str(variant['chromosome'])
            logging.CRITICAL('variant: %s, gene: %s (setBlockPosition) -- CHROMOSOME NUMBER DONT MATCH (gene: %s, variant: %s)' % (block['LD_block_id'], variant,str(block['chromosome']) ,str(variant['chromosome'])))
            sys.exit('variant: %s, gene: %s (setBlockPosition) -- CHROMOSOME NUMBER DONT MATCH (gene: %s, variant: %s)' % (block['LD_block_id'], variant,str(block['chromosome']) ,str(variant['chromosome'])))
        return block, chromosome
    ### END OF setBlockPosition

    ### getGeneBlockFromNames ###
    # recover list of unique ethnicity studied for this disease
    # called by getUniqueGenesPerEthnicity
    ###
    def getGeneBlockFromNames(self,glist, gname):
        logging.debug(' Function (getGeneBlockFromNames):  gene_name: %s' % ( gname))
        
        ### Loop on genes to add new genes in the list
        g = None # # gene counter
        for gene in gname:
            ### new gene in the list
            if gene not in glist:
                logging.info(' gene_name: %s (getGeneBlockFromNames): I am new ' % ( gene))
                if g is None: 
                    logging.info(' gene_name: %s (getGeneBlockFromNames): g none' % ( gene))
                    glist[gene] = gname
                    g = gene
                else:
                    logging.info(' gene_name1: %s (getGeneBlockFromNames): g exist %s' % ( gene,g))
                    if gene not in glist[g]:
                        glist[g].append(gene)
                    glist[gene] = glist[g]
            ### already existing gene: update block of genes
            else:   
                logging.info(' gene_name: %s (getGeneBlockFromNames): I already exist ' % ( gene))             
                if g is None: 
                    logging.info(' gene_name: %s (getGeneBlockFromNames): g none' % ( gene))
                    g = gene
                else:
                    logging.info(' gene_name2: %s (getGeneBlockFromNames): g exist %s' % ( gene,g))
                    for g2 in glist[gene] :
                        if g2 not in glist[g]:
                            glist[g].append(g2)
                    glist[gene] = glist[g] 
        return glist
    ### END OF getGeneBlockFromNames 
    
    ### defineGeneBlocks ###
    # define LD block for each ethnicity
    # Called from main
    ###
    def defineGeneBlocks(self,disease):
        logging.debug(' Function (define_gene_block): disease: %s' % (disease))
        self.getEthnicityGeneList()
        ### Loop on ethnicities
        for eth in self.unique_list:
            logging.debug(' ethnicity: %s' % (eth))
            ### get unique SNP list in this Gene
            variant_list = self.getVariantGeneData(eth, self.unique_list[eth])
            ### get variants in respective gene blocks
            self.unique_list[eth], chromosome = self.assignVariantToGeneBlock(self.unique_list[eth], variant_list)
            ### aggregate gene close to each other
            self.unique_list[eth], chromosome = self.aggregateGeneBlocks(self.unique_list[eth],chromosome)      
    
        return self.unique_list
    ### END OF defineGeneBlocks
      
    ################ MAIN ###############
    def main(self):
        logging.info(' ******************* Function: main AggregateLdBlockData -- VERSION: %s ****************' % self.version )
        self.gene_block = self.defineGeneBlocks(self.disease)
