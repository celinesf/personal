#!/usr/bin/env python
"""
   Check that snp and/or merged snps are in the 23andme and illumina genotyping
   10/29/13 - 0.0.1 : creation
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Viragene Inc."
__maintainer__ = "Celine Becquet"
__email__ = "celine@vitagene.com"
__status__ = "dev" 


PATH = "/Users/Celine/vitagene/script_in_out/pharmGKB/check_genotyping"

import logging,os, csv
from SharedFunctions import SharedFunctions

class CheckInGenotyping():
    def __init__(self):
        self.util = SharedFunctions(PATH)
        self.dbsnp_map = self.util.get_map('dbsnp_map')
        
    def check_genotyping(self):
        logging.debug(' Function:  check_genotyping -') 
        batches_list = [i for i in os.listdir('%s/input' % (PATH))] 

        for filename in batches_list:
            self.util.warnMe('info', ' check_genotyping file %s' % filename)
            f = open('%s/input/%s' % (PATH,filename), 'rb')
            if '.csv' in filename:
                data = csv.reader(f, delimiter=',')
            elif '.txt' in filename:
                data = csv.reader(f, delimiter='\t')
            snps = [row[0] for row in data ]
            for snp in self.dbsnp_map:
                if snp not in snps:
                    if 'merged' in self.dbsnp_map[snp]:
                        found = None
                        for m in self.dbsnp_map[snp]['merged']:
                            if m in snps:
                                found = m
                                self.util.warnMe('info', ' found merged of %s - %s (check_genotyping)' % (snp, m))
                        if found is None:
                            self.util.warnMe('warning', ' %s and merged not found (check_genotyping)' % (snp))
            
                    else:
                        self.util.warnMe('warning', ' %s not found (check_genotyping)' % (snp))
                
### main function
if __name__ == "__main__":
    logging.basicConfig(filename='%s/%s' % (PATH,'log'), filemode='w',
                        level=logging.DEBUG,format='%(asctime)s - %(levelname)s -%(message)s')
    logging.debug(' Function: __main__ input' )
    d = CheckInGenotyping()
    d.check_genotyping()

    print(' DONE -- Function: __main__' )
    logging.debug(' DONE -- Function: __main__' )
""" END OF main """
