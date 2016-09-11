#!/usr/bin/env python
"""
   Adjusted odds ratio calculations
   solving non linear equations
   10/15/12 - 0.0.1 : 
"""


"""
   Algorithm & Flow
   1. lay out formulas
   2. solveEquations
   3. log and potentially log the error in calculations (fit)
   
   
   and log properly
"""

from scipy.optimize import fsolve
import logging
import pandas as pd   #@UnresolvedImport
pd.set_printoptions(max_colwidth=200)
import json

class SolveEquations():
    def __init__(self, df, lifetime_prevalence):
        self.variants_data_frame = df # pandas dataframe that contains all the prameters for each variant
        self.variants_data_frame["ref_adj"] = None
        self.variants_data_frame["gt2_adj"] = None
        self.variants_data_frame["gt3_adj"] = None
        self.variants_data_frame["error"]   = None
        self.pd = lifetime_prevalence # prob of disease
        
        # function calls
        self.prepareToSolveEquation()
    
    def initializeEquations(self, initial_guess):
        (x, y, z) = initial_guess

        return (self.freqs[0]*x + self.freqs[1]*y + self.freqs[2]*z - self.pd,
                self.ors[1]*x + (1-self.ors[1])*x*y - y,
                self.ors[2]*x + (1-self.ors[2])*x*z - z
                )
    
    ''' solve equations
    return OR values  '''
    def solveEquations(self):
        ''' solveEquations the non-linear equation '''
        x, y, z = fsolve(self.initializeEquations,(1, 1, 1) )
        
        logging.info(" equation Error (solveEquations) -  x: %s \t y: %s \t z: %s" % self.initializeEquations((x, y, z)))
        return  (x/(1-x))/(self.pd/(1-self.pd)), (y/(1-y))/(self.pd/(1-self.pd)), (z/(1-z))/(self.pd/(1-self.pd)), -1 #returning odds ratios
        
    def prepareToSolveEquation(self):
        logging.debug(' Function: prepareToSolveEquation')
        for i in self.variants_data_frame.index:
            s = self.variants_data_frame.ix[i]
            self.ors     = [ s["ref_or"]  , s["gt2_or"]  , s["gt3_or"]   ]
            self.freqs   = [ s["ref_freq"], s["gt2_freq"], s["gt3_freq"] ]
            
            (x, y, z, err) = self.solveEquations()
            if x == float('inf') or y == float('inf') or z == float('inf'):
                logging.warning(' found infinity (prepareToSolveEquation) - bioset: %s, snp: %s' %(s['var'], s['bioset']))
                x = y = z = -1
            s["ref_adj"] = x
            s["gt2_adj"] = y
            s["gt3_adj"] = z
            s["error"]   = err
            self.variants_data_frame.ix[i] = s
        
        logging.info("solving finished.")
        self.variant_dict ={}
        for key in self.variants_data_frame.T.to_dict():
            key1 = str(key)
            self.variant_dict[key1] = self.variants_data_frame.T.to_dict()[key]
        logging.info(json.dumps(self.variant_dict, indent=4, sort_keys=True) )
        
        
        