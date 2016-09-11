#!/usr/bin/python

'''
Created on October 10, 2012
@author: celine

Module that calculates the relative position of a variant on a chromosome arm
'''
import AlgoGeneConfig as config

class RelativePosition(object):
    
    def __init__(self):
        '''
        Constructor
        '''

    ### get_relative_position ###
    ### 
    def get_relative_position(self,chromosome, position):
        info = self.get_chromosome_info(chromosome)
        centromere = int(info['centromere'])
        lenght = int(info['base_pairs'])
        if position < centromere:
            return round(round(centromere - position) /round(centromere), 1000)
        if position > centromere:
            return -round(round(position - centromere) /round(lenght - centromere), 1000)
    ### END OF get_relative_position
    
    ### get_chromosome_info ###
    ###
    def get_chromosome_info(self,chromosome):
        chr_info = None
        finddoccursor= config.GENOPHEN_DB.general.find({'name':'chromosome_size'});        
        for c in finddoccursor:
            chr_info = c['chromosome_info'][chromosome]
        return chr_info
        


