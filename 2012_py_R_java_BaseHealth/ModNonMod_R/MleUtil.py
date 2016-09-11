#!/usr/bin/env python

"""
    Function to get MLE format citations
    06/11/13: version 1.0
"""
__author__ = "Celine Becquet"
__copyright__ = "Copyright 2013, Genophen.com"
__maintainer__ = "Celine Becquet"
__email__ = "becquet@genophen.com"
__status__ = "dev" 
"""
    1) get citation data functions called in GetModNonModRiskFactors : getAuthors, getContent
    2) utility functions for citation: getInfo, formatName
"""

import logging
import datetime
from InputOutputUtil import InputOutputUtil

DATE = datetime.datetime.now().strftime("%m-%d-%Y")

class MleUtil():
    def __init__(self,path, disease_name):
        self.disease_name = disease_name
        self.path = path
        self.util = InputOutputUtil(self.path, self.disease_name)
    
##################### utility functions for citation ###################
    '''  getInfo '''
    def getInfo(self, reference, key,content, start_sep, end_sep):
        logging.debug(' Function:  getInfo- pmid %s, key:%s' % (reference['pmid'],key))
        if key in reference:
            info = "%s%s%s%s" % (content, start_sep, reference[key],end_sep)
        else:
            info = content
        return info
     
    ''' formatName '''
    def formatName(self, author):
        logging.debug(' Function:  formatName- author %s' % author)   
        author_name = author.split(', ')
        if len(author_name) >1:
            last = author_name[0]
            if len(author_name)>2:
                self.util.warnMe('warning', ' TOO MANY FIRST NAMES (formatName) - %s' % (author_name))
            middle = author_name[1].split(' ')
            nname = 0
            for name in middle:
                if nname == 0 and len(name) >1:
                    first = name
                    nname += 1
                elif nname > 0 and len(name) >1:
                    first = "%s %s" %(first, name)
                elif nname > 0 and len(name)  == 1:
                    first = "%s %s." %(first, name)
                elif nname == 0 and len(name)  == 1:
                    first = "%s." %( name)
                else:
                    self.util.warnMe('warning', ' WEIRD FIRST NAME (formatName) - %s' % (author_name))
            return  "%s, %s" % (last, first)   
        else:
            return "%s" % (author_name[0])  
#####################


##################### get citation data functions called in GetModNonModRiskFactors ###################
    '''  getAuthors '''
    def getAuthors(self, reference):
        logging.debug(' Function:  getAuthors- pmid %s' % reference['pmid'])
        if reference['author'] is not None:
            authors = reference['author'].split(' and ')
            if len(authors) > 2:# first author last, first init., et al.
                first_author = self.formatName(authors[0])
                author_info = "%s, et al." % first_author
            elif len(authors) == 2:
                first_author = self.formatName(authors[0])
                second_author = self.formatName(authors[1])
                author_info = "%s, and %s." % ( first_author,second_author)
            elif len(authors) == 1:
                author_info = authors[0]
            return author_info
        else: return None

    '''  getContent '''
    def getContent(self, reference):
        logging.debug(' Function:  getContent- pmid %s' % reference['pmid'])
        content = " \'%s.\' " % reference['title']
        content = self.getInfo(reference,'journal', content, "", " ") 
        content = self.getInfo(reference,'volume', content, "", "") 
        content = self.getInfo(reference,'number', content, ".", "") 
        content = self.getInfo(reference,'year', content, " (", ")") 
        content = self.getInfo(reference,'pages', content, ": ", "") 
        
#         date = reference['date-added'].split(' ')
#         date = "%s %s" %(date[0],date[1])
#         date = datetime.datetime.strptime(date , "%Y-%m-%d %H:%M:%S").strftime("%d %b. %Y.")
#         content = "%s. Web. %s" % (content, date)        
        return content
#####################
