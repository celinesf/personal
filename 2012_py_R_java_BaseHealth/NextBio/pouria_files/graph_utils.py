#!/usr/bin/env python

"""
   Graph analysis
   09/22/12 - 0.0.2 : second version after the first was unceremoniously deleted
"""
__author__ = "Pouria Mojabi"
__copyright__ = "Copyright 2012, Genophen.com"
__maintainer__ = "Pouria Mojabi"
__email__ = "mojabi@genophen.com"
__status__ = "dev" #

"""
   purpose is to come up with a graph to be able to track things in nextbio and mesure overlaps of
   different batches
"""

#import matplotlib.pyplot as plt
#import networkx as nx
import os
import random


class nextbio_graph():
    def __init__(self):
        self.file = "nextbio_graph.dot"
        self.gr     = nx.Graph() # our graph structure
        
        # functions
        #self.add_nodes_edges()
        #self.write_file()
            
    
    def add_headnode(self, node):
        self.gr.add_node(node, label=node, color='green')
        
    
    def add_triple(self, head, sub, subid, subsub, subsubsub):
        nodeids = {}
        self.gr.add_node(subid, label=sub, color='red')
        nodeids[sub] = subid
        for node in [subsub, subsubsub]:
            nodeid = 0
            while self.gr.has_node(nodeid):
                nodeid = random.randint(1,10000)
            self.gr.add_node(nodeid, label=node, color='brown')
            nodeids[node] = nodeid
            
        self.gr.add_edge(head, nodeids[sub], color='blue')
        self.gr.add_edge(nodeids[sub], nodeids[subsub], color='blue')
        self.gr.add_edge(nodeids[subsub], nodeids[subsubsub], color='blue')
        

    def add_node_middle(self, nodeid1, nodeid2, newnode):
        ''' add a new node in the middle of nodes '''
        for node in [newnode]:
            nodeid = 0
            while self.gr.has_node(nodeid):
                nodeid = random.randint(1,10000)
            self.gr.add_node(nodeid, label=node, color='orange')
        
        self.gr.add_edge(nodeid1, nodeid, color='blue')
        self.gr.add_edge(nodeid, nodeid2, color='blue')

    
    def add_quad(self, head, parent, child, gchild):
        nodeids = {}
        
        #parent node is tricky:
        self.gr.add_node(head+parent, label=parent, color='pink')
        colors = {child:'pink', gchild:'pink'}
        for node in [child, gchild]:
            # find a unique nodeid
            nodeid = 0
            while self.gr.has_node(nodeid):
                nodeid = random.randint(1,10000)           
            # found a unique id, lets add nodes and edges
            self.gr.add_node(nodeid, label=node, color=colors[node])
            nodeids[node] = nodeid
        
            
        self.gr.add_edge(head             , head+parent , color='pink')
        self.gr.add_edge(head+parent      , nodeids[child]  , color='pink')
        self.gr.add_edge(nodeids[child]   , nodeids[gchild] , color='pink')
            
    
    
    
    def write_file(self):
        nx.write_dot(self.gr, self.file)
        
    
    def draw_graph(self):
        nx.draw(self.gr)  # networkx draw()
        plt.draw()        # pyplot draw()
        
    