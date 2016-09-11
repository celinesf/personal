#!/usr/bin/python
'''
Created on November 13, 2012
@author: celine

transform former json format into new
catergori: key: gender: age min max upper lower rangeee


'''
import json
import copy


################ MAIN ###############
if __name__ == "__main__":

    inputfile = open('lab_range_3.0.json')
    data = json.loads(inputfile.read())
    
    new_data = {}
    new_data['dic'] = []
    for c in data:
        print 'c',c
        for key in data[c]:
            new_doc = {}
            new_doc['c'] = c
            new_doc['v'] = key
            new_doc['u'] = data[c][key]['u']
            new_doc['d'] = data[c][key]['d']
            #new_doc['range'] = data[c][key]['d']['male']
            if 'male' in data[c][key]:
                new_doc['range'] = data[c][key]['male'][0]['d']
            
            new_data['dic'].append(new_doc)
    output =open('old_lab_ranges_out.json','w')
    output.write(json.dumps(new_data, sort_keys= True, indent=4))
    print 'DONE'

