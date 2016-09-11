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

    inputfile = open('lab_ranges.json')
    data = json.loads(inputfile.read())
    
    new_data = {}

    for doc in data['dic']:
        c = doc['c'] 
        key = doc['v'] 
        gender = doc['gender'] 
        lower = doc['lower']
        upper = doc['upper']
        d = '%s - %s' % (lower,upper)
        t = 1
        high = None
        if type(upper) == list:
            t = 2
        
        if lower is None:
            d ='<= %s' % (upper)
            if len(str(upper).split('.')) == 1:
                lower = 0
            else:
                lower = 0.0                
        if upper is None:
            d ='>= %s' % (lower)
        
        if c not in new_data:
            new_data[c] = {}
        if key not in new_data[c] :
            new_data[c][key] = {}
        if gender == "both":
            new_data[c][key]['female'] = {'min_age':0, 'max_age':100, 'lower':lower, 'upper':upper, 'd':d}
            new_data[c][key]['male'] = {'min_age':0, 'max_age':100, 'lower':lower, 'upper':upper, 'd': d }
        elif gender == "male":
            new_data[c][key]['female'] = None
            new_data[c][key]['male'] = {'min_age':0, 'max_age':100, 'lower':lower, 'upper':upper, 'd':d }  
        new_data[c][key]['type'] = t
        new_data[c][key]['u'] = doc['u']
        new_data[c][key]['d'] = doc['d']
        new_data[c][key]['min'] = doc['min']
        new_data[c][key]['max'] = doc['max']
        
    output =open('lab_ranges_out.json','w')
    output.write(json.dumps(new_data, sort_keys= True, indent=4))
    print 'DONE'

