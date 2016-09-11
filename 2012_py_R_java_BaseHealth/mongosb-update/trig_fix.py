#!/usr/bin/python
'''
Created on November 13, 2012
@author: celine

fix the triglyceride levels in the database for t2d and stroke

'''
import json
import copy



################ MAIN ###############
if __name__ == "__main__":

    filename = 'stk'
    inputfile = open('%s.json' % filename)
    data = json.loads(inputfile.read())
    newdata = copy.deepcopy(data)
    trig = []
    for doc  in data['risk_factors']:
        if doc['attr_name'] == 'triglyceride_level':
            if doc['attr_value'] != -999:
                trig.append({'attr_value':doc['attr_value']*10,'attr_name':'triglyceride_level','factor_value':doc['factor_value']})
                print doc
                newdata['risk_factors'].remove(doc)


    newtrig = []
    index = 0
    for i in range(0,len(trig)):
        newtrig.append(trig[i])
        if i+1 < len(trig):
            val = trig[i]['attr_value']
            OR = trig[i]['factor_value']
            for j in range(1,5):
                newtrig.append({'attr_value':val + j,'attr_name':'triglyceride_level','factor_value':OR})
            OR = trig[i+1]['factor_value']
            for j in range(5,10):
                newtrig.append({'attr_value':val + j,'attr_name':'triglyceride_level','factor_value':OR})
  
    
    newdata['risk_factors'].extend(newtrig)
    output =open('%s_new.json' %filename,'w')
    output.write(json.dumps(newdata, sort_keys= True, indent=4))
    print 'DONE'

