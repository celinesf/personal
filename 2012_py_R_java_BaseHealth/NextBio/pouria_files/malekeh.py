import json
import pandas as pd

json_data = open('/Users/pouria3/Downloads/GRE_PA_Rules-2.json')
data      = json.load(json_data)
data_df   = pd.DataFrame(data)

risks     = selected_variant_df.groupby('risk_factor').groups.keys()
for risk in risks:
    df = data_df[data_df['risk_factor'] == risk]
    #print ">"*40 + "\t" + risk + "\t" + "<"*40
    #print df[df['met']=="0-10"][['g_major', 'g_minor', 'duration', 'duration_units', 'frequency', 'frequency_units']]
    gmaj = df.groupby('g_major').groups.keys()
    for i in range(len(gmaj)):
        tmp = df[df['g_major'] == gmaj[i]]
        print tmp[['g_major', 'g_minor', 'duration', 'duration_units', 'frequency', 'frequency_units', 'met_inc', 'risk_factor']]
    
    
