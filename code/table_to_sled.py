import numpy as np
from astropy.table import Table


tbl = Table.read('../tables/linetable.csv', format='ascii.csv')

sleds = {}
for row in tbl:
    sleds[row['Source Name']] = {'data':[], 'error':[], 'ju':[]}

for row in tbl:
    if row['Chemical Name'] == 'HC3N':
        if row['Jupper'] > 30:
            row['Jupper'] = 30
        sleds[row['Source Name']]['ju'].append(row['Jupper'])
        sleds[row['Source Name']]['data'].append(row['TMB'])
        sleds[row['Source Name']]['error'].append(row['eTMB'])

for obj in sleds:
    if 30 not in sleds[obj]['ju']:
        sleds[obj]['ju'].append(30)
        sleds[obj]['data'].append(0)
        sleds[obj]['error'].append(0.05) # highest 1-sigma error is about 0.05
    for elt in sleds[obj]:
        sleds[obj][elt] = np.array(sleds[obj][elt])
