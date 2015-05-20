from astropy.table import Table


tbl = Table.read('../tables/linetable.csv')

sleds = {}
for row in tbl:
    sleds[row['Source Name']] = {'data':[], 'error':[], 'ju':[]}

for row in tbl:
    if row['Chemical Name'] == 'HC3N':
        sleds[row['Source Name']]['ju'].append(row['Jupper'])
        sleds[row['Source Name']]['data'].append(row['TMB'])
        sleds[row['Source Name']]['error'].append(row['eTMB'])
