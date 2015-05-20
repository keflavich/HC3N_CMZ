import fitter
from table_to_sled import sleds


for source,data in sleds.items():
    if len(set(data['ju'])) > 1:
        print("Fitting {0}".format(source))
        fitter.fit_a_sled(data['ju'], data['data'], data['error'], source)
    else:
        print("Skipping {0}".format(source))
