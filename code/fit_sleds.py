import numpy as np
import fitter
from table_to_sled import sleds


for source,data in sleds.items():
    if source not in ('M-0.02-0.07', 'M-0.13-0.08 (b1)', 'M+0.25+0.01'):
        continue

    # Filter out the high-J lines for SLED fitting (we'll deal with them separately)
    #high = data['ju'] > 15

    #data['data'] = data['data'][~high]
    #data['error'] = data['error'][~high]
    #data['ju'] = data['ju'][~high]

    ok = data['data'] >= 0
    uplim = data['data'] == 0

    if len(set(data['ju'][ok])) > 1:
        print("Fitting {0}".format(source))
        # uses the measurement error: fitter.fit_a_sled(data['ju'], data['data'], data['error'], source)

        errors = np.array(data['data'][ok])*0.2
        errors[uplim[ok]] = data['error'][uplim & ok] * 5 # conservative: 5-sigma upper limits

        # make sure we don't pass in any negative data, because that's silly.
        mod1, mod2, sp = fitter.fit_a_sled(data['ju'][ok], data['data'][ok], errors, source)
    else:
        print("Skipping {0}".format(source))
