import numpy as np
import os
import pyspeckit
import glob

def fit_a_spectrum(fn):
    sp = pyspeckit.Spectrum(fn)

    guesses = sp.slice(-100, 100, unit='km/s').specfit.moments()
    error_region = (sp.xarr.as_unit('km/s').value < -70) | (sp.xarr.as_unit('km/s').value > 100)
    error = sp.data[error_region].std()
    sp.error[:] = error

    sp.specfit(fittype='gaussian', guesses=guesses[1:])
    if sp.specfit.parinfo.AMPLITUDE0 < 3*error:
        print("Upper limit for {0}: {1}".format(sp.specname, error))

    else:
        sp.baseline(exclude_fitregion=True)
        sp.specfit(fittype='gaussian', guesses=sp.specfit.parinfo.values)

        sp.plotter()
        sp.specfit.plot_fit()
        sp.specfit.annotate()
        outfn = os.path.splitext(os.path.split(fn)[1])[0]
        sp.plotter.savefig('../figures/spectra/apex/{0}.png'.format(outfn))

    return sp

if __name__ == "__main__":

    for fn in glob.glob("../data/spectra/apex/*averaged*fits"):
        sp = fit_a_spectrum(fn)
