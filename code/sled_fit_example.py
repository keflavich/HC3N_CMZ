import numpy as np
from sled_fitter import pyspeckit, temdencol, temdenabund
import pylab as pl

eta_apex = 0.75
eta_mopra = 0.65
#eta_gbt_5 = 0.68 * 1.32
#eta_gbt_4 = 0.70 * 1.32
eta_gbt_3 = 0.68 * 1.32

# Should probably increase all of these errors by a large margin due to calibration
# uncertainty
sp = pyspeckit.Spectrum(xarr=[3, 4, 5, 10,],# 24],
                        data=[3.42307/eta_gbt_3, 5.67, 6.6, 1.811/eta_mopra],#,0.219/eta_apex],
                        error=[0.0278719/eta_gbt_3*10, 0.32, 0.4, 0.015/eta_mopra*10],#,0.0077/eta_apex],
                        unit='$T_{MB}$ (K)',
                       )

sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
           figure=pl.figure(1))
sp.specfit(fittype='hc3n_temdencol',
           guesses=[50, np.log10(5e4), 12.],
           fixed=[False, False, False],
           plot=False, use_lmfit=True)
inds = np.arange(1,30)
sp.plotter.axis.plot(inds, sp.specfit.get_model(inds), 'o', alpha=0.5)
parinfo = sp.specfit.parinfo
for column in np.linspace(parinfo.COLUMN0-parinfo.COLUMN0.error,
                          parinfo.COLUMN0+parinfo.COLUMN0.error,
                          5):
    sp.plotter.axis.plot(inds, sp.specfit.get_model_frompars(inds,
                                                             [parinfo.TEMPERATURE0.value,
                                                              parinfo.DENSITY0.value,
                                                              column]),
                         'rs', alpha=0.5)
sp.specfit.annotate()

# this is not a well-constrained problem; probably we should throw an mcmc fitter
# at it to see the (awful) degeneracies
#sp.specfit(fittype='hc3n_temdenabund',
#           # include the model from last time, then add a second high-density,
#           # low-column component
#           guesses=[52, 4e4, 1e-10] + [100,5e4,1e-11],
#           fixed=[True,True,False] + [True, False, False])
# this doesn't seem to work.

#sp.specfit.parinfo.fixed[2] = False
#mc = sp.specfit.get_emcee()
#mc.run_mcmc(mc.p0*(1+np.random.randn(*mc.p0.shape)/5.), 10)
#pl.figure(2)
#pl.clf()
#pl.subplot(2,2,1)
#pl.loglog(mc.flatchain[100:,0], mc.flatchain[100:,1], '.', alpha=0.5)
#pl.subplot(2,2,2)
#pl.loglog(mc.flatchain[100:,1], mc.flatchain[100:,2], '.', alpha=0.5)
#pl.subplot(2,2,3)
#pl.loglog(mc.flatchain[100:,0], mc.flatchain[100:,2], '.', alpha=0.5)

from constrain_parameters import HC3Nmodel

data = {(int(j),int(j)-1): (d,e) for j,d,e in zip(sp.xarr, sp.data, sp.error)}
mod = HC3Nmodel()
mod.set_constraints(line_brightnesses=data)
mod.get_parconstraints()
pl.figure(2).clf()
mod.denstemplot()
pl.figure(3).clf()
mod.denscolplot()
pl.figure(4).clf()
mod.coltemplot()

from scipy import stats

cdf = stats.chi2.cdf(mod.chi2 - mod.chi2.min(), 3)
sortinds = np.argsort(cdf.ravel())
sorted_cdf = cdf.flat[sortinds]

sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
           figure=pl.figure(5))

for dummy in range(100):
    val = np.random.random()
    sample_pos = np.argmin(np.abs(cdf - val))
    pars = [mod.temparr.flat[sample_pos], mod.densityarr.flat[sample_pos],
            mod.columnarr.flat[sample_pos]]
    #print "{0:0.3f}: {1:12d}, {2:7.2f}, {3:6.3f}, {4:7.3f}".format(val,
    #                                                               sample_pos,
    #                                                               *pars)
    sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                         'r-', alpha=0.2)


pl.draw()
pl.show()
