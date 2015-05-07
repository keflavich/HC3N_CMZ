import numpy as np
from sled_fitter import pyspeckit
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
           guesses=[50, 5e4, 1e12],
           fixed=[False, False, False],
           plot=False)
inds = np.arange(1,30)
sp.plotter.axis.plot(inds, sp.specfit.get_model(inds), 'o', alpha=0.5)
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
mc = sp.specfit.get_emcee()
mc.run_mcmc(mc.p0*(1+np.random.randn(*mc.p0.shape)/5.), 10)
pl.figure(2)
pl.clf()
pl.subplot(2,2,1)
pl.loglog(mc.flatchain[100:,0], mc.flatchain[100:,1], '.', alpha=0.5)
pl.subplot(2,2,2)
pl.loglog(mc.flatchain[100:,1], mc.flatchain[100:,2], '.', alpha=0.5)
pl.subplot(2,2,3)
pl.loglog(mc.flatchain[100:,0], mc.flatchain[100:,2], '.', alpha=0.5)


pl.draw()
pl.show()
