from sled_fitter import pyspeckit

eta_apex = 0.75
eta_mopra = 0.65
eta_gbt = 0.68 * 1.32

# Should probably increase all of these errors by a large margin due to calibration
# uncertainty
sp = pyspeckit.Spectrum(xarr=[3, 10, 24],
                        data=[0.6097/eta_gbt,1.811/eta_mopra,0.219/eta_apex],
                        error=[0.0057/eta_gbt,0.015/eta_mopra,0.0077/eta_apex],
                        unit='$T_{MB}$ (K)',
                       )

sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30)
sp.specfit(fittype='hc3n_temdenabund',
           guesses=[50, 5e4, 1e-10],
           fixed=[False, False, True])

# this is not a well-constrained problem; probably we should throw an mcmc fitter
# at it to see the (awful) degeneracies
#sp.specfit(fittype='hc3n_temdenabund',
#           # include the model from last time, then add a second high-density,
#           # low-column component
#           guesses=[52, 4e4, 1e-10] + [100,5e4,1e-11],
#           fixed=[True,True,False] + [True, False, False])
# this doesn't seem to work.

sp.specfit.parinfo.fixed[2] = False
mc = sp.specfit.get_emcee()
mc.run_mcmc(mc.p0, 100)
