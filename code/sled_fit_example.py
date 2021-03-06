import numpy as np
from sled_fitter import pyspeckit, temdencol, temdenabund
import pylab as pl

# from the APEX website, eta_mb
eta_apex = 0.75
# Using the "extended beam efficiency" from the Jones paper
eta_mopra = 0.65
# eta_gbt = eta_aperture * 1.32 to get to main beam efficiency
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
inds = np.arange(1,31)
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
constraints = mod.get_parconstraints()
pl.figure(2).clf()
mod.denstemplot()
pl.savefig('../figures/fitted_sleds/50kms_SLED_fit_lowJ_denstem/50kms_SLED_fit_lowJ_denstem.png')
pl.figure(3).clf()
mod.denscolplot()
pl.savefig('../figures/fitted_sleds/50kms_SLED_fit_lowJ_denscol/50kms_SLED_fit_lowJ_denscol.png')
pl.figure(4).clf()
mod.coltemplot()
pl.savefig('../figures/fitted_sleds/50kms_SLED_fit_lowJ_coltem/50kms_SLED_fit_lowJ_coltem.png')

from scipy import stats

cdf = stats.chi2.cdf(mod.chi2 - mod.chi2.min(), 3)
sortinds = np.argsort(cdf.ravel())
sorted_cdf = cdf.flat[sortinds]

sp.specfit.parinfo.TEMPERATURE0.value = constraints['temperature_chi2']
sp.specfit.parinfo.TEMPERATURE0.error = constraints['tmax1sig_chi2'] - constraints['temperature_chi2']
sp.specfit.parinfo.COLUMN0.value = constraints['column_chi2']
sp.specfit.parinfo.COLUMN0.error = constraints['cmax1sig_chi2'] - constraints['column_chi2']
sp.specfit.parinfo.DENSITY0.value = constraints['density_chi2']
sp.specfit.parinfo.DENSITY0.error = constraints['dmax1sig_chi2'] - constraints['density_chi2']
sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
           figure=pl.figure(5), zorder=5)

for dummy in range(100):
    val = np.random.random()
    sample_pos = np.argmin(np.abs(cdf - val))
    pars = [mod.temparr.flat[sample_pos], mod.densityarr.flat[sample_pos],
            mod.columnarr.flat[sample_pos]]
    #print "{0:0.3f}: {1:12d}, {2:7.2f}, {3:6.3f}, {4:7.3f}".format(val,
    #                                                               sample_pos,
    #                                                               *pars)
    sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                         'r-', alpha=0.1, zorder=-10)

sp.specfit.annotate()
sp.plotter.axis.set_xlabel("Rotational Level $J_U$")
sp.plotter.savefig('../figures/fitted_sleds/50kms_SLED_fit_lowJ/50kms_SLED_fit_lowJ.png')

sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
           figure=pl.figure(6), zorder=5)

temperatures = [20,50,100,300]
for t in temperatures:
    sample_pos_t = np.argmin(np.abs(mod.tarr - t))
    best_pos = np.argmin(mod.chi2[sample_pos_t,:,:].flat)
    pars = [mod.tarr[sample_pos_t],
            mod.densityarr[sample_pos_t,:,:].flat[best_pos],
            mod.columnarr[sample_pos_t,:,:].flat[best_pos]]
    label = ("T={0} K, n=$10^{{{1:0.2f}}}$ cm$^{{-3}}$, N$=10^{{{2:0.2f}}}$ "
             "cm$^{{-2}}$".format(t, pars[1], pars[2]))
    sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                         '-', alpha=0.5, linewidth=2, zorder=-10,
                         label=label)

sp.plotter.axis.set_xlabel("Rotational Level $J_U$")
sp.plotter.axis.legend(loc='upper right', fontsize=16)
sp.plotter.savefig('../figures/fitted_sleds/50kms_SLED_fit_lowJ_exampletems/50kms_SLED_fit_lowJ_exampletems.png')

pl.draw()
pl.show()
