import os
import numpy as np
from sled_fitter import pyspeckit, temdencol, temdenabund
from constrain_parameters import HC3Nmodel
from scipy import stats
import pylab as pl

def savefig(path, saver=pl):
    dir = os.path.split(path)[0]
    if not os.path.isdir(dir):
        os.mkdir(dir)
    saver.savefig(path)

def within(limits, value):
    return limits[0] < value < limits[1]

def fit_a_sled(juppers, data, error, sourcename):
    sp = pyspeckit.Spectrum(xarr=juppers, data=data, error=error, unit='$T_{MB}$ (K)',)

    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(1))
    sp.specfit(fittype='hc3n_temdencol',
               guesses=[50, np.log10(5e4), 12.],
               fixed=[False, False, False],
               plot=False, use_lmfit=True)
    inds = np.arange(1,30)
    sp.plotter.axis.plot(inds, sp.specfit.get_model(inds), 'o', alpha=0.5)
    parinfo = sp.specfit.parinfo
    # Skip this, it's just bonus viz
    #for column in np.linspace(parinfo.COLUMN0-parinfo.COLUMN0.error,
    #                          parinfo.COLUMN0+parinfo.COLUMN0.error,
    #                          5):
    #    sp.plotter.axis.plot(inds, sp.specfit.get_model_frompars(inds,
    #                                                             [parinfo.TEMPERATURE0.value,
    #                                                              parinfo.DENSITY0.value,
    #                                                              column]),
    #                         'rs', alpha=0.5)
    #sp.specfit.annotate()

    data = {(int(j),int(j)-1): (d,e) for j,d,e in zip(sp.xarr, sp.data, sp.error)}
    mod = HC3Nmodel()
    mod.set_constraints(line_brightnesses=data)
    constraints = mod.get_parconstraints()
    pl.figure(2).clf()
    mod.denstemplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_lowJ_denstem/{sourcename}_SLED_fit_lowJ_denstem.png'.format(sourcename=sourcename))
    pl.figure(3).clf()
    mod.denscolplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_lowJ_denscol/{sourcename}_SLED_fit_lowJ_denscol.png'.format(sourcename=sourcename))
    pl.figure(4).clf()
    mod.coltemplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_lowJ_coltem/{sourcename}_SLED_fit_lowJ_coltem.png'.format(sourcename=sourcename))


    cdf = stats.chi2.cdf(mod.chi2 - mod.chi2.min(), 3)
    sortinds = np.argsort(cdf.ravel())
    sorted_cdf = cdf.flat[sortinds]

    T0 = sp.specfit.parinfo.TEMPERATURE0
    if not within(T0.limits, constraints['temperature_chi2']):
        T0.limits = (2.73, 350)
    T0.value = constraints['temperature_chi2']
    T0.error = constraints['tmax1sig_chi2'] - constraints['temperature_chi2']
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
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_lowJ/{sourcename}_SLED_fit_lowJ.png'.format(sourcename=sourcename),
            sp.plotter)

    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(6), zorder=5)

    label = ("T={temperature:0.1f} K,"
             " n=$10^{{{density:0.2f}}}$ cm$^{{-3}}$,"
             " N$=10^{{{column:0.2f}}}$ "
             "cm$^{{-2}}$".format(temperature=constraints['temperature_chi2'], 
                                  density=constraints['density_chi2'],
                                  column=constraints['column_chi2'],)
            )
    sp.plotter.axis.plot(inds, temdencol(inds, 
                                         temperature=constraints['temperature_chi2'], 
                                         density=constraints['density_chi2'],
                                         column=constraints['column_chi2'],
                                        ),
                         'k-', alpha=0.5, zorder=-11,
                         label=label)

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
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_lowJ_exampletems/{sourcename}_SLED_fit_lowJ_exampletems.png'.format(sourcename=sourcename),
            sp.plotter)

    pl.draw()
    pl.show()
