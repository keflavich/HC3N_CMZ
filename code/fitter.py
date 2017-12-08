import os
import numpy as np
import warnings
from astropy import log
from sled_fitter import pyspeckit, temdencol, temdenabund, temdencolmod
from constrain_parameters import HC3Nmodel
from scipy import stats
import pylab as pl
import string
valid_fn_chars = "._-+"+string.ascii_letters+string.digits

mod = HC3Nmodel()
mod2 = HC3Nmodel()

colors = {50: 'r',
          100: 'b',
          150: 'm',}

def savefig(path, saver=pl):
    dir = os.path.split(path)[0]
    if not os.path.isdir(dir):
        os.mkdir(dir)
    saver.savefig(path, bbox_inches='tight')

def within(limits, value):
    return limits[0] < value < limits[1]

def fit_a_sled(juppers, data, error, sourcename, maxj=15):
    sourcename = "".join(x for x in sourcename.replace(" ","_") if x in valid_fn_chars)

    # HACK: J=31 and 32 are not allowed, use J=30 instead
    if max(juppers) > 30:
        warnings.warn("Changing some J levels from >30 to 30")
        juppers[juppers>30] = 30

    sp = pyspeckit.Spectrum(xarr=juppers, data=data, error=error, unit='$T_{MB}$ (K)',)

    # Plotting setup
    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(1, figsize=(12,8)))
    inds = np.arange(1,31)
    if len(sp.data) > 2:
        sp.specfit(fittype='hc3n_temdencol',
                   guesses=[50, np.log10(5e4), 12.],
                   fixed=[False, False, False],
                   plot=False, use_lmfit=True)
        sp.plotter.axis.plot(inds, sp.specfit.get_model(inds), 'o', alpha=0.5)

    # Fitting is done here
    linedata = {(int(j),int(j)-1): (d,e) for j,d,e in zip(juppers, data, error) if j<maxj}
    print("Low-J data: {0}".format(linedata))
    print("juppers: {0}".format(juppers))
    mod.set_constraints(line_brightnesses=linedata, logabundance=-8, elogabundance=2, logh2column=23, elogh2column=2, linewidth=20)
    assert hasattr(mod.chi2_h2, 'ndim')
    mod.chi2 += (mod.temparr < 50) * 1e10 + (mod.temparr > 200) * 1e10
    constraints = mod.get_parconstraints()
    print("Low-J parameter constraints: {0}".format(constraints))

    # Add a separate constraint for high-j components
    linedata2 = {(int(j),int(j)-1): (0,d) # upper limits
            if j<maxj
            else (d,e/5) # measurements.  increase weight on non-upper-limit
            for j,d,e in zip(juppers, data, error)}
    print("linedata2: {0}".format(linedata2))
    mod2.set_constraints(line_brightnesses=linedata2, logabundance=-8, elogabundance=2, logh2column=23, elogh2column=2, linewidth=20)
    mod2.chi2 += (mod2.temparr < 50) * 1e10 + (mod.temparr > 200) * 1e10
    constraints2 = mod2.get_parconstraints()
    print("Upper-limits for low-J plus fits from high-J: {0}".format(constraints2))

    # Plotting parameter constraints
    pl.figure(2, figsize=(12,8)).clf()
    mod.denstemplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_denstem.png'.format(sourcename=sourcename))
    pl.figure(3, figsize=(12,8)).clf()
    mod.denscolplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_denscol.png'.format(sourcename=sourcename))
    pl.figure(4, figsize=(12,8)).clf()
    mod.coltemplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_coltem.png'.format(sourcename=sourcename))

    pl.figure(2, figsize=(12,8)).clf()
    mod2.denstemplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_highJ_denstem.png'.format(sourcename=sourcename))
    pl.figure(3, figsize=(12,8)).clf()
    mod2.denscolplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_highJ_denscol.png'.format(sourcename=sourcename))
    pl.figure(4, figsize=(12,8)).clf()
    mod2.coltemplot()
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_highJ_coltem.png'.format(sourcename=sourcename))

    # Some statistics - used for overplotting a sample of models
    cdf = stats.chi2.cdf(mod.chi2 - mod.chi2.min(), 3)
    #sortinds = np.argsort(cdf.ravel())
    #sorted_cdf = cdf.flat[sortinds]
    cdf2 = stats.chi2.cdf(mod2.chi2 - mod2.chi2.min(), 3)

    # Plotting: show 100 randomly sampled model overlays selected from the
    # best-fit parameter space
    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(5, figsize=(12,8)), zorder=5)

    for dummy in range(100):
        val = np.random.random()
        sample_pos = np.argmin(np.abs(cdf - val))
        pars = [mod.temparr.flat[sample_pos],
                mod.densityarr.flat[sample_pos],
                mod.columnarr.flat[sample_pos]]
        if mod.temparr.flat[sample_pos] < 10 or mod.temparr.flat[sample_pos] > 100:
            continue
        #print("{0:0.3f}: {1:12d}, {2:7.2f}, {3:6.3f}, {4:7.3f}".format(val,
        #                                                               sample_pos,
        #                                                               *pars))
        sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                             'r-', alpha=0.1, zorder=-10)

        val = np.random.random()
        sample_pos2 = np.argmin(np.abs(cdf2 - val))
        pars2 = [mod2.temparr.flat[sample_pos2],
                 mod2.densityarr.flat[sample_pos2],
                 mod2.columnarr.flat[sample_pos2]]
        if mod.temparr.flat[sample_pos2] < 10 or mod.temparr.flat[sample_pos2] > 100:
            continue
        sp.plotter.axis.plot(inds, temdencol(inds, *pars2),
                             'b-', alpha=0.1, zorder=-10)


    # Plotting / annotation:
    if len(sp.data) > 2:
        parinfo = sp.specfit.parinfo
    else:
        sp.specfit.parinfo = parinfo = temdencolmod.make_parinfo()
        sp.specfit.fitter = temdencolmod

    T0 = parinfo.TEMPERATURE0
    if not within(T0.limits, constraints['temperature_chi2']):
        T0.limits = (2.73, 350)
    T0.value = constraints['temperature_chi2']
    T0.error = constraints['tmax1sig_chi2'] - constraints['temperature_chi2']
    parinfo.COLUMN0.value = constraints['column_chi2']
    parinfo.COLUMN0.error = constraints['cmax1sig_chi2'] - constraints['column_chi2']
    parinfo.DENSITY0.value = constraints['density_chi2']
    parinfo.DENSITY0.error = constraints['dmax1sig_chi2'] - constraints['density_chi2']
    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(5, figsize=(12,8)), zorder=5, clear=False)
    #sp.specfit.annotate()

    sp.plotter.axis.set_xlabel("Rotational Level $J_U$")
    sp.plotter.axis.set_xlim(0,32)
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit.png'.format(sourcename=sourcename),
            sp.plotter)


    # Plotting: show a selection of fixed-temperature models at best fit
    # column, volume density
    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(6, figsize=(12,8)), zorder=5)

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
                         'k--', alpha=0.8, zorder=1,
                         label=label)
    if 24 in juppers:
        label = ("T={temperature:0.1f} K,"
                 " n=$10^{{{density:0.2f}}}$ cm$^{{-3}}$,"
                 " N$=10^{{{column:0.2f}}}$ "
                 "cm$^{{-2}}$".format(temperature=constraints2['temperature_chi2'],
                                      density=constraints2['density_chi2'],
                                      column=constraints2['column_chi2'],)
                )
        sp.plotter.axis.plot(inds, temdencol(inds,
                                             temperature=constraints2['temperature_chi2'],
                                             density=constraints2['density_chi2'],
                                             column=constraints2['column_chi2'],
                                            ),
                             'k:', alpha=0.8, zorder=2,
                             label=label)

        temperatures = [50,100,150]
        for t in temperatures:
            sample_pos_t = np.argmin(np.abs(mod2.tarr - t))
            best_pos = np.argmin(mod2.chi2[sample_pos_t,:,:].flat)
            pars = [mod2.tarr[sample_pos_t],
                    mod2.densityarr[sample_pos_t,:,:].flat[best_pos],
                    mod2.columnarr[sample_pos_t,:,:].flat[best_pos]]
            label = ("T={0} K, n=$10^{{{1:0.2f}}}$ cm$^{{-3}}$, N$=10^{{{2:0.2f}}}$ "
                     "cm$^{{-2}}$".format(t, pars[1], pars[2]))
            sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                                 '-.', alpha=0.5, linewidth=1, zorder=-11,
                                 color=colors[t],
                                 label=label)

    temperatures = [50,100,150]
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
                             color=colors[t],
                             label=label)

    sp.plotter.axis.set_xlabel("Rotational Level $J_U$")
    sp.plotter.axis.legend(loc='upper right', fontsize=12)
    sp.plotter.axis.set_xlim(0,32)
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_exampletems.png'.format(sourcename=sourcename),
            sp.plotter)

    # Plotting: show a selection of fixed-density models at best fit
    # column, temperature
    sp.plotter(marker='s', linestyle='none', errstyle='bars', ymin=0, xmin=0, xmax=30,
               figure=pl.figure(6, figsize=(12,8)), zorder=5)

    label = ("T={temperature:0.1f} K,"
             " n=$10^{{{density:0.2f}}}$ cm$^{{-3}}$,"
             " N$=10^{{{column:0.2f}}}$ "
             "cm$^{{-2}}$".format(temperature=constraints['temperature_chi2'],
                                  density=constraints['density_chi2'],
                                  column=constraints['column_chi2'],)
            )
    log.info("Constraints being plotted: {0}".format(constraints))
    sp.plotter.axis.plot(inds, temdencol(inds,
                                         temperature=constraints['temperature_chi2'],
                                         density=constraints['density_chi2'],
                                         column=constraints['column_chi2'],
                                        ),
                         'k--', alpha=0.8, zorder=1,
                         label=label)
    if 24 in juppers:
        label = ("T={temperature:0.1f} K,"
                 " n=$10^{{{density:0.2f}}}$ cm$^{{-3}}$,"
                 " N$=10^{{{column:0.2f}}}$ "
                 "cm$^{{-2}}$".format(temperature=constraints2['temperature_chi2'],
                                      density=constraints2['density_chi2'],
                                      column=constraints2['column_chi2'],)
                )
        sp.plotter.axis.plot(inds, temdencol(inds,
                                             temperature=constraints2['temperature_chi2'],
                                             density=constraints2['density_chi2'],
                                             column=constraints2['column_chi2'],
                                            ),
                             'k:', alpha=0.8, zorder=2,
                             label=label)

        densities = np.log10([5e4,1e5,5e5])
        for d in densities:
            sample_pos_d = np.argmin(np.abs(mod2.darr - d))
            best_pos = np.argmin(mod2.chi2[:,sample_pos_d,:].flat)
            pars = [mod2.temparr[:,sample_pos_d,:].flat[best_pos],
                    mod2.darr[sample_pos_d],
                    mod2.columnarr[:,sample_pos_d,:].flat[best_pos]]
            label = ("T={0:0.1f} K, n=$10^{{{1:0.2f}}}$ cm$^{{-3}}$, N$=10^{{{2:0.2f}}}$ "
                     "cm$^{{-2}}$".format(pars[0], d, pars[2]))
            if pars[0] > 100:
                pars[0] = 100 # highest tem in grid
            sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                                 ':', alpha=0.5, linewidth=1, zorder=-11,
                                 label=label)

    densities = np.log10([5e3,1e4,5e4,1e5])
    for d in densities:
        sample_pos_d = np.argmin(np.abs(mod.darr - d))
        best_pos = np.argmin(mod.chi2[:,sample_pos_d,:].flat)
        pars = [mod.temparr[:,sample_pos_d,:].flat[best_pos],
                mod.darr[sample_pos_d],
                mod.columnarr[:,sample_pos_d,:].flat[best_pos]]
        if pars[0] > 100:
            pars[0] = 100 # highest tem in grid
        label = ("T={0:0.1f} K, n=$10^{{{1:0.2f}}}$ cm$^{{-3}}$, N$=10^{{{2:0.2f}}}$ "
                 "cm$^{{-2}}$".format(pars[0], d, pars[2]))
        sp.plotter.axis.plot(inds, temdencol(inds, *pars),
                             '-', alpha=0.5, linewidth=2, zorder=-10,
                             label=label)

    sp.plotter.axis.set_xlabel("Rotational Level $J_U$")
    sp.plotter.axis.legend(loc='upper right', fontsize=12)
    sp.plotter.axis.set_xlim(0,32)
    savefig('../figures/fitted_sleds/{sourcename}_SLED_fit_exampledens.png'.format(sourcename=sourcename),
            sp.plotter)


    pl.draw()
    pl.show()

    return mod, mod2, sp
