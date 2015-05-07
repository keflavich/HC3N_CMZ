"""
Functions for fitting temperature (and density and column) from the line
brightnesses plus whatever other constraints are available
"""
import inspect
import time
import os

import numpy as np
from scipy.ndimage.interpolation import map_coordinates
from astropy import units as u
from astropy import log
import pylab as pl
from astropy.io import fits
from scipy import stats

from h2co_modeling import grid_fitter

def gpath(fn, gridpath='/Users/adam/work/gc/hc3n/HC3N_CMZ/modelgrids/'):
    return os.path.join(gridpath, fn)

class generic_model(object):

    def chi2_fillingfactor(self, tline, etline, lineid):
        """
        Return a chi^2 value for each model parameter treating the specified
        line brightness as a lower limit

        Parameters
        ----------
        tline : float
            The line brightness temperature
        lineid : int
            The line id
        """
        chi2 = ((self.tline[lineid] - tline)/etline)**2 * (self.tline[lineid] < tline)
        return chi2

    def chi2_column(self, logh2column, elogh2column, h2coabundance, linewidth):

        h2fromh2co = self.columnarr + np.log10(np.sqrt(np.pi) * linewidth) - h2coabundance
        chi2_h2 = ((h2fromh2co-logh2column)/elogh2column)**2

        return chi2_h2

    def chi2_abundance(self, logabundance, elogabundance):
        # Should there be a factor of np.log10(np.sqrt(np.pi) * linewidth) here?
        # Should linewidth be treated as FWHM instead of sigma, as it presently is?
        model_logabundance = self.columnarr - np.log10(u.pc.to(u.cm)) - self.densityarr
        chi2X = ((model_logabundance-logabundance)/elogabundance)**2
        return chi2X


    def get_parconstraints(self,
                           chi2level=stats.chi2.ppf(stats.norm.cdf(1)-stats.norm.cdf(-1), 3)):
        """
        If parameter constraints have been set with set_constraints or
        set_constraints_fromrow

        Parameters
        ----------
        chi2level : float
            The maximum Delta-chi^2 value to include: this will be used to
            determine the min/max 1-sigma errorbars
        """
        if not hasattr(self, 'chi2'):
            raise AttributeError("Run set_constraints first")

        row = {}

        indbest = np.argmin(self.chi2)
        deltachi2b = (self.chi2-self.chi2.min())
        for parname,pararr in zip(('temperature','column','density'),
                                  (self.temparr,self.columnarr,self.densityarr)):
            row['{0}_chi2'.format(parname)] = pararr.flat[indbest]
            OK = deltachi2b<chi2level
            if np.count_nonzero(OK) > 0:
                row['{0:1.1s}min1sig_chi2'.format(parname)] = pararr[OK].min()
                row['{0:1.1s}max1sig_chi2'.format(parname)] = pararr[OK].max()
            else:
                row['{0:1.1s}min1sig_chi2'.format(parname)] = np.nan
                row['{0:1.1s}max1sig_chi2'.format(parname)] = np.nan

        for parname in ('logh2column', 'elogh2column', 'logabundance',
                        'elogabundance'):
            if hasattr(self, parname):
                row[parname] = getattr(self, parname)

        self._parconstraints = row

        return row

    def denstemplot(self):
        self.parplot('dens','tem')

    def denscolplot(self):
        self.parplot('col','dens')

    def coltemplot(self):
        self.parplot('col','tem')

class HC3Nmodel(generic_model):

    def __init__(self, tbackground=2.73, gridsize=[250.,101.,100.],
                 lines=[(3,2),(4,3),(5,4),(10,9),(19,18),(24,23)]):
        self.lines = lines
        t0 = time.time()
        fname_template = 'fjdu_hc3n_{ju}-{jl}_{type}_5kms.fits'
        self.texgrid = {(ju,jl): fits.getdata(gpath(fname_template.format(ju=ju, jl=jl, type='tex')))
                        for ju,jl in lines}
        self.taugrid = {(ju,jl): fits.getdata(gpath(fname_template.format(ju=ju, jl=jl, type='tau')))
                        for ju,jl in lines}
        self.hdr = hdr = hdrb = fits.getheader(gpath(fname_template.format(ju=lines[0][0], jl=lines[0][1], type='tex')))

        self.Tbackground = tbackground

        self.tline_ = {trans: ((1.0-np.exp(-np.array(self.taugrid[trans]))) *
                          (self.texgrid[trans]-self.Tbackground))
                      for trans in lines}

        t1 = time.time()
        log.debug("Loading grids took {0:0.1f} seconds".format(t1-t0))

        zinds,yinds,xinds = np.indices(self.tline_[lines[0]].shape)
        upsample_factor = np.array([gridsize[0]/self.tline_[lines[0]].shape[0], # temperature
                                    gridsize[1]/self.tline_[lines[0]].shape[1], # density
                                    gridsize[2]/self.tline_[lines[0]].shape[2]], # column
                                   dtype='float')
        uzinds,uyinds,uxinds = upsinds = np.indices([x*us
                                                     for x,us in zip(self.tline_[lines[0]].shape,
                                                                     upsample_factor)],
                                                   dtype='float')
        self.tline = {trans: map_coordinates(self.tline_[trans],
                                             upsinds/upsample_factor[:,None,None,None],
                                             mode='nearest')
                      for trans in lines}

        assert self.hdr['CTYPE2'].strip() == 'LOG-DENS'
        assert self.hdr['CTYPE1'].strip() == 'LOG-COLU'

        self.columnarr = ((uxinds + self.hdr['CRPIX1']-1)*self.hdr['CDELT1'] /
                      float(upsample_factor[2])+self.hdr['CRVAL1']) # log column
        self.densityarr  = ((uyinds + self.hdr['CRPIX2']-1)*self.hdr['CDELT2'] /
                      float(upsample_factor[1])+self.hdr['CRVAL2']) # log density
        self.temparr    = ((uzinds + self.hdr['CRPIX3']-1)*self.hdr['CDELT3'] /
                      float(upsample_factor[0])+self.hdr['CRVAL3']) # lin temperature
        self.drange = [self.densityarr.min(), self.densityarr.max()]
        self.crange = [self.columnarr.min(),  self.columnarr.max()]
        self.trange = [self.temparr.min(),    self.temparr.max()]
        self.darr = self.densityarr[0,:,0]
        self.carr = self.columnarr[0,0,:]
        self.tarr = self.temparr[:,0,0]
        self.axes = {'dens': self.darr,
                     'col': self.carr,
                     'tem': self.tarr}
        self.labels = {'dens': 'Density $n(\mathrm{H}_2)$ [log cm$^{-3}$]',
                       'col': 'HC$_3$N [log cm$^{-2}$/(km s$^{-1}$ pc)]',
                       'tem': 'Temperature (K)'}

        self.model_logabundance = np.log10(10**self.columnarr / u.pc.to(u.cm) /
                                           10**self.densityarr)

        t2 = time.time()
        log.debug("Grid initialization took {0:0.1f} seconds total,"
                  " {1:0.1f} since loading grids.".format(t2-t0,t2-t1))

    def grid_getmatch_line(self, tline, etline, transition):
            match,indbest,chi2r = grid_fitter.grid_getmatch(tline, etline,
                                                            self.tline[transition])
            return chi2r

    def set_constraints(self,
                        line_brightnesses=None,
                        logabundance=None, elogabundance=None,
                        logh2column=None, elogh2column=None,
                        linewidth=None):
        """
        Set parameter constraints from a variety of inputs.  This will fill in
        a variety of .chi2_[x] values.

        All errors are 1-sigma Gaussian errors.

        Logabundance and logh2column are both log_10 values, so the errorbars
        are effectively lognormal 1-sigma errors.


        Parameters
        ----------
        line_brightnesses : dict
            A dictionary of lines like {(3,2): (5.5, 0.25)}
            where the key is the ju,jl and the value is measurement,error
        linewidth: float
            line width to be specified in km/s.
        """

        self.chi2_X = (self.chi2_abundance(logabundance, elogabundance) 
                       if not any(arg is None for arg in (logabundance,
                                                          elogabundance))
                       else 0)

        self.chi2_h2 = (self.chi2_column(logh2column, elogh2column,
                                         logabundance, linewidth) 
                        if not
                        any(arg is None for arg in (logabundance, logh2column,
                                                      elogh2column, linewidth))
                        else 0)

        self.chi2_lines = ({trans: self.grid_getmatch_line(measurement, error, trans)
                            for trans,(measurement,error) in line_brightnesses.items()}
                           if line_brightnesses is not None
                           else {})

        self.compute_chi2_fromcomponents()

    def compute_chi2_fromcomponents(self):
        """
        Compute the total chi2 from the individual chi2 components
        """
        self.chi2 = (self.chi2_X + self.chi2_h2 + np.sum(self.chi2_lines.values(), axis=0))

    def parplot(self, par1='dens', par2='col', nlevs=5, levels=None):

        xax = self.axes[par1]
        yax = self.axes[par2]
        xlabel = self.labels[par1]
        ylabel = self.labels[par2]
        axis = {('col','dens'): 0,
                ('dens','tem'): 2,
                ('col','tem'): 1}[(par1,par2)]

        if levels is None:
            levels = np.arange(nlevs)


        fig = pl.gcf()
        fig.clf()

        for index,(trans, chi2) in enumerate(self.chi2_lines.items()):
            ax = pl.subplot(2,3,index+1)
            ax.contourf(xax, yax, chi2.min(axis=axis),
                        levels=chi2.min()+levels, alpha=0.5)
            ax.contour(xax, yax, self.chi2.min(axis=axis),
                       levels=self.chi2.min()+levels)
            pl.title("J={0}-{1}".format(*trans))

        if hasattr(self, 'logabundance'):
            ax6 = pl.subplot(2,3,6)
            if hasattr(self.chi2_X, 'size'):
                ax6.contourf(xax, yax, self.chi2_X.min(axis=axis),
                             levels=self.chi2_X.min()+levels, alpha=0.5)
            ax6.contour(xax, yax, self.chi2.min(axis=axis),
                       levels=self.chi2.min()+levels)
            pl.title("log(HC$_3$N/H$_2$) "
                     "$= {0:0.1f}\pm{1:0.1f}$".format(self.logabundance,
                                                      self.elogabundance))

        fig.text(0.05, 0.5, ylabel, horizontalalignment='center',
                verticalalignment='center',
                rotation='vertical', transform=fig.transFigure)
        fig.text(0.5, 0.02, xlabel, horizontalalignment='center', transform=fig.transFigure)

        if par1 == 'col':
            for ss in range(1,5):
                ax = pl.subplot(2,3,ss)
                ax.xaxis.set_ticks(np.arange(self.carr.min(), self.carr.max()))

        pl.subplots_adjust(wspace=0.25, hspace=0.45)

    def denstemplot(self):
        self.parplot('dens','tem')

    def denscolplot(self):
        self.parplot('col','dens')

    def coltemplot(self):
        self.parplot('col','tem')

