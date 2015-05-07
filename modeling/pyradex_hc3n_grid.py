"""
Create some grids for HC3N lines

Used in the HC3N CMZ project

Grid shape is [Temperature, Density, Column Density]
"""
import pyradex
import pyradex.fjdu
import numpy as np
from astropy.utils.console import ProgressBar
from astropy.io import fits
from astropy import log
import warnings
# Make sure warnings are only shown once so the progressbar doesn't get flooded
warnings.filterwarnings('once')

ntemp,ndens,ncol = 50,20,30

transitions = [(3,2), (4,3), (5,4), (10,9), (19,18), (24,23)]

temperatures = np.linspace(10,350,ntemp)
densities = np.linspace(2.5,7,ndens)
columns = np.linspace(11, 17.1, ncol)
abundance = 1.e-9 # Johnston / Ao
opr = 0.01 # assume primarily para
fortho = opr/(1.+opr)

import os
if not os.path.exists('hc3n-h2.dat'):
    import urllib
    urllib.urlretrieve('http://home.strw.leidenuniv.nl/~moldata/datafiles/hc3n-h2.dat')

def compute_grid(densities=densities, temperatures=temperatures,
                 columns=columns, fortho=fortho, deltav=5.0,
                 escapeProbGeom='lvg', Radex=pyradex.Radex,
                 run_kwargs={'reuse_last': True, 'reload_molfile': False}):

    # Initialize the RADEX fitter with some reasonable parameters
    R = Radex(species='hc3n-h2',
              column=1e14,
              temperature=50,
              escapeProbGeom=escapeProbGeom,
              collider_densities={'oH2':2e4*fortho,'pH2':2e4*(1-fortho)})

    R.run_radex()
    R.maxiter = 200

    # get the table so we can look at the frequency grid
    table = R.get_table()

    # Get warnings about temperature early
    R.temperature = temperatures.min()
    R.temperature = temperatures.max()

    # used to assess where the grid failed
    bad_pars = []

    ntemp = len(temperatures)
    ndens = len(densities)
    ncols = len(columns)

    shape = [ntemp,ndens,ncols,]

    pars = {'taugrid_{0}-{1}'.format(ju,jl): np.full(shape, np.nan)
            for ju,jl in transitions}
    pars.update({'texgrid_{0}-{1}'.format(ju,jl): np.full(shape, np.nan)
                for ju,jl in transitions})
    pars.update({'fluxgrid_{0}-{1}'.format(ju,jl): np.full(shape, np.nan)
                 for ju,jl in transitions})

    for iTem,tt in enumerate(ProgressBar(temperatures)):
        R.temperature = tt
        for iDens,dd in enumerate(densities):
            R.density = {'oH2':10**dd*fortho,'pH2':10**dd*(1-fortho)}
            for iCol,cc in enumerate(columns):
                #R.abundance = abundance # reset column to the appropriate value
                R.column_per_bin = 10**cc
                R.deltav = deltav
                #niter = R.run_radex(reuse_last=False, reload_molfile=True)
                niter = R.run_radex(**run_kwargs)

                if niter == R.maxiter:
                    bad_pars.append([tt,dd,cc])

                TI = R.source_line_surfbrightness
                for ju,jl in transitions:
                    key = (R.upperlevelnumber == ju) & (R.lowerlevelnumber == jl)
                    pars['taugrid_{0}-{1}'.format(ju,jl)][iTem,iDens,iCol] = R.tau[key]
                    pars['texgrid_{0}-{1}'.format(ju,jl)][iTem,iDens,iCol] = R.tex[key].value
                    pars['fluxgrid_{0}-{1}'.format(ju,jl)][iTem,iDens,iCol] = TI[key].value

    return (TI, pars, bad_pars)

def makefits(data, btype, densities=densities, temperatures=temperatures,
             columns=columns, ):

    newfile = fits.PrimaryHDU(data=data)
    newfile.header.update('BTYPE' ,  btype )


    newfile.header.update('CRVAL1' ,  min(columns) )
    newfile.header.update('CRPIX1' ,  1 )
    newfile.header.update('CDELT1' , columns[1]-columns[0] )
    newfile.header.update('CTYPE1' ,  'LOG-COLU' )

    newfile.header.update('CRVAL2' ,  min(densities) )
    newfile.header.update('CRPIX2' ,  1 )
    newfile.header.update('CDELT2' , densities[1]-densities[0] )
    newfile.header.update('CTYPE2' ,  'LOG-DENS' )

    newfile.header.update('CRVAL3' ,  (min(temperatures)) )
    newfile.header.update('CRPIX3' ,  1 )
    if len(np.unique(temperatures)) == 1:
        newfile.header.update('CTYPE3' ,  'ONE-TEMP' )
        newfile.header.update('CDELT3' , temperatures[0])
    else:
        newfile.header.update('CTYPE3' ,  'LIN-TEMP' )
        newfile.header.update('CDELT3' , (np.unique(temperatures)[1]) - (np.unique(temperatures)[0]) )
    return newfile

if __name__ == "__main__":
    import re
    from paths import gpath
    import os
    if not os.path.isdir(gpath('')):
        os.mkdir(gpath(''))
    bt = re.compile("tex|tau|flux")

    (fTI, fpars, fbad_pars) = compute_grid(Radex=pyradex.fjdu.Fjdu,
                                           run_kwargs={})
    
    for pn in fpars:
        btype = bt.search(pn).group()
        ff = makefits(fpars[pn], btype, densities=densities,
                      temperatures=temperatures, columns=columns)
        outfile = 'fjdu_hc3n_{line}_{type}_{dv}.fits'.format(line=pn.split("_")[-1],
                                                              type=btype,
                                                              dv='5kms')
        ff.writeto(gpath(outfile),
                   clobber=True)
        print outfile

    #ff = makefits(fpars['fluxgrid_321']/fpars['fluxgrid_303'], 'ratio',
    #              densities=densities, temperatures=temperatures,
    #              columns=columns)
    #outfile = 'fjdu_hc3n_{line}_{type}_{dv}.fits'.format(line='321to303',
    #                                                      type='ratio',
    #                                                      dv='5kms')
    #ff.writeto(gpath(outfile), clobber=True)

    (TI, pars, bad_pars) = compute_grid()

    for pn in pars:
        btype = bt.search(pn).group()
        ff = makefits(pars[pn], btype, densities=densities,
                      temperatures=temperatures, columns=columns)
        outfile = 'hc3n_{line}_{type}_{dv}.fits'.format(line=pn.split("_")[-1],
                                                          type=btype,
                                                          dv='5kms')
        ff.writeto(gpath(outfile),
                   clobber=True)
        print outfile

    #ff = makefits(pars['fluxgrid_321']/pars['fluxgrid_303'], 'ratio',
    #         densities=densities, temperatures=temperatures, columns=columns)
    #outfile = 'hc3n_{line}_{type}_{dv}.fits'.format(line='321to303',
    #                                                      type='ratio',
    #                                                      dv='5kms')
    #ff.writeto(gpath(outfile), clobber=True)

    log.info("FJDU had {0} bad pars".format(len(fbad_pars)))
    log.info("RADEX had {0} bad pars".format(len(bad_pars)))
    

    # look at differences
    for pn in pars:
        btype = bt.search(pn).group()
        outfile = 'hc3n_{line}_{type}_{dv}.fits'.format(line=pn.split("_")[-1],
                                                          type=btype,
                                                          dv='5kms')
        header = fits.getheader(gpath(outfile))
        im1 = fits.getdata(gpath('fjdu_'+outfile))
        im2 = fits.getdata(gpath(outfile))
        hdu = fits.PrimaryHDU(data=im1-im2, header=header)
        hdu.writeto(gpath('diff_fjdu-radex_'+outfile), clobber=True)

