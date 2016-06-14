import numpy as np
import os
import pyregion
import pyspeckit
import paths
from spectral_cube import SpectralCube
from astropy import coordinates
from astropy import units as u
import pylab as pl

regions = pyregion.open(paths.rpath('apertures.reg')) + pyregion.open(paths.rpath('gbt_pointings.reg'))
if not os.path.isdir(paths.specpath('')):
    os.mkdir(paths.specpath(''))
if not os.path.isdir(paths.specpath('figures')):
    os.mkdir(paths.specpath('figures'))

cubes = (
         ('HC3N_18-17', '/Users/adam/work/gc/hc3n/HC3N_CMZ/dropbox/50kms-HC3N_163_-50-150kms.fits'),
         ('HC3N_19-18', '/Users/adam/work/gc/hc3n/HC3N_CMZ/dropbox/50kms-HC3N_172_-50-150kms.fits'),
         ('HC3N_21-20', '/Users/adam/work/gc/hc3n/HC3N_CMZ/dropbox/50kms-HC3N_191_-50-150kms.fits'),
         ('HC3N_24-23b', '/Users/adam/work/gc/hc3n/HC3N_CMZ/dropbox/50kms-HC3N_218_-50-150kms.fits'),
         ('HC3N_24-23', '/Volumes/passport/apex/reduced/molecule_cubes/APEX_HC3N_24-23.fits'),
         ('HC3N_10-9', '/Volumes/passport/mopra/CMZ_3mm_HC3N.fits'),
         ('HC3N_3-2', '/Volumes/passport/gc_hc3n/cube_M-0.02-0.07_HC3N.fits'),
         ('HC3N_3-2', '/Volumes/passport/gc_hc3n/cube_G0.18-0.04_HC3N.fits'),
        )

for linename, cubename in cubes:
    cube = SpectralCube.read(cubename)
    for reg in regions:

        # Ignore regions that don't overlap with the cube
        try:
            subcube = cube.subcube_from_ds9region(pyregion.ShapeList([reg]))
        except ValueError:
            continue
        if any(dim == 0 for dim in subcube.shape):
            continue

        sp = subcube.mean(axis=(1,2))
        if np.all(np.isnan(sp)):
            continue

        coord = coordinates.SkyCoord(reg.coord_list[0], reg.coord_list[1],
                                     frame='galactic', unit=('deg','deg'))
        name = "G{0:06.3f}{1:+06.3f}".format(coord.l.deg, coord.b.deg)
        sp.hdu.writeto(paths.specpath("{name}_{line}.fits".format(name=name,
                                                                  line=linename)),
                       clobber=True)

        sp = pyspeckit.Spectrum.from_hdu(sp.hdu)
        sp.data -= np.median(sp.data)
        sp.plotter(figure=pl.figure(1), xmin=-1e5, xmax=1e5)
        sp.specfit(guesses=sp.specfit.moments(vheight=True, negamp=False)[1:])
        sp.baseline(excludefit=True, order=1)
        sp.specfit(guesses=sp.specfit.parinfo.values)
        sp.plotter.savefig(paths.specpath('figures/{name}_{line}.png'.format(name=name,
                                                                             line=linename)))
