import numpy as np
import pyregion
import string
import os
import glob
from spectral_cube import SpectralCube
import FITS_tools
import pyspeckit
valid_fn_chars = "._-+"+string.ascii_letters+string.digits

regions = pyregion.open('../tables/pointingtable.reg')
radius = 20. # arcsec
for reg in regions:
    #circle = pyregion.parser_helper.Shape(shape_name='circle', shape_params=reg.coord_list+[radius/3600.])
    #circle.coord_format = reg.coord_format
    reg.coord_list = reg.coord_list + [radius/3600.]
    reg.name = 'circle'

cubefilenames = glob.glob("../data/*/*31-30*-150kms.lmv")

if not os.path.isdir('../data/spectra/'):
    os.mkdir('../data/spectra/')
if not os.path.isdir('../data/spectra/apex/'):
    os.mkdir('../data/spectra/apex/')

for cubefilename in cubefilenames:
    cube31 = SpectralCube.read(cubefilename)
    cube32filename = cubefilename.replace("31-30","32-31")
    if os.path.isfile(cube32filename):
        cube32 = SpectralCube.read(cube32filename)
    else:
        print("{0} does not exist, skipping".format(cube32filename))

    if cube31.shape != cube32.shape:
        cube32data = FITS_tools.regrid_cube_hdu(cube32.hdu, cube31.header).data
        if np.all(np.isnan(cube32data)):
            print("{0} and {1} do not overlap".format(cubefilename, cube32filename))
            continue
    else:
        cube32data = cube32.filled_data[:]

    meancubedata = (cube31.filled_data[:] + cube32data)/2.
    hdu = cube31.hdu
    hdu.data = meancubedata.value
    outcubename = cubefilename.replace("31-30", "averaged").replace(".lmv",".fits")
    hdu.writeto(outcubename, clobber=True)

    meancubedata = pyspeckit.cubes.spectral_smooth(meancubedata.value, 5, downsample=True)
    hdu.data = meancubedata
    hdu.header['CDELT3'] = hdu.header.get('CDELT3') * float(5)
    hdu.header['CRPIX3'] = hdu.header.get('CRPIX3') / float(5)
    outcubename = cubefilename.replace("31-30", "averaged_ds5").replace(".lmv",".fits")
    hdu.writeto(outcubename, clobber=True)

    cube = SpectralCube.read(hdu)

    cubename = os.path.splitext(os.path.split(outcubename)[1])[0]

    for reg in regions:
        try:
            subcube = cube.subcube_from_ds9region(pyregion.ShapeList([reg]))
        except ValueError:
            # there is probably no overlap
            continue

        name = reg.attr[1]['text']
        name = "".join(x for x in name.replace(" ","_") if x in valid_fn_chars)

        meanspec = subcube.mean(axis=(1,2))
        hdu = meanspec.hdu
        hdu.header['OBJECT'] = name
        hdu.header['ORIGFILE'] = os.path.split(cubefilename)[1]

        hdu.writeto('../data/spectra/apex/{name}_{cubename}.fits'.format(name=name,
                                                                       cubename=cubename),
                    clobber=True)
