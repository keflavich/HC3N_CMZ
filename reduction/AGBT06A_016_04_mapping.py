import os
import astropy.io.fits as pyfits
import itertools
import sys
from sdpy import makecube,make_off_template,calibrate_map_scans
import numpy as np
# to ignore div-by-zero errors?
np.seterr(all='ignore')
import FITS_tools

prefix = 'AGBT06A_016_{session:02d}'
outpath = rootpath = '/Volumes/passport/gc_hc3n/'
velocityrange = [-500,500]

sampler_lines = {'D7': 'HC3N',
                 'D8': 'HC3N',
                 'C5': 'HC3N',
                 'C6': 'HC3N',
                 'B3': 'NH388',
                 'B4': 'NH388',
                 'A1': 'NH388',
                 'A2': 'NH388',
                }
all_species = {'NH388':26.51896e9,
               'HC3N':27.29429e9}

sessions = (4,4,4,)
scans =     ([10, 34],
             [44, 93],
             [96, 114])

samplers_feeds_mapping = {
                 'D7': 'F1',
                 #'D8': 'F2',
                 'C5': 'F1',
                 #'C6': 'F2',
                 'B3': 'F1',
                 #'B4': 'F2',
                 'A1': 'F1',
                 #'A2': 'F2',
                }

names =     ('G0.18-0.04',
             'M-0.02-0.07',
             'G0.18-0.04',
            )
shapes = ([0.20298, -0.03243,40,40],
          [-0.01221,-0.07189,40,40],
          [0.20298, -0.03243,40,40])


for species in all_species:
    for name,shape in zip(names,shapes):
        cubename = os.path.join(outpath,'cube_{name}_{species}'.format(name=name, species=species))
        print name,shape,cubename
        cd3 = 2.0 # km/s
        naxis3 = (max(velocityrange)-min(velocityrange)) / cd3
        crval3 = (max(velocityrange)+min(velocityrange)) / 2.
        makecube.generate_header(shape[0], shape[1], naxis1=shape[2], naxis2=shape[3], pixsize=15,
                                 naxis3=naxis3, cd3=cd3, clobber=True, restfreq=all_species[species],
                                 crpix3=naxis3/2.,
                                 crval3=crval3, cunit3='km/s',
                                 output_cubeheader="{name}_{species}_cubeheader.txt".format(name=name, species=species),
                                 output_flatheader="{name}_{species}_header.txt".format(name=name, species=species),
                                 bmaj=28/3600., bmin=28/3600.)
        makecube.make_blank_images(cubename,clobber=True,
                                   cubeheader="{name}_{species}_cubeheader.txt".format(name=name, species=species),
                                   flatheader="{name}_{species}_header.txt".format(name=name, species=species),
                                  )

    log.info("Completed blank making {0}".format(species))

    for scan,name,session in zip(scans,names,sessions):
        cubename = os.path.join(outpath,'cube_{name}_{species}'.format(name=name, species=species))
        for sm,fd in samplers_feeds_mapping.items():
            if sampler_lines[sm] != species:
                print "Skipping {0} because it's not in this map. ({1},{2}): {3}".format(species, sm, fd, sampler_lines[sm])
                continue
            else:
                print name,cubename

            fn = os.path.join(outpath,
                              "{prefix}_{scan1}to{scan2}_{sampler}_{feed}.fits".format(scan1=scan[0],
                                                                                       scan2=scan[1],
                                                                                       prefix=prefix.format(session=session),
                                                                                       feed=fd,
                                                                                       sampler=sm)
                             )
            print "Loading file %s" % fn, " for cubename ", cubename
            data = pyfits.getdata(fn)
            fileheader = pyfits.getheader(fn)

            makecube.add_data_to_cube(
                                      cubefilename=cubename+'.fits',
                                      data=data,
                                      fileheader=fileheader,
                                      filename=fn,
                                      nhits=cubename+'_nhits.fits',
                                      cubeheader="{name}_{species}_cubeheader.txt".format(name=name, species=species),
                                      flatheader="{name}_{species}_header.txt".format(name=name, species=species),
                                      chmod=True,
                                      add_with_kernel=True,
                                      kernel_fwhm=15./3600.,
                                      velocityrange=velocityrange,
                                      diagnostic_plot_name=fn.replace('.fits','_data_scrubbed.png'),
                                      linefreq=all_species[species],
                                      smoothto=2)

    #os.system(os.path.join(outpath,'LimaBean_H2CO33_cube_starlink.sh'))


    for name in names:
        cubename = os.path.join(outpath,'cube_{name}_{species}'.format(name=name, species=species))
        sub = fits.getdata(cubename+".fits") - fits.getdata(cubename+"_continuum.fits")
        fits.PrimaryHDU(data=sub, header=fits.getheader(cubename+".fits")).writeto(cubename+"_sub.fits", clobber=True)
        makecube.make_flats(cubename,vrange=[-20,60],noisevrange=[150,200])
