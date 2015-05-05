import astropy.io.fits as pyfits
from sdpy import makecube,make_off_template,calibrate_map_scans
import numpy as np
from astropy import units as u

# Set up the prefix for the input filename and the output directory here
prefix = 'AGBT06A_016_04'
rootpath = '/Volumes/passport/gc_hc3n/'
# Put in the observed zenith optical depth here
tauz = 0.03
filename = '{rootpath}{prefix}.raw.acs.fits'.format(prefix=prefix, rootpath=rootpath)
filepyfits = pyfits.open(filename,memmap=True)
datapfits = filepyfits[1].data
dataarr = datapfits.DATA

"""
These are some examples of how you'd do "introspection" on the data (i.e.,
figure out what feeds correspond to what frequencies)

In [22]: set(zip(datapfits['SAMPLER'], datapfits['FEED']))
Out[22]:
{('A1', 1),
 ('A2', 2),
 ('B3', 1),
 ('B4', 2),
 ('C5', 1),
 ('C6', 2),
 ('D7', 1),
 ('D8', 2)}

In [105]: set(zip(datapfits['SAMPLER'], datapfits['FEED'], datapfits['CRVAL1']//1e6))
Out[105]:
{('A1', 1, 26259.0), # 
 ('A1', 1, 26262.0),
 ('A2', 2, 26259.0),
 ('A2', 2, 26262.0),
 ('B3', 1, 26549.0), # NH3 8-8 26.51898
 ('B3', 1, 26552.0),
 ('B4', 2, 26549.0),
 ('B4', 2, 26552.0),
 ('C5', 1, 27029.0), # ???
 ('C5', 1, 27032.0),
 ('C6', 2, 27029.0),
 ('C6', 2, 27032.0),
 ('D7', 1, 27379.0), # NH3 9-9 27.47794e9
 ('D7', 1, 27382.0), # HC3N 27.29429e9
 ('D8', 2, 27379.0),
 ('D8', 2, 27382.0)}
"""
samplers = {
        0: ["A1","A2",],
        1: ["B3","B4",],
        2: ["C5","C6",],
        3: ["D7","D8",],
        }

feeds = {
        0: [1,2],
        1: [1,2],
        2: [1,2],
        3: [1,2]
        }
"""
scans = list(set(zip(datapfits['SCAN'], datapfits['OBJECT'])))
sorted([x for x,y in scans if y =='off1'])
sorted([x for x,y in scans if y =='off2'])
"""

# Unfortunately, none of this can be determined automatically:
# you need to give a list of the reference scans (which can be acquired
# following the approach above) and the target names and the off position names
for obsmode,refscans,scanrange,sourcename,offname in zip(
            ('RALongMap','DecLatMap','DecLatMap',),
            ([10, 13, 16, 19, 22, 25, 28, 31, 34],
             [44, 47, 50, 53, 56, 59, 62, 65, 67, 70, 73, 76, 79, 82, 85, 88, 91, 93],
             [96, 99, 102, 105, 108, 111, 114]),
            ([8,35],
             [42,92],
             [94,115],
            ),
            ('G0.18-0.04',
             'M-0.02-0.07',
             'G0.18-0.04',
            ),
            ('off2',
             'off1',
             'off2',
            ),
           ):

    # You should not need to edit anything below here

    ref1,ref2 = sorted(refscans)[0],sorted(refscans)[-1]

    for ifnum in samplers:
        for sampler,feednum in zip(samplers[ifnum],feeds[ifnum]):

            savefile = "{rootpath}/{prefix}_{0}_fd{1}_if{2}_sr{3}-{4}".format(sampler,feednum,ifnum,ref1,ref2,
                                                                  prefix=prefix,
                                                                  rootpath=rootpath)

            # This is entirely optional: leave it commented out
            # off_template = make_off_template.make_off(filename,
            #                                           scanrange=scanrange,
            #                                           sourcename=offname,
            #                                           sampler=sampler,
            #                                           feednum=feednum,
            #                                           savefile=savefile,
            #                                           clobber=True,
            #                                           extension=1,
            #                                           exclude_spectral_ends=10)

            # data = datapfits
            # meanspec = data['DATA'][(data['SAMPLER']==sampler) &
            #                              (data['OBJECT']==sourcename)].mean(axis=0)
            # meanspec = meanspec - off_template/np.median(off_template) * np.median(meanspec)
            # hdu = pyfits.PrimaryHDU(data=meanspec, header=pyfits.getheader(savefile+"_offspectra.fits"))
            # hdu.writeto(savefile+"_meanspectrum.fits", clobber=True)

            outfn = rootpath+'/'+prefix+'_%ito%i_%s_F%i.fits' % (ref1,ref2,sampler,feednum)
            calibrate_map_scans.calibrate_cube_data(filename,
                                                    outfn,
                                                    scanrange=scanrange,
                                                    #min_scale_reference=10,
                                                    feednum=feednum,
                                                    refscans=refscans,
                                                    sampler=sampler,
                                                    filepyfits=filepyfits,
                                                    datapfits=datapfits,
                                                    tauz=tauz,
                                                    dataarr=dataarr,
                                                    obsmode=obsmode,
                                                    sourcename=sourcename,
                                                    verbose=True,
                                                    #off_template=off_template,
                                                   )



