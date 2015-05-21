import numpy as np
import re
from astropy.table import Table
from astropy import coordinates
import string
import pyregion

splitrera = re.compile("\$\^\{\\\\mathrm\{[hmsd]\}\}\$")
splitredec = re.compile(r'\\(degr|arcsec|arcmin)')

with open('HC3N_posns.tex') as f:
    lines = f.readlines()

def floatify(lst):
    return [float(x) for x in lst]
def intify(lst):
    return [int(x) for x in lst]

suffix = ""
transition = ""

datalines = []
regionlines = []

datastart = False
dataend = False
for line in lines:
    if "&" not in line:
        continue
    if 'Source' in line.split("&")[0]:
        datastart = True
        print "Starting data reading at ",line
        continue
    if dataend or not datastart:
        continue
    if 'footnotetext' in line:
        dataend = True

    spl = [x.strip() for x in line.split("&")]
    if len(spl) == 4:
        if spl[0]:
            source_name = spl[0].replace(r"{\bf","").replace("}","").strip()
        else:
            continue

        ra = "{0}:{1}:{2}".format(*splitrera.split(spl[1].strip()))
        dec = "{0}:{1}:{2}".format(*[x.strip() for x in splitredec.split(spl[2].strip())[::2]])
        coord = coordinates.SkyCoord(ra+" "+dec, unit=('hour','deg'), frame='fk5')
        
        datalines.append([source_name, coord.ra.deg, coord.dec.deg])
        regionlines.append("point({ra}, {dec}) # text={{{name}}}".format(ra=coord.ra.deg, dec=coord.dec.deg, name=source_name))

table = Table(data=zip(*datalines),
              names=['Source Name', 'RA', 'Dec'],
              dtype=[(str, 20), float, float])

table.write("pointingtable.csv", format='ascii.csv')

with open("pointingtable.reg",'w') as f:
    f.write("fk5\n")
    for line in regionlines:
        f.write(line+"\n")
