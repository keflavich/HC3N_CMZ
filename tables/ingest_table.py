import numpy as np
import re
from astropy.table import Table
import string

replacedict = {r'\cyano':'HC3N',
               r'\isoa':'H13CCCN',
               r'\isob':'HC13CCN',
               r'\isoc':'HCC13CN',}

pattern = re.compile('(' + '|'.join(re.escape(x) for x in replacedict.keys()) + ')')

removelist = [r'\footnotemark[a]',
              r'\footnotemark[b]',
              r'\footnotemark{a}',
             ]
rmpattern = re.compile('(' + '|'.join(re.escape(x) for x in removelist) + ')')

with open('table.tex') as f:
    lines = f.readlines()

def floatify(lst):
    return [float(x) for x in lst]
def intify(lst):
    return [int(x) for x in lst]

suffix = ""
transition = ""
chemical_name = ""

datalines = []

datastart = False
dataend = False
for line in lines:
    if not line.strip():
        continue
    if line.split()[0] == 'Source':
        datastart = True
        continue
    if dataend or not datastart:
        continue
    if 'footnotetext' in line:
        dataend = True

    spl = [x.strip().replace("$","") for x in line.split("&")]
    if len(spl) == 9:
        if spl[0]:
            source_name = spl[0].replace(r"{\bf","").replace("}","").strip()
            suffix = ""
            transition = ""
        if spl[1]:
            if chemical_name != pattern.sub(lambda x: replacedict[x.group()], spl[1]):
                if suffix:
                    source_name = source_name[:-1]
                suffix = ""
                transition = ""
            chemical_name = pattern.sub(lambda x: replacedict[x.group()], spl[1])
        if spl[2]:
            if rmpattern.sub("", spl[2]) == transition:
                # Same transition
                if suffix == "":
                    # Not same species
                    suffix = "a"
                    source_name = source_name + suffix
                else:
                    suffix = chr(ord(suffix) + 1)
                    source_name = source_name[:-1] + suffix
            else:
                # Different transition: reset
                if suffix:
                    source_name = source_name[:-1]
                suffix = ""
            transition = rmpattern.sub("", spl[2])
            ju,jl = intify(transition.strip("()"+string.ascii_letters).split("-"))
        else:
            continue

        if spl[3]:
            aperture = int(spl[3].strip("'").replace(r'\arcsec',''))

        if "<" in spl[4]:
            peaktmb = 0
            epeaktmb = float(spl[4].split()[-1])
        elif spl[4]:
            peaktmb,epeaktmb = floatify(spl[4].split("\pm"))

        if spl[5]:
            try:
                vlsr,evlsr = floatify(spl[5].split("\pm"))
            except ValueError:
                vlsr,evlsr = np.nan,np.nan

        if spl[6]:
            width,ewidth = floatify(spl[6].split("\pm"))
        
        datalines.append([source_name, chemical_name, ju, jl, peaktmb,
                          epeaktmb, vlsr, evlsr, width, ewidth, aperture])

table = Table(data=list(zip(*datalines)),
              names=['Source Name', 'Chemical Name', 'Jupper', 'Jlower', 'TMB', 'eTMB', 'VLSR', 'eVLSR', 'width', 'ewidth', 'aperture'],
              dtype=[(str, 20), (str, 10), int, int, float, float, float, float, float, float, int])

table.write("linetable.csv", format='ascii.csv', overwrite=True)
