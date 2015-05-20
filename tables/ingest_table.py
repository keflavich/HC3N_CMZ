import numpy as np
import re
from astropy.table import Table
import string

replacedict = {r'\cyano':'HC3N',
               r'\isoa':'H13CCCN',
               r'\isob':'HC13CCN',
               r'\isoc':'HCC13CN',}

pattern = re.compile('(' + '|'.join(re.escape(x) for x in replacedict.keys()) + ')')

removelist = r'\footnotemark[a]'
rmpattern = re.compile('(' + '|'.join(re.escape(x) for x in removelist) + ')')

with open('table.tex') as f:
    lines = f.readlines()

def floatify(lst):
    return [float(x) for x in lst]
def intify(lst):
    return [int(x) for x in lst]

datalines = []

datastart = False
dataend = False
for line in lines:
    if line.split()[0] == 'Source':
        datastart = True
        continue
    if dataend or not datastart:
        continue
    if 'footnotetext' in line:
        dataend = True

    spl = [x.strip().replace("$","") for x in line.split("&")]
    if len(spl) == 8:
        if spl[0]:
            source_name = spl[0].replace(r"{\bf","").replace("}","")
        if spl[1]: 
            chemical_name = pattern.sub(lambda x: replacedict[x.group()], spl[1])
        if spl[2]:
            transition = rmpattern.sub("", spl[2])
            ju,jl = intify(transition.strip("()"+string.letters).split("-"))
        else:
            continue

        if "<" in spl[3]:
            peaktmb = np.nan
            epeaktmb = float(spl[3].split()[-1])
        elif spl[3]:
            peaktmb,epeaktmb = floatify(spl[3].split("\pm"))

        if spl[4]:
            try:
                vlsr,evlsr = floatify(spl[4].split("\pm"))
            except ValueError:
                vlsr,evlsr = np.nan,np.nan

        if spl[5]:
            width,ewidth = floatify(spl[5].split("\pm"))
        
        datalines.append([source_name, chemical_name, ju, jl, peaktmb,
                          epeaktmb, vlsr, evlsr, width, ewidth])

table = Table(data=zip(*datalines),
              names=['Source Name', 'Chemical Name', 'Jupper', 'Jlower', 'TMB', 'eTMB', 'VLSR', 'eVLSR', 'width', 'ewidth'],
              dtype=[(str, 20), (str, 10), int, int, float, float, float, float, float, float])

table.write("linetable.csv", format='ascii.csv')
