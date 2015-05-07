import os

def rpath(x, basepath=__file__):
    return os.path.join(os.path.split(__file__)[0], '..', 'regions', x)

def specpath(x, basepath=__file__):
    return os.path.join(os.path.split(__file__)[0], '..', 'spectra', x)

def gpath(x, basepath=__file__):
    return os.path.join(os.path.split(__file__)[0], '..', 'modelgrids', x)
