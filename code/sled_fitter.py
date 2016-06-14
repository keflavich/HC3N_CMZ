import numpy as np
import pyradex
from pyradex import fjdu
import pyspeckit
from pyspeckit.spectrum.models.model import SpectralModel

#R = pyradex.Radex(species='hc3n-h2', abundance=1e-9, density=1e4, temperature=20)
R = pyradex.Radex(species='hc3n-h2', column=1e12, density=1e4, temperature=20)
#R = fjdu.Fjdu(species='hc3n-h2', column=1e12, density=1e4, temperature=20)
R.run_radex()

def temdencol(xarr, temperature, density, column, fortho=0, deltav=5.0):
    if density < 10:
        density = 10.**density
    if column < 30:
        column = 10.**column
    elif column < 1000:
        raise ValueError('column is absurd')
    table = R(collider_densities={'oH2':density*fortho,
                                  'pH2':density*(1-fortho)},
              temperature=temperature,
              column=column,
              deltav=5.0,
             )
    assert R.column.value == column
    np.testing.assert_almost_equal(R.total_density.value, density)
    assert R.temperature.value == temperature

    my_xarr = np.array(table['upperlevel'], dtype='int')
    result = [table['T_B'][my_xarr == ii][0] for ii in xarr]
    return np.array(result)

def temdenabund(xarr, temperature, density, abundance, fortho=0, deltav=5.0):
    table = R(collider_densities={'oH2':density*fortho,
                                  'pH2':density*(1-fortho)}, temperature=temperature,
              abundance=abundance, deltav=deltav)
    assert R.abundance == abundance
    np.testing.assert_almost_equal(R.total_density.value, density)
    assert R.temperature.value == temperature
    my_xarr = np.array(table['upperlevel'], dtype='int')
    result = [table['T_B'][my_xarr == ii][0] for ii in xarr]
    return np.array(result)

temdencolmod = SpectralModel(temdencol, npars=3,
                             parnames=['temperature','density','column'],
                             parlimits=[(2.73, 100), (0, 10), (8, 20)],
                             parlimited=[(True,True)]*3,
                             shortvarnames=['T','n(H_2)','N'],
                            )

temdenabundmod = SpectralModel(temdenabund, npars=3,
                               parnames=['temperature','density','abundance'],
                               parlimits=[(2.73, 100), (1, 1e10), (1e-16, 1e-4)],
                               parlimited=[(True,True)]*3,
                               shortvarnames=['T','n(H_2)','X'],
                              )

pyspeckit.spectrum.fitters.default_Registry.add_fitter(name='hc3n_temdencol',
                                                       function=temdencolmod,
                                                       npars=3)

pyspeckit.spectrum.fitters.default_Registry.add_fitter(name='hc3n_temdenabund',
                                                       function=temdenabundmod,
                                                       npars=3)
