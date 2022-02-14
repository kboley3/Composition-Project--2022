from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from . import equation_of_state as eos
from ..tools import bracket
import warnings
from pdb import set_trace as b
import numpy as np

def bulk_modulus(rho, params):
    """
    compute the bulk modulus as per the third order
    birch-murnaghan equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in :math:`[Pa]`.
    """

    x = rho/params['rho_0']  

    f = 0.5 * (pow(x, 2. / 3.) - 1.0)

    K = pow(1. + 2. * f, 5. / 2.) * (params['K_0'] + (3. * params['K_0'] * params['Kprime_0'] -
                                                      5 * params['K_0']) * f + 27. / 2. * (params['K_0'] * params['Kprime_0'] - 4. * params['K_0']) * f * f)

    return K


def Density(pressure, params):
    func = lambda x: birch_murnaghan(x , params) - pressure
    
    try:
        sol = bracket(func,4.5* params['rho_0'],0.1 * params['rho_0'])
        if sol[0] <0 or  sol[1] <0:
            print(sol[0], sol[1])
    except:
        raise ValueError(
            'Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?')
    return opt.brentq(func, sol[0], sol[1])


def birch_murnaghan(rho, params):
    """
    equation for the Third order birch-murnaghan equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    """
    x = rho/params['rho_0']
    if x<0:
        x=-x
    
    a=3. * params['K_0'] / 2. * (pow(x, 7. / 3.) - pow(x, 5. / 3.)) \
        * (1. - .75 * (4. - params['Kprime_0']) * (pow(x, 2. / 3.) - 1.)) + params['P_0']
    if np.isnan(a):
       print(rho, x*params['rho_0'])
    
    

    return 3. * params['K_0'] / 2. * (pow(x, 7. / 3.) - pow(x, 5. / 3.)) \
        * (1. - .75 * (4. - params['Kprime_0']) * (pow(x, 2. / 3.) - 1.)) + params['P_0']


class BM3D(eos.EquationOfState):

    """
    Base class for the isothermal Birch Murnaghan equation of state.  This is fourth order in strain, and
    has no temperature dependence.
    """


    def density(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return Density(pressure, params)

    def pressure(self, temperature, volume, params):
        return birch_murnaghan(volume / params['V_0'], params)

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus(volume, params)

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus :math:`K_s` of the mineral. :math:`[Pa]`.
        """
        return bulk_modulus(volume, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        return 0.
    def volume(self, pressure, temperature, params):
        """
        Returns density :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        
        return params["molar_mass"] / self.density(pressure, temperature, params)

    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[1/K]`
        """
        return 0.

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[unitless]`
        """
        return 0.

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        if 'P_0' not in params:
            params['P_0'] = 0.

        # If G and Gprime are not included this is presumably deliberate,
        # as we can model density and bulk modulus just fine without them,
        # so just add them to the dictionary as nans
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_0' not in params:
            params['Gprime_0'] = float('nan')

        # Check that all the required keys are in the dictionary
        expected_keys = ['rho_0', 'K_0', 'Kprime_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Finally, check that the values are reasonable.
        if params['P_0'] < 0.:
            warnings.warn('Unusual value for P_0', stacklevel=2)
        if params['rho_0'] < 1.e-2 or params['rho_0'] > 1.e5:
            warnings.warn('Unusual value for rho_0', stacklevel=2)
        if params['K_0'] < 1.e9 or params['K_0'] > 1.e13:
            warnings.warn('Unusual value for K_0', stacklevel=2)
        if params['Kprime_0'] < 0. or params['Kprime_0'] > 10.:
            warnings.warn('Unusual value for Kprime_0', stacklevel=2)
       
