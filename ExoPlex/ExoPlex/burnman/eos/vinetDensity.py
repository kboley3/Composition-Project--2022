# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from . import equation_of_state as eos
import warnings
from math import exp
from pdb import set_trace as b
import numpy as np


def bulk_modulus(rho, params):
    """
    compute the bulk modulus as per the third order
    Vinet equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in :math:`[Pa]`.
    """
    
    x =  params['rho_0']/rho 
    eta = (3. / 2.) * (params['Kprime_0'] - 1.)

    K = (params['K_0'] * pow(x, -2. / 3.)) * \
        (1 + ((eta * pow(x, 1. / 3.) + 1.) * (1. - pow(x, 1. / 3.)))) * \
        exp(eta * (1. - pow(x, 1. / 3.)))
    return K


def vinet(T,rho, params):
    """
    equation for the third order Vinet equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    """
    x=params['rho_0']/rho
    eta = (3. / 2.) * (params['Kprime_0'] - 1.)
    P300=3. * params['K_0'] * (x** (-2. / 3.)) * (1. - (x**( 1. / 3.))) * np.exp(eta * (1. - (x**( 1. / 3.))))
    
    gamma=mie(x,params)
    Theta=DebyeT(x,params,gamma)
    E=Energy(T, Theta)
    E300=Energy(300,Theta)
    
    Pth=(gamma/volume(rho,params))*(E-E300)
    Pressure=P300+Pth
    return Pressure


def mie(x,params):
    
        """
        Returns density :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return params['gamma0']*(x)**params['q']
    
def DebyeT(x,params,gamma):
    
        """
        Returns density :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        
        """
        
        return params['Theta0']*np.exp((params['gamma0']-gamma)/ params['q'])  
    
    
def Energy(T,theta):
    
        """
        Returns density :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        R=8.3144621 # m^3 Pa K−1 mol−1
        n=3
        x=theta/T
        E=3*(1./3.-x/8.+(x**2)/60.-(x**4)/5040.)*3*n*R*T
        
        
        return E
    
def volume(rho,params):
        """
        Returns density :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        
        return params['molar_mass'] / rho


def Density(temperature,pressure, params):
    """
    Get the Vinet density at a reference temperature for a given
    pressure :math:`[Pa]`. Returns density in :math:`[g/cm^3]`
    """

    func = lambda x: vinet(temperature,x, params) - pressure
    Rho = opt.brentq(func, 1e-6 * params['rho_0'], 1e4* params['rho_0'])
    
    #print(type(params['rho_0']))
    return Rho


class VinetD(eos.EquationOfState):
    
    """
    Base class for the isothermal Vinet equation of state.  This is third order in strain, and
    has no temperature dependence.
    """

    def density(self, pressure, temperature, params):
        """
        Returns density :math:`[kg/m^3]` as a function of pressure :math:`[Pa]`.
        """
        return Density(temperature,pressure, params)

    def pressure(self, temperature, rho, params):
        return vinet(temperature,rho , params)

    def isothermal_bulk_modulus(self, pressure, temperature, rho, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus(rho, params)

    def adiabatic_bulk_modulus(self, pressure, temperature, rho, params):
        """
        Returns adiabatic bulk modulus :math:`K_s` of the mineral. :math:`[Pa]`.
        """
        return bulk_modulus(rho, params)

    def shear_modulus(self, pressure, temperature, rho, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        Currently not included in the Vinet EOS, so omitted.
        """
        return 0.

    def heat_capacity_v(self, pressure, temperature, rho, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def heat_capacity_p(self, pressure, temperature, rho, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def thermal_expansivity(self, pressure, temperature, rho, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[1/K]`
        """
        return 0.

    def grueneisen_parameter(self, pressure, temperature, rho, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[unitless]`
        """
        return 0.
    def volume(self, pressure, temperature, params):
        """
        Returns density :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        
        return params["molar_mass"] / self.density(pressure, temperature, params)
    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        
        # G is not included in the Vinet EOS so we shall set them to NaN's
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_0' not in params:
            params['Gprime_0'] = float('nan')

        # check that all the required keys are in the dictionary
        expected_keys = ['rho_0', 'K_0', 'Kprime_0', 'gamma0','q', 'Theta0','molar_mass']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # now check that the values are reasonable.  I mostly just
        # made up these values from experience, and we are only
        # raising a warning.  Better way to do this? [IR]
        if params['rho_0'] < 1.e-2 or params['rho_0'] > 1.e5:
            warnings.warn('Unusual value for rho_0', stacklevel=2)
        if params['K_0'] < 1.e9 or params['K_0'] > 1.e13:
            warnings.warn('Unusual value for K_0', stacklevel=2)
        if params['Kprime_0'] < -5. or params['Kprime_0'] > 10.:
            warnings.warn('Unusual value for Kprime_0', stacklevel=2)
