#!/usr/bin/env python

"""
Density profiles need to be set by:
1) rhos and rs
2) Integrated J-factor at a given radius and scale radius
3) Ingegrated J-factor only (need to approximate scale radius)

It would be good to generalize profile shapes, i.e., following the prescription of Zhou (1996).

Careful `rhos` is a normalization parameter and is *NOT* the same as rho(rs).

"""

import copy
from collections import OrderedDict as odict

import numpy as np
import scipy.special as spfn

from pymodeler import Model, Param
from dmsky.utils.units import Units

class DensityProfile(Model):
    _params = odict([
        ('rs',     Param(1.0)   ),
        ('rhos',   Param(1.0)   ),
        ('rmin',   Param(0.0)   ),
        ('rmax',   Param(np.inf)),
        ('rhomax', Param(np.inf)),
    ])

    def __call__(self,r):
        return self.rho(r)

    def rho(self, r):
        scalar = np.isscalar(r)
        r = np.atleast_1d(r)

        rho = self._rho(r)

        if self.rmin:   rho[r < self.rmin] = self._rho(self.rmin)
        if self.rmax:   rho[r > self.rmax] = 0
        if self.rhomax: rho[rho > self.rhomax] = self.rhomax

        if scalar:
            return np.asscalar(rho)
        else:
            return rho

    def _rho(self, r):
        msg = "%s._rho not implemented"%(self.__class__.__name__)
        raise Exception(msg)

    def mass(self,r = None):
        if r is None: r = self.rmax
        scalar = np.isscalar(r)
        r = np.atleast_1d(r)
        mass = self._mass(r) 
        if scalar: return mass[0]
        else: return mass

    def _mass(self, r):
        msg = "%s._mass not implemented"%(self.__class__.__name__)
        raise Exception(msg)
        
    def set_rho_r(self,rho,r):
        """Fix the density normalization at a given radius."""
        rhor = self._rho(r)
        self.rhos *= (rho/rhor)

    def set_mvir_c(self,mvir,c):
        rhoc = 9.9E-30*Units.g_cm3
        rvir = np.power(mvir*3.0/(177.7*4*np.pi*rhoc*0.27),1./3.)
        self.rs = rvir/c

        mrvir = self.mass(rvir)
        self.rhos *= mvir/mrvir

    def _cache(self,name=None):
        pass

class UniformProfile(DensityProfile):
    """ Uniform spherical profile
    rho(r) = rhos for r <= rs
    rho(r) = 0    otherwise
    """
    
    def _rho(self,r):
        x = r/self.rs
        return np.where(x<=1,self.rhos,0.0)

    def _mass(self,r):
        return 4*np.pi/3 * self.rhos * np.where(r < rs, r**3, self.rs**3)
   
class IsothermalProfile(DensityProfile):
    """ Non-Singular Isothermal Profile:
    Begeman et al. MNRAS 249, 523 (1991)
    http://adsabs.harvard.edu/full/1991MNRAS.249..523B
    rho(r) = rhos/(1+(r/rs))**2
    """

    def _rho(self,r):
        x = r/self.rs
        return self.rhos*(1+x)**(-2)

    def _mass(self,r):
        x = r/self.rs
        return x - (x+1)**-1 - 2*np.log(x+1)
        

class BurkertProfile(DensityProfile):
    """Burkert ApJ 447, L25 (1995) [Eqn. 2]
    http://arxiv.org/abs/astro-ph/9504041
    rho(r) = rho0 * r0**3 / ( (r + r0)*(r**2+r0**2) )
    ==>
    rho(r) = rhos / ( (1+r/rs)*(1+(r/rs)**2) )
    """
    def _rho(self,r):
        x = r/self.rs
        return self.rhos * ((1+x) * (1+x**2))**(-1)

    def mass(self,r):
        x = r/self.rs     
        return np.pi*self.rhos*(np.log(x**2+1)+2*np.log(x+1)-2*np.arctan(x))
   
class NFWProfile(DensityProfile):
    """Navarro, Frenk, and White, ApJ 462, 563 (1996)
    http://arxiv.org/abs/astro-ph/9508025
    rho(r) = rhos / ((r/rs) * (1+r/rs)**2)
    """

    #def set_jval(self,jval,rs,dist):
    #    rhos = np.sqrt(3./(4.*np.pi)*jval*dist**2/rs**3)
    #    self.rs = rs
    #    self.rhos = rhos

    def _mass(self,r):
        """ Analytic integrated mass """
        x = r/self.rs
        return 4*np.pi * self.rhos * self.rs**3 * (np.log(1+x)-x/(1+x))

    def jvalue_fast(self,r=None, dist=None):
        """Fast integrated J-factor"""
        if r is None: r = self.rmax
        x = r/self.rs
        return (4*np.pi/3.)*self.rhos**2*self.rs**3*(1-(1+x)**-3)

    def _rho(self,r):
        x = r/self.rs
        return self.rhos * x**-1 * (1+x)**-2
    
class EinastoProfile(DensityProfile):
    """ Einasto profile
    Einasto Trudy Inst. Astrofiz. Alma-Ata 5, 87 (1965) (Russian) [Eqn. 4]
    http://adsabs.harvard.edu/abs/1965TrAlm...5...87E
    rho(r) = rhos*exp(-2*((r/rs)**alpha-1)/alpha)
    ==>
    
    """
    _params = odict(
        DensityProfile._params.items() + 
        [
            ('alpha',     Param(0.17)),
        ])
    
    def _mass(self,r):
        """ Analytic mass calculation.
        FIXME: It'd be good to have a reference for this...
        """
        x = r/self.rs
        gamma = spfn.gamma(3./self.alpha)
        gammainc = spfn.gammainc(3. * self.alpha**-1,(2. * self.alpha**(-1) * x**self.alpha))
        alphainv = self.alpha**-1

        return 4*np.pi*self.rhos*self.rs**3 * alphainv * \
            np.exp(2.*alphainv)* \
            np.power(2.*alphainv, -3.*alphainv)* \
            gamma * gammainc

    def _rho(self,r):
        x = r/self.rs
        return self.rhos*np.exp(-2. * self.alpha**-1 * (x**(self.alpha) - 1) )

class GNFWProfile(DensityProfile):
    """ Generalized NFW Profile
    Strigari et al. ApJ 678, 614 (2008) [Eqn. 3]
    http://arxiv.org/abs/0709.1510
    rho(r) = rhos / ( (r/rs)**gamma * (1+r/rs)**(3-gamma))
    """
    _params = odict(
        DensityProfile._params.items() + 
        [
            ('gamma',     Param(1)),
        ])


    def _rho(self,r):
        x = r/self.rs
        return self.rhos * x**(-self.gamma) * (1+x)**(self.gamma-3)

class ZhouProfile(DensityProfile):
    """Generalized double power-law models
    Zhou MNRAS 278, 488 (1996) [Eqn. 1]
    http://arxiv.org/abs/astro-ph/9509122
    rho(r) = C * (r/rs)**-gamma * (1 + (r/rs)**1/alpha))**-(beta-gamma)*alpha
    C = 4 * rhos

    also see...
    Zhou MNRAS 287, 525 (1997) [Eqn. 2]
    http://arxiv.org/abs/astro-ph/9605029
    Strigari et al., Nature 454 (2008) [Eqn. 8]
    http://arxiv.org/abs/0808.3772
    """
    _params = odict(
        DensityProfile._params.items() + 
        [
            ('alpha',     Param(1)),
            ('beta',      Param(3)),
            ('gamma',     Param(1)),
        ])

    def _rho(self,r):
        x = r/self.rs
        return self.rhos * x**-self.gamma * (1+x**(1/self.alpha))**(-(self.beta-self.gamma)*self.alpha)

Uniform = UniformProfile
Isothermal = IsothermalProfile
Burkert = BurkertProfile
NFW = NFWProfile
Einasto = EinastoProfile
gNFW = GNFWProfile
Zhou = ZhouProfile

def factory(type, **kwargs):
    import dmsky.factory
    return dmsky.factory.factory(type, module=__name__,**kwargs)
