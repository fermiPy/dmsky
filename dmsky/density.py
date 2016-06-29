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

from pymodeler.model import Model
from pymodeler.parameter import *

from dmsky.utils.units import Units

class DensityProfile(Model):
    _params = odict([
        ('rs',     Parameter(default=1.0)   ),
        ('rhos',   Parameter(default=1.0)   ),
        ('rmin',   Parameter(default=0.0)   ),
        ('rmax',   Parameter(default=np.inf)),
        ('rhomax', Parameter(default=np.inf)),
        ('covar',  Derived(dtype=np.ndarray     ,comment='Covariance matrix for parameters')),
    ])

    def __call__(self,r):
        """ return the denisty at a given radius
        """
        return self.rho(r)

    @property 
    def deriv_params(self):
        return ["rs","rhos"]


    def rho(self, r):
        """ return the denisty at a given radius
        """
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

    def rho_deriv(self, r, paramNames):
        """ return the deirvatives of the density as a function of radius, 
            w.r.t. a list of parameters

            This will return an n x m array, where
            n is the number of radii
            m is the number of parameters              
        """
        if np.isscalar(r):
            nr = 1
        else:
            nr = len(r)
        
        npar = len(paramNames)
        
        initParVals = np.array([ self.__getattr__(pName) for pName in paramNames])
        deltaParVals = initParVals * 0.001
        
        init_r = self._rho(r)

        derivs = []

        # loop over parameters and take the numerical derivatives
        for initPar,deltaPar,parName in zip(initParVals,deltaParVals,paramNames):
            par = self.getp(parName)
            newParVal = initPar+deltaPar
            par.set_value(newParVal)
            new_r = self._rho(r)
            dr_dp = (new_r - init_r) / ( newParVal -  initPar)
            derivs.append(dr_dp)
            par.set_value(initPar)

        ret = np.vstack(derivs)
        return ret

    def rho_uncertainty(self, r):
        """
        """
        cov_mat = self.covar
        if np.isscalar(r):
            nr = 1
            deriv_vect = np.matrix(self.rho_deriv(r,self.deriv_params))
            err2 = (deriv_vect.T * cov_mat * deriv_vect)[0,0]
        else:
            nr = len(r)
            err2 = np.zeros((nr))
            for i,r_i in enumerate(r):
                deriv_vect = np.matrix(self.rho_deriv(r_i,self.deriv_params))
                err2[i] = deriv_vect * cov_mat * deriv_vect.T
        return np.sqrt(err2)

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
    
    def _covar(self):
        """ Default implementation of covariance matrix, 
        
        This just uses the parameter errors and ignores the off-diagonal terms
        """
        npar = len(self.deriv_params)
        m = np.matrix(np.zeros((npar,npar)))
        for i,pname in enumerate(self.deriv_params):
            par_err = self.getp(pname).symmetric_error
            m[i,i] = par_err*par_err
        return m 

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
        return 4*np.pi/3 * self.rhos * np.where(r < self.rs, r**3, self.rs**3)
   
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
            ('alpha',     Parameter(default=0.17)),
        ])

    @property 
    def deriv_params(self):
        return ["rs","rhos","alpha"]
    
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
            ('gamma',     Parameter(default=1.)),
        ])

    @property 
    def deriv_params(self):
        return ["rs","rhos","gamma"]  

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
            ('alpha',     Parameter(default=1.)),
            ('beta',      Parameter(default=3.)),
            ('gamma',     Parameter(default=1.)),
        ])

    @property 
    def deriv_params(self):
        return ["rs","rhos","alpha","beta","gamma"]

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



def scale_list(l,scale_value):
    for i,v in enumerate(l):
        l[i] = v*scale_value
    return l

def scale_dict(d,scale_value):
    for k,v in d.items():
        if isinstance(v,list):
            d[k] = scale_list(v,scale_value)
        else:
            d[k] = v*scale_value
    return d

def scale_param(p,scale_value):
    if isinstance(p,dict):
        return scale_dict(p,scale_value)
    elif isinstance(p,list):
        return scale_list(p,scale_value)
    else:
        return p*scale_value

def scale_dict_param(d,k,scale_value, default_value):
    try: 
        d[k] = scale_param( d[k], scale_value )
    except KeyError:
        d[k] = scale_value * default_value


def factory(ptype, **kwargs):
    import dmsky.factory

    prof_copy = kwargs.copy()
    units = prof_copy.pop('units',None)
    if units:
        density,distance = units.rsplit('_',1)      
        scale_density = getattr(Units,density)
        scale_distance = getattr(Units,distance)

        scale_dict_param( prof_copy, 'rhos', scale_density, DensityProfile._params['rhos'].default )
        scale_dict_param( prof_copy, 'rs',   scale_distance, DensityProfile._params['rs'].default )
        scale_dict_param( prof_copy, 'rmin', scale_distance, DensityProfile._params['rmin'].default )
        scale_dict_param( prof_copy, 'rmax', scale_distance, DensityProfile._params['rmax'].default )
        scale_dict_param( prof_copy, 'rhomax', scale_density, DensityProfile._params['rhomax'].default )
            
    return dmsky.factory.factory(ptype, module=__name__,**prof_copy)

