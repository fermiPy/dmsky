#!/usr/bin/env python
"""
Python module for computing the line-of-sight integral over a
spherically symmetric distribution.

"""
__author__   = "Matthew Wood"
__date__     = "12/01/2011"

import copy
import numpy as np

from scipy.integrate import quad
from scipy.interpolate import bisplrep,bisplev
from scipy.interpolate import interp1d,interp2d
from scipy.interpolate import UnivariateSpline,SmoothBivariateSpline
import scipy.special as spfn
import scipy.optimize as opt

class LoSFn(object):
    """Integrand function (luminosity density) for LoS integration.  The
    parameter `alpha` introduces a change of variables:

    x' = x^(1/alpha). 

    A value of alpha > 1 samples the integrand closer to x = 0
    (distance of closest approach). The change of variables requires
    the substitution:

    dx = alpha * (x')^(alpha-1) dx'

    """

    def __init__(self,dp,d,xi,alpha=3.0):
        """
        Parameters
        ----------
        dp:    Density profile.
        d:     Distance to halo center.
        xi:    Offset angle in radians.
        alpha: Rescaling exponent for line-of-sight coordinate.
        """
        self.dp = dp
        self.d2 = d**2
        self.sinxi2 = np.sin(xi)**2
        self.alpha = alpha

    def __call__(self, xp, alpha=None):
        """
        xp: Distance along the LoS.
        """
        if alpha is None: alpha = self.alpha
        x = xp**alpha
        r = np.sqrt(x**2+self.d2*self.sinxi2)
        return self.func(r)*alpha*xp**(alpha-1.0)

    def func(self, r):
        return self.dp.rho(r)**2

class LoSAnnihilate(LoSFn):
    """Integrand function for LoS annihilation (J-factor)."""

    def __init__(self,dp,d,xi,alpha=3.0):
        super(LoSAnnihilate,self).__init__(dp,d,xi,alpha)

    def func(self, r):
        return self.dp.rho(r)**2

class LoSDecay(LoSAnnihilate):
    """Integrand function for LoS decay (D-factor)."""

    def __init__(self,dp,d,xi,alpha=1.0):
        super(LoSDecay,self).__init__(dp,d,xi,alpha)
        
    def func(self,r):
        return self.dp.rho(r)

class LoSIntegral(object):
    """Slowest (and most accurate?) LoS integral. Uses
    scipy.integrate.quad with a change of variables to better sample
    the LoS close to the halo center.
    """
    
    def __init__(self, density, dhalo, alpha=3.0, ann=True):
        """
        Parameters
        ----------
        density: Density profile.
        dhalo:   Distance to halo center.
        alpha:   Parameter determining the integration variable: x' = x^(1/alpha)
        ann:     Annihilation or decay
        """
        self.dp = density
        self.dhalo = dhalo
        self.alpha = float(alpha)
        self.ann = ann

    @property
    def rmax(self):
        return self.dp.rmax

    @property
    def name(self):
        return self.__class__.__name__

    def __call__(self, psi, dhalo=None, degrees=False):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------
        psi     : Array of offset angles (in radians by default)
        dhalo   : Array of halo distances
        degrees : Interpret `psi` in degrees
        """
        scalar = np.isscalar(psi)
        if degrees: psi = np.deg2rad(psi)
        psi = np.atleast_1d(psi)
        
        if dhalo is None: 
            if self.dhalo is None:
                msg = "Halo distance must be specified"
                raise Exception(msg)
            dhalo = self.dhalo
        dhalo = np.atleast_1d(dhalo)
            
        if dhalo.size == 1:
            dhalo = dhalo * np.ones_like(psi)
        else: 
            if dhalo.size != psi.size:
                msg = "Array sizes must match"
                raise ValueError(msg)

        v = self._integrate(psi,dhalo)

        if scalar: return v[0]
        else:      return v
            

    def _integrate(self, psi, dhalo):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------
        psi:   Array of offset angles (in radians)
        dhalo: Array of halo distances
        """

        # Arrays must be the same shape
        assert (psi.shape == dhalo.shape)

        # Output Array
        v = np.zeros_like(psi)

        # Closest approach to halo center
        rmin = dhalo*np.sin(psi)
        # Maximum extent of halo
        rmax = self.rmax

        for i, t in np.ndenumerate(psi):
            s0,s1 = 0,0
            
            if self.ann:
                losfn = LoSAnnihilate(self.dp,dhalo[i],psi[i],self.alpha)
            else:
                losfn = LoSDecay(self.dp,dhalo[i],psi[i],self.alpha)

            # Closest approach to halo center
            #rmin = dhalo[i]*np.sin(psi[i])

            # If observer inside the halo...
            if rmax > dhalo[i]:

                if psi[i] < np.pi/2.:
                    x0 = np.power(dhalo[i]*np.cos(psi[i]),1./self.alpha)
                    s0 = 2*quad(losfn,0.0,x0)[0]

                    x1 = np.power(np.sqrt(rmax**2 - rmin[i]**2),1./self.alpha)
                    s1 = quad(losfn,x0,x1)[0]
                else:
                    x0 = np.power(np.abs(dhalo[i]*np.cos(psi[i])),1./self.alpha)

                    x1 = np.power(np.sqrt(rmax**2 - rmin[i]**2),1./self.alpha)
                    s1 = quad(losfn,x0,x1)[0]

            # If observer outside the halo...
            elif (rmax > rmin[i]) and (psi[i] < np.pi/2.):
                x0 = np.power(np.sqrt(rmax**2 - rmin[i]**2),1./self.alpha)
                s0 = 2*quad(losfn,0.0,x0)[0]
                
            v[i] = s0+s1

        return v


class LoSIntegralFast(LoSIntegral): 
    """
    Vectorized version of LoSIntegral that performs midpoint
    integration with a fixed number of steps.
    """

    def __init__(self, density, dhalo, alpha=3.0, ann=True, nsteps=400):
        """
        Parameters
        ----------
        density: Density profile.
        dhalo:   Distance to halo center.
        alpha:   Parameter determining the integration variable: x' = x^(1/alpha)
        ann:     Annihilation or decay
        nsteps:  Number of integration steps.  Increase this parameter to
                 improve the accuracy of the LoS integral.
        """
        super(LoSIntegralFast,self).__init__(density,dhalo,alpha,ann)

        self.nsteps = nsteps
        xedge = np.linspace(0,1.0,self.nsteps+1)
        self.x = 0.5*(xedge[1:] + xedge[:-1])

    @property
    def rmax(self):
        if self.dp.rmax < np.inf:
            return self.dp.rmax
        else:
            return 1000*self.dp.rs

    def _integrate(self,psi,dhalo):
        # Arrays must be the same shape
        assert (psi.shape == dhalo.shape)

        v = np.zeros(psi.shape + (1,))

        dhalo = np.ones(v.shape)*dhalo[...,np.newaxis]
        psi = np.ones(v.shape)*psi[...,np.newaxis]

        # Closest approach to halo center
        rmin = dhalo*np.sin(psi)

        if self.ann: losfn = LoSAnnihilate(self.dp,dhalo,psi,self.alpha)
        else:        losfn = LoSDecay(self.dp,dhalo,psi,self.alpha)

        msk0 = self.rmax > dhalo
        msk1 = self.rmax > rmin

        # Distance between observer and point of closest approach
        xlim0 = np.abs(dhalo*np.cos(psi))**(1./self.alpha)

        # Distance from point of closest approach to maximum
        # integration radius
        xlim1 = np.zeros(shape=psi.shape)
        xlim1[msk1] = np.sqrt(self.rmax**2 - rmin[msk1]**2)**(1./self.alpha)

        # If observer inside the halo...
        if np.any(msk0):

            msk01 = msk0 & (psi < np.pi/2.)
            msk02 = msk0 & ~(psi < np.pi/2.)

            if np.any(msk01):
                dx0 = xlim0/float(self.nsteps)
                dx1 = (xlim1-xlim0)/float(self.nsteps)

                x0 = self.x*xlim0
                x1 = xlim0 + self.x*(xlim1-xlim0)

                s0 = 2*np.apply_over_axes(np.sum,losfn(x0)*dx0,axes=[-1])
                s1 = np.apply_over_axes(np.sum,losfn(x1)*dx1,axes=[-1])

                v[msk01] = s0[msk01]+s1[msk01]

            if np.any(msk02):
                dx1 = (xlim1-xlim0)/float(self.nsteps)

                x1 = xlim0 + self.x*(xlim1-xlim0)
                s0 = np.apply_over_axes(np.sum,losfn(x1)*dx1,axes=[-1])
            
                v[msk02] = s0[msk02]
                
        # Observer outside the halo...
        # Only calculate integral for psi < pi/2
        msk11 = (~msk0 & msk1) & (psi < np.pi/2.)
        if np.any(msk11):
            dx0 = xlim1/float(self.nsteps)
            x0 = xlim1*self.x
            s0 = 2*np.apply_over_axes(np.sum,losfn(x0)*dx0,axes=[-1])
                
            v[msk11] = s0[msk11]

        return v.reshape(v.shape[:-1])
        

class LoSIntegralInterp(LoSIntegralFast):
    """ Interpolate fast integral a for even faster look-up. """
    
    def __init__(self, density, dhalo, alpha=3.0, ann=True, nsteps=400):
        """
        Parameters
        ----------
        density: Density profile.
        dhalo:   Distance to halo center.
        alpha:   Parameter determining the integration variable: x' = x^(1/alpha)
        ann:     Annihilation or decay
        nsteps:  Number of integration steps.  Increase this parameter to
                 improve the accuracy of the LoS integral.
        """
        super(LoSIntegralInterp,self).__init__(density, dhalo, alpha, ann, nsteps)
        self.func = self.create_func(self.dhalo)

    def create_profile(self, dhalo, nsteps=None):
        if not nsteps: nsteps = self.nsteps
        dhalo = np.unique(np.atleast_1d(dhalo))
        dp = self.dp

        psi = np.logspace(np.log10(1e-7),np.log10(np.pi),nsteps)

        _dhalo, _psi = np.meshgrid(dhalo,psi)
        _jval = super(LoSIntegralInterp,self)._integrate(_psi,_dhalo)
        return np.log10([_dhalo, _psi, _jval])

    def create_func(self, dhalo):
        """Create the spline function
        """
        log_dhalo,log_psi,log_jval = self.create_profile(dhalo)

        zeroval = -99
        log_jval[np.where(log_jval==-np.inf)] = zeroval

        if log_dhalo.shape[-1] == 1:
            #print 'interp1d'
            #spline=UnivariateSpline(log_psi.flat,log_jval.flat,k=2,s=0)
            #fn = lambda psi: 10**(spline(np.log10(psi)))
            interp = interp1d(log_psi.flat,log_jval.flat,kind='linear')
            def fn(psi,dhalo):
                log_jval = interp(np.log10(psi))
                log_jval[np.where(log_jval < zeroval+1)] = -np.inf
                return 10**log_jval
        else:
            #print 'interp2d'
            #spline = bisplrep(log_psi,log_dhalo,log_jval,s=0.0,kx=2,ky=2)
            #fn = lambda psi,dhalo: 10**bisplev(np.log10(psi[:,0]),np.log10(dhalo[0,:]),spline)
            interp = interp2d(log_psi,log_dhalo,log_jval,kind='linear')
            def fn(psi,dhalo):
                log_jval = interp(np.log10(psi[:,0]),np.log10(dhalo[0,:])).T
                log_jval[np.where(log_jval < zeroval+1)] = -np.inf
                return 10**log_jval

        return fn

    def _integrate(self,psi,dhalo):
        # Arrays must be the same shape
        if  psi.shape != dhalo.shape:
            msg = "Shape of psi and dhalo must match"
            raise ValueError(msg)

        if psi.ndim > 1 and not (np.unique(psi) == psi[:,0]).all():
            msg = "np.unique(psi) != psi[:,0]"
            raise ValueError(msg)

        if dhalo.ndim >1 and not (np.unique(dhalo) == dhalo[0,:]).all():
            msg = "np.unique(dhalo) != dhalo[0,:]"
            raise ValueError(msg)
            
        # All halo distances within pre-existing interpolation range
        if ((np.max(self.dhalo)>=dhalo) & (np.min(self.dhalo)<=dhalo)).all():
            func = self.func
        else:
            func = self.create_func(dhalo)

        v = func(psi,dhalo)

        if v.shape != psi.shape:
            msg = "Input and output shape do not match"
            raise ValueError(msg)

        return v

class LoSIntegralFile(LoSIntegralInterp): 
    """
    Interpolate over a pre-generated file.

    NOT IMPLEMENTED YET
    """
    def __init__(self, dp, dist, filename, ann=True):
        """Create a fast look-up interpolation"""
        super(LoSIntegralInterp,self).__init__(dp, dist, ann=ann)
        self.filename = filename

    def create_profile(self, dhalo, npsi=300):
        log_psi,log_jval = np.loadtxt(filename,unpack=True)
        return self.dhalo, log_psi, log_jval


if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
