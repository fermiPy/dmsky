#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"

import copy
import numpy as np

from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
import scipy.special as spfn
import scipy.optimize as opt

        
class LoSAnnihilate(object):
    """Integrand function for LoS parameter (J).  The parameter alpha
    introduces a change of coordinates x' = x^(1/alpha).  The change
    of variables means that we need make the substitution:

    dx = alpha * (x')^(alpha-1) dx'

    A value of alpha > 1 weights the points at which we sample the
    integrand closer to x = 0 (distance of closest approach).
    
    ADW: I think this is a luminosity density...
    """

    def __init__(self,dp,d,xi,alpha=4.0):
        """
        Parameters
        ----------
        d: Distance to halo center.
        xi: Offset angle in radians.
        dp: Density profile.
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


class LoSDecay(LoSAnnihilate):
    def __init__(self,dp,d,xi,alpha=1.0):
        super(LoSDecay,self).__init__(dp,d,xi,alpha)
        
    def func(self,r):
        return self.dp.rho(r)


class LoSIntegral(object):
    """
    Slowest (and most accurate) LoS integral.
    """
    def __init__(self, density, distance, alpha=3.0, ann=True):
        self.dp = density
        self.dhalo = distance
        self.alpha = alpha
        self.ann = ann

    @property
    def rmax(self):
        return self.dp.rmax

    def __call__(self, psi, dhalo=None, degrees=False):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------
        psi : array_like 
        Array of offset angles (in radians)
        """
        scalar = np.isscalar(psi)
        if degrees: psi = np.deg2rad(psi)
        psi = np.atleast_1d(psi)
        
        if dhalo is None: dhalo = self.dhalo
        dhalo = np.atleast_1d(dhalo)
            
        if dhalo.size == 1:
            dhalo = dhalo * np.ones_like(psi)
        else: 
            if dhalo.size != psi.size:
                msg = "Array sizes must match"
                raise ValueError(msg)

        v = self._integrate(psi,dhalo)

        if scalar:
            return v[0]
        else:
            return v

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
        #rmin = np.where(psi < np.pi/2, dhalo*np.sin(psi), dhalo)

        for i, (_psi,_dhalo,_rmin) in enumerate(zip(psi,dhalo,rmin)):
            s0,s1 = 0,0
            
            if self.ann:
                losfn = LoSAnnihilate(self.dp,_dhalo,_psi,self.alpha)
            else:
                losfn = LoSDecay(self.dp,_dhalo,_psi,self.alpha)

            # Observer inside the halo...
            if self.rmax > _dhalo:
                if _psi < np.pi/2.:
                    # Distance between observer and point of closest approach
                    x0 = (self.dhalo*np.cos(_psi))**(1./self.alpha)
                    s0 = 2*quad(losfn,0.0,x0)[0]

                    # Distance from point of closest approach to maximum radius
                    x1 = np.sqrt(self.rmax**2 - _rmin**2)**(1./self.alpha)
                    s1 = quad(losfn,x0,x1)[0]
                else:
                    x0 = np.abs(self.dhalo*np.cos(_psi))**(1./self.alpha)
                    x1 = np.sqrt(self.rmax**2 - _rmin**2)**(1./self.alpha)
                                          
                    s1 = quad(losfn,x0,x1)[0]
            # Observer outside the halo...
            elif (self.rmax > _rmin) & (_psi < np.pi/2.):
                x0 = np.sqrt(self.rmax**2 - _rmin**2)**(1./self.alpha)
                s0 = 2*quad(losfn,0.0,x0)[0]
                
            v[i] = s0+s1
        return v


class LoSIntegralFast(LoSIntegral): 
    """Vectorized version of LoSIntegral that performs midpoint
    integration with a fixed number of steps.
    """

    def __init__(self, density, distance, alpha=3.0, ann=True, nsteps=400):
        """
        Parameters
        ----------
        dist: Distance to halo center.
        dp:   Density profile.
        alpha: Parameter determining the integration variable: x' = x^(1/alpha)
        rmax: Radius from center of halo at which LoS integral is truncated.
        nstep: Number of integration steps.  Increase this parameter to
        improve the accuracy of the LoS integral.
        """
        super(LoSIntegralFast,self).__init__(density,distance,alpha,ann)

        self.nsteps = nsteps
        xedge = np.linspace(0,1.0,self.nsteps+1)
        self.x = 0.5*(xedge[1:] + xedge[:-1])

    @property
    def rmax(self):
        if self.dp.rmax < np.inf:
            return self.dp.rmax
        else:
            return 100*self.dp.rs

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

                #s0 = 2*np.sum(losfn(x0)*dx0,axis=0)
                #s1 = np.sum(losfn(x1)*dx1,axis=0)
                s0 = 2*np.apply_over_axes(np.sum,losfn(x0)*dx0,axes=[-1])
                s1 = np.apply_over_axes(np.sum,losfn(x1)*dx1,axes=[-1])

                v[msk01] = s0[msk01]+s1[msk01]

            if np.any(msk02):
            
                dx1 = (xlim1-xlim0)/float(self.nsteps)

                x1 = xlim0 + self.x*(xlim1-xlim0)
                #s0 = np.sum(losfn(x1)*dx1,axis=0)
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
        

class LoSIntegralFunc(LoSIntegralFast):
    """ Subclass for interpolation, spline, and file functions.
    """
    def __init__(self, dp, dist, rmax=None, alpha=3.0, ann=True, nstep=400):
        """Create a fast lookup function"""
        super(LoSIntegralFunc,self).__init__(dp, dist, rmax, alpha, ann, nstep)
        self._func = self.create_func(dist)

    def create_profile(self, dhalo, npsi=150):
        dhalo = np.unique(np.atleast_1d(dhalo))
        dp = self._dp

        rmin = dp.rmin if dp.rmin else dp.rs*1e-6
        psimin = np.arctan(rmin/dhalo.max())
        rmax = dp.rmax if dp.rmax < np.inf else 100*dp.rs
        psimax = np.arctan(rmax/dhalo.min()) if rmax > dhalo.max() else np.pi
        psi = np.logspace(np.log10(psimin),np.log10(psimax),npsi)
        # Pin the last point to 180 deg
        #psi[-1] = np.pi
        
        _dhalo, _psi = np.meshgrid(dhalo,psi)
        _jval = super(LoSIntegralFunc,self).__call__(_psi,_dhalo)
        return _dhalo, _psi, _jval

    def create_func(self, dhalo):
        """Create the spline function
        """
        _dhalo,_psi,_jval = self.create_profile(dhalo)
        log_dhalo,log_psi,log_jval = np.clip(np.log10([_dhalo,_psi,_jval]),-32,None)
        scalar = (_dhalo.shape[-1] == 1)

        #kx = ky = k = 2
        #s=0
        if scalar:
            print "Univariate"
            spline=UnivariateSpline(log_psi.flat,log_jval.flat,k=2,s=0)
            fn = lambda psi: 10**(spline(np.log10(psi)))
        else:
            print "Bivariate"
            spline = bisplrep(log_psi,log_dhalo,log_jval,s=0.0,kx=2,ky=2)
            fn = lambda psi,dhalo: 10**bisplev(np.log10(psi[:,0]),np.log10(dhalo[0,:]),spline)
            #spline=SmoothBivariateSpline(log_dhalo.flat,log_psi.flat,log_jval.flat,
            #                             kx=1,ky=1)
            #fn = lambda dhalo, psi: 10**(spline.ev(np.log10(dhalo),np.log10(psi)))

        return fn,spline

    def __call__(self,psi,dhalo=None):
        """Compute the LoS integral from an interpolating function.

        Returns
        -------
        vals: LoS amplitude per steradian.
        """
        if dhalo is None or np.all(np.in1d(dhalo,self._dist)):
            func = self._func
        else:
            func = self.create_func(dhalo)

        dhalo = np.atleast_1d(dhalo)
        psi = np.atleast_1d(psi)
        if len(dhalo) == 1:
            return func(psi)
        else:
            return func(psi,dhalo)

    def _integrate(self,psi,dhalo):
        # Arrays must be the same shape
        if  psi.shape != dhalo.shape:
            msg = "Shape of psi and dhalo must match"
            raise ValueError(msg)

        if (np.unique(psi) != psi[:,0]).all():
            msg = "np.unique(psi) != psi[:,0]"
            raise ValueError(msg)

        if (np.unique(dhalo) != dhalo[0,:]).all():
            msg = "np.unique(dhalo) != dhalo[0,:]"
            raise ValueError(msg)
            
        if ((np.max(self.dhalo)>=dhalo) & (np.min(self.dhalo)<=dhalo)).all():
            func = self._func
        else:
            func = self.create_func(dhalo)

        if dhalo.size == 1:
            v = func(psi)
        else:
            v = func(psi,dhalo)

        if v.shape != psi.shape:
            msg = "Input and output shape do not match"
            raise ValueError(msg)

        return v

class LoSIntegralSpline(LoSIntegralFast): 
    pass


class LoSIntegralFile(LoSIntegral): 
    def __init__(self, filename):
        pass

    def __call__(self, psi):
        pass



if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()
