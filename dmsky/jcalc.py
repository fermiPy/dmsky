#!/usr/bin/env python

"""
@file   jcalc.py

@brief Python modules that are used to compute the line-of-sight
integral over a spherically symmetric DM distribution.

@author Matthew Wood       <mdwood@slac.stanford.edu>
@author Alex Drlica-Wagner <kadrlica@fnal.gov>
"""

__author__   = "Matthew Wood"
__date__     = "12/01/2011"

import copy
import numpy as np

from scipy.integrate import quad
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from scipy.interpolate import interp1d, UnivariateSpline
import scipy.special as spfn
import scipy.optimize as opt

class LoSFn(object):
    """Integrand function for LoS parameter (J).  The parameter alpha
    introduces a change of coordinates x' = x^(1/alpha).  The change
    of variables means that we need make the substitution:

    dx = alpha * (x')^(alpha-1) dx'

    A value of alpha > 1 weights the points at which we sample the
    integrand closer to x = 0 (distance of closest approach).

    Parameters
    ----------
    d: Distance to halo center.
    xi: Offset angle in radians.
    dp: Density profile.
    alpha: Rescaling exponent for line-of-sight coordinate.
    """

    def __init__(self,d,xi,dp,alpha=4.0):
        self._d = d
        self._d2 = d*d
        self._xi = xi
        self._sinxi = np.sin(xi)
        self._sinxi2 = np.power(self._sinxi,2)
        self._dp = dp
        self._alpha = alpha

    def __call__(self,xp):

        x = np.power(xp,self._alpha)
        r = np.sqrt(x*x+self._d2*self._sinxi2)
        rho2 = np.power(self._dp.rho(r),2)
        return rho2*self._alpha*np.power(xp,self._alpha-1.0)
    
class LoSFnDecay(LoSFn):
    def __init__(self,d,xi,dp,alpha=1.0):
        super(LoSFnDecay,self).__init__(d,xi,dp,alpha)
        
    def __call__(self,xp):
        #xp = np.asarray(xp)
        #if xp.ndim == 0: xp = np.array([xp])

        x = np.power(xp,self._alpha)
        r = np.sqrt(x*x+self._d2*self._sinxi2)
        rho = self._dp.rho(r)
        return rho*self._alpha*np.power(xp,self._alpha-1.0)

class LoSIntegralFn(object):
    """Object that computes integral over DM density squared along a
    line-of-sight offset by an angle psi from the center of the DM
    halo.  We introduce a change of coordinates so that the integrand
    is more densely sampled near the distance of closest of approach
    to the halo center.

    Parameters
    ----------
    dist: Distance to halo center.
    dp: Density profile.
    alpha: Parameter determining the integration variable: x' = x^(1/alpha)
    rmax: Radius from center of halo at which LoS integral is truncated.
    """
    def __init__(self, dp, dist, rmax=None, alpha=3.0,ann=True):
        if rmax is None: rmax = np.inf

        self._dp = dp
        self._dist = dist
        self._rmax = rmax
        self._alpha = alpha
        self._ann = ann

    def __call__(self,psi,dhalo=None):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------
        psi : array_like 
        Array of offset angles (in radians)

        dhalo : array_like
        Array of halo distances.
        """

        if dhalo is None: dhalo = np.array(self._dist,ndmin=1)
        else: dhalo = np.array(dhalo,ndmin=1)

        psi = np.array(psi,ndmin=1)

        if dhalo.shape != psi.shape:
            dhalo = dhalo*np.ones(shape=psi.shape)

        v = np.zeros(shape=psi.shape)

        for i, t in np.ndenumerate(psi):

            s0 = 0
            s1 = 0

            if self._ann:
                losfn = LoSFn(dhalo[i],t,self._dp,self._alpha)
            else:
                losfn = LoSFnDecay(dhalo[i],t,self._dp,self._alpha)

            # Closest approach to halo center
            rmin = dhalo[i]*np.sin(psi[i])

            # If observer inside the halo...
            if self._rmax > dhalo[i]:

                if psi[i] < np.pi/2.:

                    x0 = np.power(dhalo[i]*np.cos(psi[i]),1./self._alpha)
                    s0 = 2*quad(losfn,0.0,x0)[0]

                    x1 = np.power(np.sqrt(self._rmax**2 -
                                          rmin**2),1./self._alpha)
                
                    s1 = quad(losfn,x0,x1)[0]
                else:
                    x0 = np.power(np.abs(dhalo[i]*np.cos(psi[i])),
                                  1./self._alpha)

                    x1 = np.power(np.sqrt(self._rmax**2 -
                                          rmin**2),1./self._alpha)
                    s1 = quad(losfn,x0,x1)[0]

            # If observer outside the halo...
            elif self._rmax > rmin:
                x0 = np.power(np.sqrt(self._rmax**2 -
                                      rmin**2),1./self._alpha)
                s0 = 2*quad(losfn,0.0,x0)[0]
                
            v[i] = s0+s1

        return v

class LoSIntegralFnFast(LoSIntegralFn):
    """Vectorized version of LoSIntegralFn that performs midpoint
    integration with a fixed number of steps.

    Parameters
    ----------
    dist: Distance to halo center.
    dp:   Density profile.
    alpha: Parameter determining the integration variable: x' = x^(1/alpha)
    rmax: Radius from center of halo at which LoS integral is truncated.
    nstep: Number of integration steps.  Increase this parameter to
    improve the accuracy of the LoS integral.
    """
    def __init__(self, dp, dist, rmax=None, alpha=3.0,ann=True,nstep=400):
        super(LoSIntegralFnFast,self).__init__(dp,dist,rmax,alpha,ann)

        self._nstep = nstep
        xedge = np.linspace(0,1.0,self._nstep+1)
        self._x = 0.5*(xedge[1:] + xedge[:-1])

    def __call__(self,psi,dhalo=None):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------
        psi : array_like 
        Array of offset angles (in radians)

        dhalo : array_like
        Array of halo distances.
        """

        if dhalo is None: dhalo = np.array(self._dist,ndmin=1)
        else: dhalo = np.array(dhalo,ndmin=1)

        psi = np.array(psi,ndmin=1)
        v = psi*dhalo
        v.fill(0)
        v = v.reshape(v.shape + (1,))

        dhalo = np.ones(v.shape)*dhalo[...,np.newaxis]
        psi = np.ones(v.shape)*psi[...,np.newaxis]

        if self._ann: losfn = LoSFn(dhalo,psi,self._dp,self._alpha)
        else: losfn = LoSFnDecay(dhalo,psi,self._dp,self._alpha)

        # Closest approach to halo center
        rmin = dhalo*np.sin(psi)

        msk0 = self._rmax > dhalo
        msk1 = self._rmax > rmin

        # Distance between observer and point of closest approach
        xlim0 = np.power(np.abs(dhalo*np.cos(psi)),1./self._alpha)

        # Distance from point of closest approach to maximum
        # integration radius
        xlim1 = np.zeros(shape=psi.shape)
        xlim1[msk1] = np.power(np.sqrt(self._rmax**2 - rmin[msk1]**2),
                               1./self._alpha)

        # If observer inside the halo...
        if np.any(msk0):

            msk01 = msk0 & (psi < np.pi/2.)
            msk02 = msk0 & ~(psi < np.pi/2.)

            if np.any(msk01):

                dx0 = xlim0/float(self._nstep)
                dx1 = (xlim1-xlim0)/float(self._nstep)

                x0 = self._x*xlim0
                x1 = xlim0 + self._x*(xlim1-xlim0)

#                s0 = 2*np.sum(losfn(x0)*dx0,axis=0)
#                s1 = np.sum(losfn(x1)*dx1,axis=0)
                s0 = 2*np.apply_over_axes(np.sum,losfn(x0)*dx0,axes=[-1])
                s1 = np.apply_over_axes(np.sum,losfn(x1)*dx1,axes=[-1])

                v[msk01] = s0[msk01]+s1[msk01]

            if np.any(msk02):
            
                dx1 = (xlim1-xlim0)/float(self._nstep)

                x1 = xlim0 + self._x*(xlim1-xlim0)
#                s0 = np.sum(losfn(x1)*dx1,axis=0)
                s0 = np.apply_over_axes(np.sum,losfn(x1)*dx1,axes=[-1])
            
                v[msk02] = s0[msk02]
                
        # If observer outside the halo...
        if np.any(~msk0 & msk1):
            
            dx0 = xlim1/float(self._nstep)
            x0 = xlim1*self._x
            s0 = 2*np.apply_over_axes(np.sum,losfn(x0)*dx0,axes=[-1])

            v[~msk0 & msk1] = s0[~msk0 & msk1]

        v = v.reshape(v.shape[:-1])
        return v

class LoSIntegralSplineFn(object):

    def __init__(self,dp=None,nx=40,ny=20):
        self.dp = copy.copy(dp)

        if self.dp is not None:
            nx = 40
            ny = 20
            dhalo, psi = np.mgrid[1:2:ny*1j,0.001:2.0:nx*1j]
            dhalo = np.power(10,dhalo)
            psi = np.radians(psi)            
            f = LoSIntegralFn(self.dp)
            self.z = f(dhalo,psi)
            self.init_spline(dhalo,psi,self.z)

    def init_spline(self,dhalo,psi,z):
        """Compute knots and coefficients of an interpolating spline
        given a grid of points in halo distance (dhalo) and offset
        angle (psi) at which the LoS integral has been computed.
        """

        kx = 2
        ky = 2
        self._psi_min = psi.min()
        self._tck = bisplrep(dhalo,psi,np.log10(z),s=0.0,kx=kx,ky=ky,
                             nxest=int(kx+np.sqrt(len(z.flat))),
                             nyest=int(ky+np.sqrt(len(z.flat))))

    def __call__(self,dhalo,psi,rho=1,rs=1):
        """Compute the LoS integral using a 2D spline table.

        Returns
        -------

        vals: LoS amplitude per steradian.
        """

        dhalo = np.asarray(dhalo)
        psi = np.asarray(psi)

        if dhalo.ndim == 0: dhalo = np.array([dhalo])
        if psi.ndim == 0: psi = np.array([psi])

        if psi.ndim == 2 and dhalo.ndim == 2:
            v = np.power(10,bisplev(dhalo[:,0],psi[0,:],self._tck))
        else:
            v = np.power(10,bisplev(dhalo,psi,self._tck))

        v *= rho*rho*rs
        return v


def SolidAngleIntegral(psi,pdf,angle):
    """ Compute the solid-angle integrated j-value
    within a given radius

    Parameters
    ----------
    psi : array_like 
    Array of offset angles (in radians)

    pdf : array_like
    Array of j-values at angle psi
    
    angle : array_like
    Maximum integration angle (in degrees)
    """
    angle = np.asarray(angle)
    if angle.ndim == 0: angle = np.array([angle])

    scale=max(pdf)
    norm_pdf = pdf/scale
    bad = np.where(norm_pdf <= 0)[0]
    idx = bad.min() if bad.size else len(pdf)
    log_spline = UnivariateSpline(psi[:idx],np.log10(norm_pdf[:idx]),k=1,s=0)
    spline = lambda r: 10**(log_spline(r))
    integrand = lambda r: spline(r)*2*np.pi*np.sin(r)

    integral = []
    for a in angle:
        integral.append(quad(integrand, 0, np.radians(a),full_output=True)[0])
    integral = np.asarray(integral)

    return integral*scale

class JProfile(object):
    def __init__(self,losfn):

        log_psi = np.linspace(np.log10(np.radians(1E-5)),
                              np.log10(np.radians(90.)),1001)
        self._psi = np.power(10,log_psi)
        self._psi = np.insert(self._psi,0,0)

        domega = 2*np.pi*(-np.cos(self._psi[1:])+np.cos(self._psi[:-1]))
        x = 0.5*(self._psi[1:]+self._psi[:-1])

        self._jpsi = losfn(x)
        self._spline = UnivariateSpline(x,self._jpsi,s=0,k=1)
        self._jcum = np.cumsum(self._spline(x)*domega)
        self._jcum = np.insert(self._jcum,0,0)

        self._cum_spline = UnivariateSpline(self._psi,self._jcum,s=0,k=2)

    @staticmethod
    def create(dp,dist,rmax):
        losfn = LoSIntegralFn(dp,dist,rmax=rmax)        
        return JProfile(losfn)

    def __call__(self,psi):
        return self._spline(psi)

    def integrate(self,psimax):

        xedge = np.linspace(0.0,np.radians(psimax),1001)
        x = 0.5*(xedge[1:] + xedge[:-1])
        domega = 2.0*np.pi*(-np.cos(xedge[1:])+np.cos(xedge[:-1]))
        return np.sum(self._spline(x)*domega)

    def cumsum(self,psi):
        return self._cum_spline(np.radians(psi))

class DensityProfile(object):
    """ DM density profile that truncates at a maximum DM density.
    
    rho(r) = rho(r) for rho(r) < rhomax AND r > rmin
           = rhomax for rho(r) >= rhomax
           = rho(rmin) for r <= rmin
    
    Parameters
    ----------
    rhos : Density normalization parameter.
    
    rmin : Inner radius interior to which the density will be fixed to
    a constant value. (rhomax = rho(rmin)).

    rhomax : Maximum DM density.  If rhomax and rmin are both defined
    the maximum DM density will be the lesser of rhomax and rho(rmin).
    
    """
    def __init__(self,rhos,rmin=None,rhomax=None):
        self._name = 'profile'
        self._rmin=rmin
        self._rhomax=rhomax
        self._rhos = rhos

    def setMassConcentration(self,mvir,c):

        rhoc = 9.9E-30*Units.g_cm3
        rvir = np.power(mvir*3.0/(177.7*4*np.pi*rhoc*0.27),1./3.)
        rs = rvir/c

        self._rs = rs
        mrvir = self.mass(rvir)
        self._rhos = self._rhos*mvir/mrvir

    def rho(self,r):

        r = np.asarray(r)
        if r.ndim == 0: r = r.reshape((1))

        if self._rhomax is None and self._rmin is None: 
            return self._rho(r)
        elif self._rhomax is None:
            rho = self._rho(r)        
            rho[r<self._rmin] = self._rho(self._rmin)
            return rho
        elif self._rmin is None:
            rho = self._rho(r)        
            rho[rho>self._rhomax] = self._rhomax
            return rho
        else:
            rho = self._rho(r) 
            rhomax = min(self._rho(self._rmin),self._rhomax)
            rho[rho>rhomax] = rhomax
            return rho

#            return np.where(rho>self._rhomax,[self._rhomax],rho)
        
    def set_rho(self,rho,r):
        """Fix the density normalization at a given radius."""
        rhor = self._rho(r)
        self._rhos = rho*self._rhos/rhor
        
    def name(self):
        return self._name

    @staticmethod
    def create(**kwargs):
        """Method for instantiating a density profile object given the
        profile name and a dictionary."""

        kwargs.setdefault('type','nfw')
        kwargs.setdefault('rhos','1.0')

        o = dict(kwargs)
        name = o['type'].lower()

        def extract(keys,d):
            od = {}
            for k in keys: 
                if k in d: od[k] = d[k]
            return od

        if name == 'nfw':
            dp = NFWProfile(**extract(['rhos','rs','rmin'],o))
        elif name == 'gnfw':
            dp = GNFWProfile(**extract(['rhos','rs','rmin','gamma'],o))
        elif name == 'isothermal':
            dp = IsothermalProfile(**extract(['rhos','rs','rmin'],o))
        elif name == 'einasto':
            dp = EinastoProfile(**extract(['rhos','rs','rmin','alpha'],o))
        elif name == 'burkert':
            dp = BurkertProfile(**extract(['rhos','rs','rmin'],o))
        else:
            print 'No such halo type: ', name
            #sys.exit(1)

        if 'rhor' in o:
            dp.set_rho(o['rhor'][0]*Units.gev_cm3,
                       o['rhor'][1]*Units.kpc)
        elif 'jval' in o:
            dp.set_jval(o['jval']*Units.gev2_cm5,
                        o['rs'],
                        o['dist'])

        return dp


class BurkertProfile(DensityProfile):
    """ Burkert (1995)
        rho(r) = rhos/( (1+r/rs)(1+(r/rs)**2) )
    """
    def __init__(self,rhos=1,rs=1,rmin=None,rhomax=None):        
        self._rs = rs
        super(BurkertProfile,self).__init__(rhos,rmin,rhomax)
        self._name = 'burkert'

    def _rho(self,r):
        x = r/self._rs
        return self._rhos*np.power(1+x,-1)*np.power(1+x*x,-1)

    def mass(self,r):
        x = r/self._rs        
        return np.pi*self._rhos*(np.log(x**2+1)+2*np.log(x+1)-2*np.arctan(x))

class IsothermalProfile(DensityProfile):
    """ Isothermal Profile
        rho(r) = rhos/(1+(r/rs))**2
    """

    def __init__(self,rhos=1,rs=1,rmin=None,rhomax=None):        
        self._rs = rs
        super(IsothermalProfile,self).__init__(rhos,rmin,rhomax)
        self._name = 'isothermal'

    def _rho(self,r):
        x = r/self._rs
        return self._rhos*np.power(1+x,-2)

    def mass(self,r):
        raise Exception("IsothermalProfile mass not implemented")
    
class NFWProfile(DensityProfile):
    """ Navarro, Frenk, and White (1996)
        rho(r) = rhos/( (r/rs)(1+r/rs)**2)
    """
    def __init__(self,rhos=1,rs=1,rmin=None,rhomax=None):
        self._rs = rs
        super(NFWProfile,self).__init__(rhos,rmin,rhomax)
        self._name = 'nfw'

    def set(self,rhos,rs):
        self._rs = rs
        self._rhos = rhos

    def set_jval(self,jval,rs,dist):
        rhos = np.sqrt(3./(4.*np.pi)*jval*dist**2/rs**3)
        self._rs = rs
        self._rhos = rhos

    def mass(self,r):
        x = r/self._rs
        return 4*np.pi*self._rhos*np.power(self._rs,3)*(np.log(1+x)-x/(1+x))

    def jval(self,r=None,rhos=None,rs=None):

        if rhos is None: rhos = self._rhos
        if rs is None: rs = self._rs

        if r is not None:
            x = r/rs
            return (4*np.pi/3.)*rhos**2*rs**3*(1.-np.power(1.+x,-3))
        else:
            return (4*np.pi/3.)*rhos**2*rs**3

#(4*M_PI/3.)*std::pow(a(0),2)*std::pow(a(1),3)*(1.-std::pow(1+x,-3));

    def _rho(self,r):
        x = r/self._rs
        return self._rhos*np.power(x,-1)*np.power(1+x,-2)        
    
class EinastoProfile(DensityProfile):
    """ Einasto profile
        rho(r) = rhos*exp(-2*((r/rs)**alpha-1)/alpha)
    """
    def __init__(self,rhos=1,rs=1,alpha=0.17,rmin=None,rhomax=None):
        self._rs = rs
        self._alpha = alpha
        super(EinastoProfile,self).__init__(rhos,rmin,rhomax)
        self._name = 'einasto'

    def set(self,rhos,rs):
        self._rs = rs
        self._rhos = rhos

    def mass(self,r):

        x = r/self._rs
        gamma = spfn.gamma(3./self._alpha)

        return 4*np.pi*self._rhos*np.power(self._rs,3)/self._alpha* \
            np.exp(2./self._alpha)* \
            np.power(2./self._alpha,-3./self._alpha)* \
            gamma*spfn.gammainc(3./self._alpha,
                                (2./self._alpha)*np.power(x,self._alpha))

    def _rho(self,r):
        x = r/self._rs
        return self._rhos*np.exp(-2./self._alpha*(np.power(x,self._alpha)-1))

class GNFWProfile(DensityProfile):
    """ Generalized NFW Profile
        rho(r) = rhos/( (r/rs)^g(1+r/rs)**(3-g))
    """
    def __init__(self,rhos=1,rs=1,gamma=1.0,rmin=None,rhomax=None):
        self._rs = rs
        self._gamma = gamma
        super(GNFWProfile,self).__init__(rhos,rmin,rhomax)
        self._name = 'nfw'

    def set(self,rhos,rs):
        self._rs = rs
        self._rhos = rhos

    def mass(self,r):
#        x = r/self._rs
#        return 4*np.pi*self._rhos*np.power(self._rs,3)*(np.log(1+x)-x/(1+x))
        return 0

    def _rho(self,r):
        x = r/self._rs
        return self._rhos*np.power(x,-self._gamma)* \
            np.power(1+x,-(3-self._gamma))    

class GeneralNFWProfile(DensityProfile):
    """ Strigari et al. (2007)
        rho(r) = rhos/( (r/rs)**a (1+(r/rs)**b )**(c-a)/b
        Default: NFW profile
    """
    def __init__(self,rhos=1,rs=1,a=1,b=1,c=3,rmin=None,rhomax=None):
        self._rs = rs
        self._a = a
        self._b = b
        self._c = c
        super(GeneralNFWProfile,self).__init__(rhos,rmin,rhomax)
        self._name = 'general_nfw'

    def _rho(self,r):
        x = r/self._rs
        return self._rhos/(x**self._a*(1+x**self._b)**((self._c-self._a)/self._b))

class NFWcProfile(GNFWProfile):
    """ Contracted NFW profile
    ADW: UNTESTED
    """
    def __init__(self,rhos=1,rs=1,rmin=None,rhomax=None):
        gamma = 1.3
        super(NFWcProfile,self).__init__(rhos,rs,gamma,rmin,rhomax)
        self._name = 'nfwc'
    

class UniformProfile(object):
    """ Uniform spherical profile
        rho(r) = rhos for r < rs
        rho(r) = 0    otherwise
    """
    def __init__(self,rhos=1,rs=1):
        self._name = 'uniform'
        self._rhos = rhos
        self._rs = rs

    def _rho(self,r):
        return np.where(r<rs,rhos,0)

class Units(object):
    pc = 3.08568e18   # pc to cm
    kpc = pc*1e3      # kpc to cm
    msun = 1.98892e33 # solar mass to g
    gev = 1.78266e-24 # gev to g
    g = 1.0
    m2 = 1E4
    hr = 3600.
    deg2 = np.power(np.pi/180.,2)

    msun_pc3 = msun*np.power(pc,-3) 
    msun_kpc3 = msun*np.power(kpc,-3)
    msun2_pc5 = np.power(msun,2)*np.power(pc,-5)
    msun2_kpc5 = np.power(msun,2)*np.power(kpc,-5)
    gev2_cm5 = np.power(gev,2)
    gev_cm3 = np.power(gev,1)
    gev_cm2 = np.power(gev,1)
    g_cm3 = 1.0
    cm3_s = 1.0

if __name__ == '__main__':
    print "Line-of-sight Integral Package..."

    import matplotlib.pyplot as plt

    psi = np.linspace(0.01,0.1,500)
    dp = NFWProfile(1,1)

    fn0 = LoSIntegralFnFast(dp,100,10)
    fn1 = LoSIntegralFn(dp,100,10)

    dhalo = np.linspace(100,100,500)
    v0 = fn0(dhalo,psi)
    v1 = fn1(dhalo,psi)

    delta = (v1-v0)/v0
    print delta

    plt.ion()
    plt.hist(delta,bins=100)
    plt.show()
    
