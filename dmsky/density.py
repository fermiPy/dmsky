#!/usr/bin/env python

"""
Density profiles need to be set by:
1) rhos and rs
2) Integrated J-factor at a given radius and scale radius
3) Ingegrated J-factor only (need to approximate scale radius)

It would be good to generalize profile shapes, i.e., following the prescription of Zhou (1996).

Careful `rhos` is a normalization parameter and is *NOT* the same as rho(rs).

"""

from collections import OrderedDict as odict

import numpy as np
import scipy.special as spfn

from pymodeler.model import Model
from pymodeler.parameter import Parameter, Derived

from dmsky.utils.units import Units


class DensityProfile(Model):
    """A DM density profile"""
    _params = odict([
        ('rs', Parameter(default=1.0)),
        ('rhos', Parameter(default=1.0)),
        ('rmin', Parameter(default=0.0)),
        ('rmax', Parameter(default=np.inf)),
        ('rhomax', Parameter(default=np.inf)),
        ('covar', Derived(dtype=np.ndarray, help='Covariance matrix for parameters')),
    ])

    def __call__(self, r):
        """Return the density for given radii.

        Parameters
        ----------

        r : `numpy.array` or float
            The radii

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input radii

        """
        return self.rho(r)

    @property
    def deriv_params(self):
        """Return the list of paramters we can take derivatives w.r.t.
        """
        return ["rs", "rhos"]

    def rho(self, r):
        """Return the density for given radii.

        Parameters
        ----------

        r : `numpy.array` or float
            The radii

        Returns
        -------

        values : `numpy.array`
            Return values, same shape as the input radii

        """
        scalar = np.isscalar(r)
        r = np.atleast_1d(r)

        rho = self._rho(r)

        if self.rmin:
            rho[r < self.rmin] = self._rho(self.rmin)
        if self.rmax:
            rho[r > self.rmax] = 0
        if self.rhomax:
            rho[rho > self.rhomax] = self.rhomax

        if scalar:
            return np.asscalar(rho)
        return rho

    def rho_deriv(self, r, paramNames):
        """Return the derivatives of the density as a function of radius,
           w.r.t. a list of parameters

        Parameters
        ----------

        r : `numpy.array` or float
            The radii

        paramNames : list
            The names of the parameters to differentiation w.r.t.

        Returns
        -------

        matrix : `numpy.array`
            An n x m array, where:
            n is the number of radii
            m is the number of parameters

        """
        initParVals = np.array([self.__getattr__(pName) for pName in paramNames])
        deltaParVals = initParVals * 0.001

        init_r = self._rho(r)

        derivs = []

        # loop over parameters and take the numerical derivatives
        for initPar, deltaPar, parName in zip(initParVals, deltaParVals, paramNames):
            par = self.getp(parName)
            newParVal = initPar + deltaPar
            par.set_value(newParVal)
            new_r = self._rho(r)
            dr_dp = (new_r - init_r) / (newParVal - initPar)
            derivs.append(dr_dp)
            par.set_value(initPar)

        ret = np.vstack(derivs)
        return ret

    def rho_uncertainty(self, r):
        """Calculate the uncertainty of the density at given radii

        Parameters
        ----------

        r : `numpy.array` or float
            The radii

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input radii

        """
        cov_mat = self.covar
        if np.isscalar(r):
            nr = 1
            deriv_vect = np.matrix(self.rho_deriv(r, self.deriv_params))
            err2 = (deriv_vect.T * cov_mat * deriv_vect)[0, 0]
        else:
            nr = len(r)
            err2 = np.zeros((nr))
            for i, r_i in enumerate(r):
                deriv_vect = np.matrix(self.rho_deriv(r_i, self.deriv_params))
                err2[i] = deriv_vect * cov_mat * deriv_vect.T
        return np.sqrt(err2)

    def _rho(self, r):
        """Internal function for sub-class to return the density at given radii
        """
        raise NotImplementedError("%s._rho not implemented" % (self.__class__.__name__))

    def mass(self, r=None):
        """Compute the mass of the object out to a particular radius.

        Parameters
        ----------

        r : `numpy.array` or float
            The radii

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input radii

        """
        if r is None:
            r = self.rmax
        scalar = np.isscalar(r)
        r = np.atleast_1d(r)
        mass = self._mass(r)
        if scalar:
            return mass[0]
        return mass

    def _mass(self, r):
        """Internal function for sub-class to compute the mass at given radii
        """
        raise NotImplementedError("%s._mass not implemented" % (self.__class__.__name__))

    def set_rho_r(self, rho, r):
        """Fix the density normalization at a given radius.

        Parameters
        ----------

        rho : float
            The normalization density

        r : float
            The corresponding radius

        """
        rhor = self._rho(r)
        rhos = self.getp('rhos')
        rhos *= (rho / rhor)


    def set_mvir_c(self, mvir, c):
        """Fix the mass inside the virial radius.

        Parameters
        ----------

        mvir : float
            The virial radius

        c : float
            Scale factor

        """
        rhoc = 9.9E-30 * Units.g_cm3
        rvir = np.power(mvir * 3.0 / (177.7 * 4 * np.pi * rhoc * 0.27), 1. / 3.)
        rs_val = rvir / c
        self.setp('rs', value=rs_val)

        mrvir = self.mass(rvir)
        rhos = self.getp('rhos')
        rhos *= mvir / mrvir

    def _covar(self):
        """Default implementation of covariance matrix,

        This just uses the parameter errors and ignores the off-diagonal terms

        Returns
        -------

        covs : `numpy.array`
             n x n matrix, where n is the number of differentiable parameters

        """
        npar = len(self.deriv_params)
        m = np.matrix(np.zeros((npar, npar)))
        for i, pname in enumerate(self.deriv_params):
            par_err = self.getp(pname).symmetric_error
            m[i, i] = par_err * par_err
        return m

    def _cache(self, name=None):
        """Cache any `Derived` paramters that are slow to compute
        """
        pass


class UniformProfile(DensityProfile):
    """ Uniform spherical profile
    rho(r) = rhos for r <= rs
    rho(r) = 0    otherwise
    """

    def _rho(self, r):
        """Internal function for sub-class to return the density at given radii
        """
        x = r / self.rs
        return np.where(x <= 1, self.rhos, 0.0)

    def _mass(self, r):
        """Internal function for sub-class to compute the mass at at given radii
        """
        return 4 * np.pi / 3 * self.rhos * np.where(r < self.rs, r**3, self.rs**3)


class IsothermalProfile(DensityProfile):
    """ Non-Singular Isothermal Profile:
    Begeman et al. MNRAS 249, 523 (1991)
    http://adsabs.harvard.edu/full/1991MNRAS.249..523B
    rho(r) = rhos/(1+(r/rs))**2
    """

    def _rho(self, r):
        """Internal function for sub-class to return the density at given radii
        """
        x = r / self.rs
        return self.rhos * (1 + x)**(-2)

    def _mass(self, r):
        """Internal function for sub-class to compute the mass at at given radii
        """
        x = r / self.rs
        return x - (x + 1)**-1 - 2 * np.log(x + 1)


class BurkertProfile(DensityProfile):
    """Burkert ApJ 447, L25 (1995) [Eqn. 2]
    http://arxiv.org/abs/astro-ph/9504041
    rho(r) = rho0 * r0**3 / ( (r + r0)*(r**2+r0**2) )
    ==>
    rho(r) = rhos / ( (1+r/rs)*(1+(r/rs)**2) )
    """

    def _rho(self, r):
        """Internal function for sub-class to return the density at given radii
        """
        x = r / self.rs
        return self.rhos * ((1 + x) * (1 + x**2))**(-1)

    def _mass(self, r):
        """Compute the mass out to given radii analytically
        """
        x = r / self.rs
        return np.pi * self.rhos * (np.log(x**2 + 1) + 2 * np.log(x + 1) - 2 * np.arctan(x))


class NFWProfile(DensityProfile):
    """Navarro, Frenk, and White, ApJ 462, 563 (1996)
    http://arxiv.org/abs/astro-ph/9508025
    rho(r) = rhos / ((r/rs) * (1+r/rs)**2)
    """

    # def set_jval(self,jval,rs,dist):
    #    rhos = np.sqrt(3./(4.*np.pi)*jval*dist**2/rs**3)
    #    self.rs = rs
    #    self.rhos = rhos

    def _mass(self, r):
        """Compute the mass out to given radii analytically
        """
        x = r / self.rs
        return 4 * np.pi * self.rhos * self.rs**3 * (np.log(1 + x) - x / (1 + x))

    def jvalue_fast(self, r=None):
        """Fast integrated J-factor computation
        """
        if r is None:
            r = self.rmax
        x = r / self.rs
        return (4 * np.pi / 3.) * self.rhos**2 * self.rs**3 * (1 - (1 + x)**-3)

    def _rho(self, r):
        """Internal function for sub-class to compute the mass at given radii
        """
        x = r / self.rs
        return self.rhos * x**-1 * (1 + x)**-2


class EinastoProfile(DensityProfile):
    """ Einasto profile
    Einasto Trudy Inst. Astrofiz. Alma-Ata 5, 87 (1965) (Russian) [Eqn. 4]
    http://adsabs.harvard.edu/abs/1965TrAlm...5...87E
    rho(r) = rhos*exp(-2*((r/rs)**alpha-1)/alpha)
    ==>

    """
    _params = odict(
        list(DensityProfile._params.items()) +
        [
            ('alpha', Parameter(default=0.17)),
        ])

    @property
    def deriv_params(self):
        """Return the list of paramters we can take derivatives w.r.t.
        """
        return ["rs", "rhos", "alpha"]

    def _mass(self, r):
        """Compute the mass out to given radii analytically
        """
        x = r / self.rs
        gamma = spfn.gamma(3. / self.alpha)
        gammainc = spfn.gammainc(3. * self.alpha**-1, (2. * self.alpha**(-1) * x**self.alpha))
        alphainv = self.alpha**-1

        return 4 * np.pi * self.rhos * self.rs**3 * alphainv * \
            np.exp(2. * alphainv) * \
            np.power(2. * alphainv, -3. * alphainv) * \
            gamma * gammainc

    def _rho(self, r):
        """Return the denisty at a given radius
        """
        x = r / self.rs
        return self.rhos * np.exp(-2. * self.alpha**-1 * (x**(self.alpha) - 1))


class GNFWProfile(DensityProfile):
    """ Generalized NFW Profile
    Strigari et al. ApJ 678, 614 (2008) [Eqn. 3]
    http://arxiv.org/abs/0709.1510
    rho(r) = rhos / ( (r/rs)**gamma * (1+r/rs)**(3-gamma))
    """
    _params = odict(
        list(DensityProfile._params.items()) +
        [
            ('gamma', Parameter(default=1.)),
        ])

    @property
    def deriv_params(self):
        """Return the list of paramters we can take derivatives w.r.t.
        """
        return ["rs", "rhos", "gamma"]

    def _rho(self, r):
        """Internal function for sub-class to return the density at given radii
        """
        x = r / self.rs
        return self.rhos * x**(-self.gamma) * (1 + x)**(self.gamma - 3)

    def _mass(self, r):
        """Compute the mass out to given radii analytically
        """
        raise NotImplementedError("No analytic function for mass of GNFWProfile")


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
        list(DensityProfile._params.items()) +
        [
            ('alpha', Parameter(default=1.)),
            ('beta', Parameter(default=3.)),
            ('gamma', Parameter(default=1.)),
        ])

    @property
    def deriv_params(self):
        """Return the list of paramters we can take derivatives w.r.t.
        """
        return ["rs", "rhos", "alpha", "beta", "gamma"]

    def _rho(self, r):
        """Internal function for sub-class to return the density at given radii
        """
        x = r / self.rs
        return self.rhos * x**-self.gamma * \
            (1 + x**(1 / self.alpha))**(-(self.beta - self.gamma) * self.alpha)

    def _mass(self, r):
        """Compute the mass out to given radii analytically
        """
        raise NotImplementedError("No analytic function for mass of ZhouProfile")

Uniform = UniformProfile
Isothermal = IsothermalProfile
Burkert = BurkertProfile
NFW = NFWProfile
Einasto = EinastoProfile
gNFW = GNFWProfile
Zhou = ZhouProfile


def scale_list(l, scale_value):
    """Scale all the parameters on a list by the same value
    The operates on the list in place

    Parameters
    ----------

    l : list
    scale_value : float


    Returns
    -------

    l : list

    """
    for i, v in enumerate(l):
        l[i] = v * scale_value
    return l


def scale_dict(d, scale_value):
    """Scale all the parameters on in a dict by the same value
    The operates on the dictionary in place

    Parameters
    ----------

    d : dict
    scale_value : float


    Returns
    -------

    d : dict
        Dictionary with scaled values.

    """
    for k, v in d.items():
        if isinstance(v, list):
            d[k] = scale_list(v, scale_value)
        else:
            d[k] = v * scale_value
    return d


def scale_param(p, scale_value):
    """Generic function to scale parameters.
    Works for lists and dicts as well as simple paramters
    This operates in place.

    Parameters
    ----------

    p : dict or list or `Parameter`
    scale_value : float

    Returns
    -------

    p : dict or list or `Parameter`

    """
    if isinstance(p, dict):
        return scale_dict(p, scale_value)
    elif isinstance(p, list):
        return scale_list(p, scale_value)
    return p * scale_value


def scale_dict_param(d, k, scale_value, default_value):
    """Scale a parameter in a dict, or assign a default value

    Parameters
    ----------

    d : dict
        Input dictrionary

    k : str
        Key of `Parameter` to modify

    scale_value : float
        Value to scale existing `Parameter` by

    default_value : float
        Value to assign if k is not in d

    """
    try:
        d[k] = scale_param(d[k], scale_value)
    except KeyError:
        d[k] = scale_value * default_value


def factory(ptype, **kwargs):
    """Factory method to build `DenityProfile` objects
    Keyword arguments are passed to class c'tor

    Parameters
    ----------

    ptype : str
        Density profile type

    Returns
    -------

    profile : `DensityProfile`
        Newly created object

    """
    import dmsky.factory

    prof_copy = kwargs.copy()
    units = prof_copy.pop('units', None)
    if units:
        density, distance = units.rsplit('_', 1)
        scale_density = getattr(Units, density)
        scale_distance = getattr(Units, distance)

        scale_dict_param(prof_copy, 'rhos', scale_density, DensityProfile._params['rhos'].default)
        scale_dict_param(prof_copy, 'rs', scale_distance, DensityProfile._params['rs'].default)
        scale_dict_param(prof_copy, 'rmin', scale_distance, DensityProfile._params['rmin'].default)
        scale_dict_param(prof_copy, 'rmax', scale_distance, DensityProfile._params['rmax'].default)
        scale_dict_param(prof_copy, 'rhomax', scale_density,
                         DensityProfile._params['rhomax'].default)

    return dmsky.factory.factory(ptype, module=__name__, **prof_copy)
