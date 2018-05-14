#!/usr/bin/env python
"""
Python module for computing the line-of-sight integral over a
spherically symmetric distribution.

"""
__author__ = "Matthew Wood"
__date__ = "12/01/2011"

import numpy as np

from scipy.integrate import quad
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import UnivariateSpline

from dmsky.utils.units import Units


class LoSFn(object):
    """Integrand function (luminosity density) for LoS integration.  The
    parameter `alpha` introduces a change of variables:

    x' = x^(1/alpha).

    A value of alpha > 1 samples the integrand closer to x = 0
    (distance of closest approach). The change of variables requires
    the substitution:

    dx = alpha * (x')^(alpha-1) dx'

    """

    def __init__(self, dp, d, xi, alpha=3.0):
        """C'tor

        Parameters
        ----------

        dp : `DensityProfile`
            The density profile to integrate

        d : float
            Distance to the halo center

        xi: float
            Offset angle in radians.

        alpha: float
            Rescaling exponent for line-of-sight coordinate.
        """
        self.dp = dp
        self.d2 = d**2
        self.sinxi2 = np.sin(xi)**2
        self.alpha = alpha

    def __call__(self, xp, alpha=None):
        """"Return the density along the Line-of-sight.

        Parameters
        ----------

        xp : `numpy.array`
            Distance along the line of sight

        alpha : float
            Rescaling exponent for line-of-sight coordinate

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input xp

        """
        if alpha is None:
            alpha = self.alpha
        x = xp**alpha
        r = np.sqrt(x**2 + self.d2 * self.sinxi2)
        return self.func(r) * alpha * xp**(alpha - 1.0)

    def func(self, r):
        """Function to compute the integrand
        """
        return self.dp.rho(r)**2


class LoSAnnihilate(LoSFn):
    """Integrand function for LoS annihilation (J-factor)."""

    def func(self, r):
        """Function to compute the line-of-sight integrand
        """
        return self.dp.rho(r)**2


class LoSAnnihilate_Deriv(LoSFn):
    """Integrand function for LoS annihilation (J-factor)."""

    def __init__(self, dp, d, xi, paramNames, alpha=3.0):
        """C'tor

        Parameters
        ----------

        dp : `DensityProfile`
            The density profile to integrate

        d : float
            Distance to the halo center

        xi: float
            Offset angle in radians.

        paramNames : list
            Parameters to differentiate w.r.t.

        alpha: float
            Rescaling exponent for line-of-sight coordinate.

        """
        super(LoSAnnihilate_Deriv, self).__init__(dp, d, xi, alpha)
        self.__paramNames = paramNames

    def func(self, r):
        """Function to compute the line-of-sight integrand
        """
        return 2. * self.dp.rho(r) * self.dp.rho_deriv(r, self.__paramNames)


class LoSDecay(LoSFn):
    """Integrand function for LoS decay (D-factor)."""

    def func(self, r):
        """Function to compute the line-of-sight integrand
        """
        return self.dp.rho(r)


class LoSDecay_Deriv(LoSFn):
    """Integrand function for LoS decay (D-factor)."""

    def __init__(self, dp, d, xi, paramNames, alpha=1.0):
        """C'tor

        Parameters
        ----------

        dp : `DensityProfile`
            The density profile to integrate

        d : float
            Distance to the halo center

        xi: float
            Offset angle in radians.

        paramNames : list
            Parameters to differentiate w.r.t.

        alpha: float
            Rescaling exponent for line-of-sight coordinate.

        """
        super(LoSDecay_Deriv, self).__init__(dp, d, xi, alpha)
        self.__paramNames = paramNames

    def func(self, r):
        """Function to compute the line-of-sight integrand
        """
        return self.dp.rho(r) * self.dp.rho_deriv(r, self.__paramNames)


class LoSIntegral(object):
    """Slowest (and most accurate?) LoS integral. Uses
    scipy.integrate.quad with a change of variables to better sample
    the LoS close to the halo center.
    """

    def __init__(self, density, dhalo, alpha=3.0, ann=True, derivPar=None):
        """C'tor

        Parameters
        ----------

        density : `DensityProfile`
            The density profile to integrate

        dhalo : float
            Distance to the halo center

        alpha: float
            Rescaling exponent for line-of-sight coordinate.

        ann: bool
            True for annihilation, False for decay

        derivPar : list
            Parameters to differentiate w.r.t.

        """
        self.dp = density
        self.dhalo = dhalo
        self.alpha = float(alpha)
        self.ann = ann
        self.derivPar = derivPar

    @property
    def rmax(self):
        """Return the maximum integration radius
        """
        return self.dp.rmax

    @property
    def name(self):
        """Return the name of this profile
        """
        return self.__class__.__name__

    def _make_losfn(self, dhalo, psi):
        """Build a function that returns the line-of-sight integral

        Parameters
        ----------

        dhalo : float
            Distance to the halo center

        psi : float
            Angular offsets from the halo center


        Returns
        -------

        integrator : `LoSFn`
            Object that computes the line of sight integral

        """
        if self.ann:
            if self.derivPar is None:
                losfn = LoSAnnihilate(self.dp, dhalo, psi, self.alpha)
            else:
                losfn = LoSAnnihilate_Deriv(self.dp, dhalo, psi, [self.derivPar], self.alpha)
        else:
            if self.derivPar is None:
                losfn = LoSDecay(self.dp, dhalo, psi, self.alpha)
            else:
                losfn = LoSDecay_Deriv(self.dp, dhalo, psi, [self.derivPar], self.alpha)
        return losfn

    def __call__(self, psi, dhalo=None, degrees=False):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------

        psi     : `numpy.ndarray`
            Array of offset angles (in radians by default)

        dhalo   : `numpy.ndarray`
            Array of halo distances

        degrees : bool
            Flag to interpret `psi` in degrees instead of radians

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input psi

        """
        scalar = np.isscalar(psi)
        if degrees:
            psi = np.deg2rad(psi)
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

        v = self._integrate(psi, dhalo)

        if scalar:
            return v[0]
        return v

    def _integrate(self, psi, dhalo):
        """Evaluate the LoS integral at the offset angle psi for a halo
        located at the distance dhalo.

        Parameters
        ----------

        psi     : `numpy.ndarray`
            Array of offset angles (in radians by default)

        dhalo   : `numpy.ndarray`
            Array of halo distances

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input psi

        """

        # Arrays must be the same shape
        assert psi.shape == dhalo.shape

        # Output Array
        v = np.zeros_like(psi)

        # Closest approach to halo center
        rmin = dhalo * np.sin(psi)
        # Maximum extent of halo
        rmax = self.rmax

        for i, t in np.ndenumerate(psi):
            s0, s1 = 0, 0

            losfn = self._make_losfn(dhalo[i], psi[i])

            # Closest approach to halo center
            #rmin = dhalo[i]*np.sin(psi[i])

            # If observer inside the halo...
            if rmax > dhalo[i]:

                if psi[i] < np.pi / 2.:
                    x0 = np.power(dhalo[i] * np.cos(psi[i]), 1. / self.alpha)
                    s0 = 2 * quad(losfn, 0.0, x0)[0]

                    x1 = np.power(np.sqrt(rmax**2 - rmin[i]**2), 1. / self.alpha)
                    s1 = quad(losfn, x0, x1)[0]
                else:
                    x0 = np.power(np.abs(dhalo[i] * np.cos(psi[i])), 1. / self.alpha)

                    x1 = np.power(np.sqrt(rmax**2 - rmin[i]**2), 1. / self.alpha)
                    s1 = quad(losfn, x0, x1)[0]

            # If observer outside the halo...
            elif (rmax > rmin[i]) and (psi[i] < np.pi / 2.):
                x0 = np.power(np.sqrt(rmax**2 - rmin[i]**2), 1. / self.alpha)
                s0 = 2 * quad(losfn, 0.0, x0)[0]

            v[i] = s0 + s1

        return v

    def angularIntegral(self, angle=None):
        """Compute the solid-angle integrated j-value
        within a given radius

        Parameters
        ----------
        angle : `numpy.ndarray` or None
            Maximum integration angle (in degrees)
            If None, use the 'rmax' and 'dhalo' parameters.

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input xp

        """
        if angle is None:
            angle = np.degrees(np.arctan2(self.rmax, self.dhalo))

        angle = np.asarray(angle)
        if angle.ndim == 0:
            angle = np.array([angle])

        integrand = lambda r: self(r) * 2 * np.pi * np.sin(r)

        integral = []
        for a in angle:
            integral.append(quad(integrand, 1e-7, np.radians(a), full_output=True)[0])
        integral = np.asarray(integral)

        return integral


class LoSIntegralFast(LoSIntegral):
    """Vectorized version of LoSIntegral that performs midpoint
    integration with a fixed number of steps.
    """

    def __init__(self, density, dhalo, alpha=3.0, ann=True, nsteps=400, derivPar=None):
        """C'tor

        Parameters
        ----------

        density : `DensityProfile`
            The density profile to integrate

        dhalo : float
            Distance to the halo center

        alpha: float
            Rescaling exponent for line-of-sight coordinate.

        ann: bool
            True for annihilation, False for decay

        nsteps : int
            Number of steps for

        derivPar : list
            Parameters to differentiate w.r.t.
        """
        super(LoSIntegralFast, self).__init__(density, dhalo, alpha, ann, derivPar)

        self.nsteps = nsteps
        xedge = np.linspace(0, 1.0, self.nsteps + 1)
        self.x = 0.5 * (xedge[1:] + xedge[:-1])

    @property
    def rmax(self):
        """Return the maximum integration radius
        """
        if self.dp.rmax < np.inf:
            return self.dp.rmax
        return 1000 * self.dp.rs

    def _integrate(self, psi, dhalo):
        """Internal function to do the LOS integral

        Parameters
        ----------

        psi     : `numpy.ndarray`
            Array of offset angles (in radians by default)

        dhalo   : `numpy.ndarray`
            Array of halo distances

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input psi

        """
        # Arrays must be the same shape
        assert psi.shape == dhalo.shape

        v = np.zeros(psi.shape + (1,))

        dhalo = np.ones(v.shape) * dhalo[..., np.newaxis]
        psi = np.ones(v.shape) * psi[..., np.newaxis]

        # Closest approach to halo center
        rmin = dhalo * np.sin(psi)

        losfn = self._make_losfn(dhalo, psi)

        msk0 = self.rmax > dhalo
        msk1 = self.rmax > rmin

        # Distance between observer and point of closest approach
        xlim0 = np.abs(dhalo * np.cos(psi))**(1. / self.alpha)

        # Distance from point of closest approach to maximum
        # integration radius
        xlim1 = np.zeros(shape=psi.shape)
        xlim1[msk1] = np.sqrt(self.rmax**2 - rmin[msk1]**2)**(1. / self.alpha)

        # If observer inside the halo...
        if np.any(msk0):

            msk01 = msk0 & (psi < np.pi / 2.)
            msk02 = msk0 & ~(psi < np.pi / 2.)

            if np.any(msk01):
                dx0 = xlim0 / float(self.nsteps)
                dx1 = (xlim1 - xlim0) / float(self.nsteps)

                x0 = self.x * xlim0
                x1 = xlim0 + self.x * (xlim1 - xlim0)

                s0 = 2 * np.apply_over_axes(np.sum, losfn(x0) * dx0, axes=[-1])
                s1 = np.apply_over_axes(np.sum, losfn(x1) * dx1, axes=[-1])

                v[msk01] = s0[msk01] + s1[msk01]

            if np.any(msk02):
                dx1 = (xlim1 - xlim0) / float(self.nsteps)

                x1 = xlim0 + self.x * (xlim1 - xlim0)

                s0 = np.apply_over_axes(np.sum, losfn(x1) * dx1, axes=[-1])

                v[msk02] = s0[msk02]

        # Observer outside the halo...
        # Only calculate integral for psi < pi/2
        msk11 = (~msk0 & msk1) & (psi < np.pi / 2.)
        if np.any(msk11):
            dx0 = xlim1 / float(self.nsteps)
            x0 = xlim1 * self.x
            s0 = 2 * np.apply_over_axes(np.sum, losfn(x0) * dx0, axes=[-1])

            v[msk11] = s0[msk11]

        return v.reshape(v.shape[:-1])


class LoSIntegralInterp(LoSIntegralFast):
    """ Interpolate fast integral a for even faster look-up. """

    def __init__(self, density, dhalo, alpha=3.0, ann=True, nsteps=400, derivPar=None):
        """C'tor

        Parameters
        ----------

        density : `DensityProfile`
            The density profile to integrate

        dhalo : float
            Distance to the halo center

        alpha: float
            Rescaling exponent for line-of-sight coordinate.

        ann: bool
            True for annihilation, False for decay

        nsteps : int
            Number of steps for vectorization

        derivPar : list
            Parameters to differentiate w.r.t.

        """
        super(LoSIntegralInterp, self).__init__(density, dhalo, alpha, ann, nsteps, derivPar)
        self.func = self.create_func(self.dhalo)

    def create_profile(self, dhalo, nsteps=None):
        """Create a spatial J-factor profile

        Parameters
        ----------

        dhalo   : `numpy.ndarray`
            Array of halo distances

        nsteps : int
            Number of steps for vectorization


        Returns
        -------

        dhalo, psi : `numpy.meshgrid`
            Array of halo distances and angular offsets

        jval : `numpy.array`
            Corresponding J-factors

        """
        if not nsteps:
            nsteps = self.nsteps
        dhalo = np.unique(np.atleast_1d(dhalo))
        psi = np.logspace(np.log10(1e-7), np.log10(np.pi), nsteps)

        _dhalo, _psi = np.meshgrid(dhalo, psi)
        _jval = super(LoSIntegralInterp, self)._integrate(_psi, _dhalo)
        return np.log10([_dhalo, _psi, _jval])

    def create_func(self, dhalo):
        """Create the spline function

        Parameters
        ----------

        dhalo   : `numpy.ndarray`
            Array of halo distances

        Returns
        -------

        func : function
            A function that return J-factor as a function of psi and dhalo

        """
        log_dhalo, log_psi, log_jval = self.create_profile(dhalo)

        zeroval = -99
        log_jval[np.where(log_jval == -np.inf)] = zeroval

        if log_dhalo.shape[-1] == 1:
            # print('interp1d')
            # spline=UnivariateSpline(log_psi.flat,log_jval.flat,k=2,s=0)
            #fn = lambda psi: 10**(spline(np.log10(psi)))
            interp = interp1d(log_psi.flat, log_jval.flat, kind='linear')

            def fn(psi, dhalo):
                """Function to compute the J-factor
                """
                log_jval = interp(np.log10(psi))
                log_jval[np.where(log_jval < zeroval + 1)] = -np.inf
                return 10**log_jval
        else:
            # print('interp2d')
            #spline = bisplrep(log_psi,log_dhalo,log_jval,s=0.0,kx=2,ky=2)
            #fn = lambda psi,dhalo: 10**bisplev(np.log10(psi[:,0]),np.log10(dhalo[0,:]),spline)
            interp = interp2d(log_psi, log_dhalo, log_jval, kind='linear')

            def fn(psi, dhalo):
                """Function to compute the J-factor
                """
                log_jval = interp(np.log10(psi[:, 0]), np.log10(dhalo[0, :])).T
                log_jval[np.where(log_jval < zeroval + 1)] = -np.inf
                return 10**log_jval

        return fn

    def _integrate(self, psi, dhalo):
        """Internal function to the the LOS integral

        Parameters
        ----------

        psi     : `numpy.ndarray`
            Array of offset angles (in radians by default)

        dhalo   : `numpy.ndarray`
            Array of halo distances

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input psi

        """
        # Arrays must be the same shape
        if psi.shape != dhalo.shape:
            msg = "Shape of psi and dhalo must match"
            raise ValueError(msg)

        if psi.ndim > 1 and not (np.unique(psi) == psi[:, 0]).all():
            msg = "np.unique(psi) != psi[:,0]"
            raise ValueError(msg)

        if dhalo.ndim > 1 and not (np.unique(dhalo) == dhalo[0, :]).all():
            msg = "np.unique(dhalo) != dhalo[0,:]"
            raise ValueError(msg)

        # All halo distances within pre-existing interpolation range
        if ((np.max(self.dhalo) >= dhalo) & (np.min(self.dhalo) <= dhalo)).all():
            func = self.func
        else:
            func = self.create_func(dhalo)

        v = func(psi, dhalo)

        if v.shape != psi.shape:
            msg = "Input and output shape do not match"
            raise ValueError(msg)

        return v


class LoSIntegralFile(LoSIntegralInterp):
    """Interpolate over a pre-generated file.

    NOT IMPLEMENTED YET
    """

    def __init__(self, dp, dist, filename, ann=True):
        """C'tor

        Parameters
        ----------

        dp : `DensityProfile`
            The density profile to integrate

        dist : float
            Distance to the halo center

        filename: str
            File with tabulated results

        ann: bool
           True for annihilation, False for decay

        """
        super(LoSIntegralFile, self).__init__(dp, dist, ann=ann)
        self.filename = filename

    def create_profile(self, dhalo, nsteps=300):
        """Build the profile values"""
        log_psi, log_jval = np.loadtxt(self.filename, unpack=True)
        return self.dhalo, log_psi, log_jval


class ROIIntegrator(object):
    """Class to integrate a J-factor over a region of interest
    """
    def __init__(self, jspline, lat_cut, lon_cut, source_list=None):
        """ C'tor
        """
        self._jspline = jspline
        self._lat_cut = lat_cut
        self._lon_cut = lon_cut

        nbin_thetagc = 720
        thetagc_max = 180.

        self._phi_edges = np.linspace(0., 360., 720 + 1)
        self._theta_edges = np.linspace(0., thetagc_max, nbin_thetagc + 1)

        self._sources = None

        if not source_list is None:
            source_list = np.loadtxt(opts.source_list, unpack=True, usecols=(1, 2))
            self._sources = Vector3D.createLatLon(np.radians(source_list[0]),
                                                  np.radians(source_list[1]))

        self.compute()

    def compute(self):
        """Integrate the ROI
        """
        yaxis = Vector3D(np.pi / 2. * np.array([0., 1., 0.]))

        costh_edges = np.cos(np.radians(self._theta_edges))
        costh_width = costh_edges[:-1] - costh_edges[1:]

        phi = 0.5 * (self._phi_edges[:-1] + self._phi_edges[1:])
        self._theta = 0.5 * (self._theta_edges[:-1] + self._theta_edges[1:])

        self._jv = []
        self._domega = []

        for i0, th in enumerate(self._theta):

            jtot = integrate(lambda t: self._jspline(8.5 * Units.kpc, t) * np.sin(t),
                             np.radians(self._theta_edges[i0]),
                             np.radians(self._theta_edges[i0 + 1]), 100)

#    jval = jspline(np.radians(th))*costh_width[i0]
            v = Vector3D.createThetaPhi(np.radians(th), np.radians(phi))
            v.rotate(yaxis)

            lat = np.degrees(v.lat())
            lon = np.degrees(v.phi())

            src_msk = len(lat) * [True]

            if not self._sources is None:

                for k in range(len(v.lat())):
                    p = Vector3D(v._x[:, k])

                    sep = np.degrees(p.separation(self._sources))
                    imin = np.argmin(sep)
                    minsep = sep[imin]

                    if minsep < 0.62:
                        src_msk[k] = False

            msk = ((np.abs(lat) >= self._lat_cut) |
                   ((np.abs(lat) <= self._lat_cut) & (np.abs(lon) < self._lon_cut)))

            msk &= src_msk
            dphi = 2. * np.pi * float(len(lat[msk])) / float(len(phi))
            jtot *= dphi
#            jsum += jtot
#            domegasum += costh_width[i0]*dphi

            self._jv.append(jtot)
            self._domega.append(costh_width[i0] * dphi)

        self._jv = np.array(self._jv)
        self._jv_cum = np.cumsum(self._jv)

        self._jv_cum_spline = UnivariateSpline(self._theta_edges[1:],
                                               self._jv_cum,
                                               s=0, k=1)

        self._domega = np.array(self._domega)
        self._domega_cum = np.cumsum(self._domega)

    def eval(self, rgc, decay=False):
        """Evaluate the J-factor
        """
        if decay:
            units0 = Units.gev_cm2
            units1 = (8.5 * Units.kpc * 0.4 * Units.gev_cm3)
        else:
            units0 = Units.gev2_cm5
            units1 = (8.5 * Units.kpc * np.power(0.4 * Units.gev_cm3, 2))

        rgc = [float(t) for t in rgc.split('/')]

        if len(rgc) == 1:
            jv = self._jv_cum_spline(rgc[0])
            domega = np.cos(np.radians(rgc[0])) * 2 * np.pi / Units.deg2
        else:
            jv = self._jv_cum_spline(rgc[1]) - self._jv_cum_spline(rgc[0])
            domega = -(np.cos(np.radians(rgc[1])) -
                       np.cos(np.radians(rgc[0]))) * 2 * np.pi / Units.deg2

        print('%20.6g %20.6g %20.6g %20.6g' % (jv,
                                               jv / units0,
                                               jv / units1, domega))

    def print_profile(self, decay=False):
        """Print the profile
        """
        if decay:
            units0 = Units.gev_cm2
            units1 = (8.5 * Units.kpc * 0.4 * Units.gev_cm3)
        else:
            units0 = Units.gev2_cm5
            units1 = (8.5 * Units.kpc * np.power(0.4 * Units.gev_cm3, 2))

        for i, th in enumerate(self._theta_edges[1:]):

            jv = self._jv_cum[i]

            print('%10.2f %20.6g %20.6g %20.6g %20.6g' % (th, jv,
                                                          jv / units0,
                                                          jv / units1,
                                                          self._domega_cum[i]))
