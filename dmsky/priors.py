# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Classes to manage priors
"""
from __future__ import absolute_import, division, print_function

import numpy as np
import scipy.stats as stats
from scipy.interpolate import interp1d
from scipy.integrate import quad

from dmsky.utils import stat_funcs
from dmsky.utils import tools


class PriorFunctor(object):
    """A functor class that wraps simple functions we use to
       make priors on parameters.
    """

    def __init__(self, funcname, scale=1.0):
        """C'tor

        Parameters
        ----------

        funcname : str
            Name for this function, used for bookeeping

        scale : float
            Scale factor applied to input values.

        """
        self._funcname = funcname
        self._scale = scale

    def __call__(self, x):
        """Return the Prior value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        raise NotImplementedError("PriorFunctor.__call__")

    def log_value(self, x):
        """Return the log of the function value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        return np.log(self.__call__(self._scale*x))

    def normalization(self):
        """Normalization, i.e. the integral of the function
           over the normalization_range.
        """
        norm_r = self._normalization_range()
        return quad(self, norm_r[0]*self.scale, norm_r[1]*self.scale)[0]

    def _normalization_range(self):
        """Normalization range.
        """
        return 0, np.inf

    def mean(self):
        """Mean value of the function.
        """
        raise NotImplementedError("prior_functor.mean")

    def sigma(self):
        """The 'width' of the function.
        What this means depend on the function being used.
        """
        raise NotImplementedError("prior_functor.sigma")

    @property
    def funcname(self):
        """A string identifying the function.
        """
        return self._funcname

    @property
    def scale(self):
        """The scale factor applied to input values
        """
        return self._scale

    def marginalization_bins(self):
        """Binning to use to do the marginalization integrals

        Default is to marginalize over two decades,
        centered on mean, using 1000 bins
        """
        log_mean = np.log10(self.mean())
        return np.logspace(-1. + log_mean, 1. + log_mean, 1001)/self._scale

    def profile_bins(self):
        """The binning to use to do the profile fitting

        Default is to profile over +-5 sigma,
        Centered on mean, using 100 bins
        """
        log_mean = np.log10(self.mean())
        log_half_width = max(5. * self.sigma(), 3.)
        return np.logspace(log_mean - log_half_width,
                           log_mean + log_half_width, 101)/self._scale



class FunctionPrior(PriorFunctor):
    """Implementation of a prior that simply wraps an existing function
    """

    def __init__(self, funcname, mu, sigma, fn, lnfn=None, scale=1.0):
        """C'tor

        Parameters
        ----------

        funcname : str
            Name for this function, used for bookeeping

        mu : float
            Central value of the Prior, used to get scan ranges

        sigma : float
             Width of the Prior, used to get scan ranges

        fn : function
             Function that returns the Prior value

        lnfn : function or None
             Optional function that returns the log of the Prior value

        scale : float
            Scale factor applied to input values.

        """
        # FIXME, why doesn't super(FunctionPrior, self) work here?
        PriorFunctor.__init__(self, funcname, scale)
        self._mu = mu
        self._sigma = sigma
        self._fn = fn
        self._lnfn = lnfn

    def normalization(self):
        """The normalization
        i.e., the intergral of the function over the normalization_range
        """
        norm_r = self._normalization_range()
        return quad(self, norm_r[0]*self.scale, norm_r[1]*self.scale)[0]

    def mean(self):
        """Return the mean value of the function.
        """
        return self._mu

    def sigma(self):
        """Return the 'width' of the function.
        What this means depend on the function being used.
        """
        return self._sigma

    def log_value(self, x):
        """"Return the log of the function value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        if self._lnfn is None:
            return np.log(self._fn(x*self.scale, self._mu, self._sigma))
        return self._lnfn(x*self.scale, self._mu, self._sigma)

    def __call__(self, x):
        """Return the Prior value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        return self._fn(x*self.scale, self._mu, self._sigma)


class GaussPrior(FunctionPrior):
    """Implemenation of a Prior that wraps a Gaussian
    """

    def __init__(self, mu, sigma, scale=1.0):
        """C'tor

        Parameters
        ----------

        mu : float
            Central value of the Gaussian

        sigma : float
             Sigma of the Gaussian

        scale : float
            Scale factor applied to input values.

        """
        super(GaussPrior, self).__init__("gauss", mu, sigma,
                                         fn=stat_funcs.gauss,
                                         lnfn=stat_funcs.lngauss,
                                         scale=scale)


class LGaussPrior(FunctionPrior):
    """Implemenation of a Prior that wraps a log Gaussian
    """
    def __init__(self, mu, sigma, scale=1.0):
        """C'tor

        Parameters
        ----------

        mu : float
            Central value of the log Gaussian

        sigma : float
             Sigma of the Gaussian

        """
        super(LGaussPrior, self).__init__("lgauss", mu, sigma,
                                          fn=stat_funcs.lgauss,
                                          lnfn=stat_funcs.lnlgauss,
                                          scale=scale)


class LGaussLikePrior(FunctionPrior):
    """Implemenation of a Prior that wraps the
    inverse of the log of a Gaussian (i.e., x and y axes are swapped)
    """
    def __init__(self, mu, sigma, scale=1.0):
        """C'tor

        Parameters
        ----------

        mu : float
            Central value of the underlying Gaussian

        sigma : float
             Sigma of the underlying Gaussian

        scale : float
            Scale factor applied to input values.

        """
        def fn(x, y, s):
            """Swap the axes of the lgauss function"""
            return stat_funcs.lgauss(y, x, s)
        def lnfn(x, y, s):
            """Swap the axes of the lnlgauss function"""
            return stat_funcs.lnlgauss(y, x, s)
        super(LGaussLikePrior, self).__init__("lgauss_like", mu, sigma,
                                              fn=fn, lnfn=lnfn, scale=scale)


class LGaussLogPrior(FunctionPrior):
    """Implemenation of a Prior that wraps the
    inverse of the log of a Gaussian (i.e., x and y axes are swapped)
    The prior is implemented in log-space.
    """

    def __init__(self, mu, sigma, scale=1.0):
        """C'tor

        Parameters
        ----------

        mu : float
            Central value of the underlying Gaussian

        sigma : float
             Sigma of the underlying Gaussian

        scale : float
            Scale factor applied to input values.

        """
        def fn(x, y, s):
            """Swap the axes of the lgauss function and work in log space"""
            return stat_funcs.lgauss(x, y, s, logpdf=True)

        def lnfn(x, y, s):
            """Swap the axes of the lnlgauss function and work in log space"""
            return stat_funcs.lnlgauss(x, y, s, logpdf=True)
        super(LGaussLogPrior, self).__init__("lgauss_log", mu, sigma,
                                             fn=fn, lnfn=lnfn, scale=scale)


class LognormPrior(PriorFunctor):
    """ A wrapper around the lognormal function.

    A note on the highly confusing scipy.stats.lognorm function...
    The three inputs to this function are:

    s : This is the variance of the underlying
    gaussian distribution

    scale = 1.0 : This is the mean of the linear-space
    lognormal distribution.
    The mean of the underlying normal distribution
    occurs at ln(scale)

    loc = 0 : This linearly shifts the distribution in x (DO NOT USE)


    The convention is different for numpy.random.lognormal

    mean : This is the mean of the underlying
    normal distribution (so mean = log(scale))

    sigma : This is the standard deviation of the
    underlying normal distribution (so sigma = s)


    For random sampling:
    numpy.random.lognormal(mean, sigma, size)

    mean : This is the mean of the underlying
    normal distribution (so mean = exp(scale))

    sigma : This is the standard deviation of the
    underlying normal distribution (so sigma = s)

    scipy.stats.lognorm.rvs(s, scale, loc, size)

    s : This is the standard deviation of the
    underlying normal distribution

    scale : This is the mean of the generated
    random sample scale = exp(mean)

    Remember, pdf in log space is
    plot( log(x), stats.lognorm(sigma,scale=exp(mean)).pdf(x)*x )

    Parameters
    ----------
    mu : float
        Mean value of the function
    sigma : float
        Variance of the underlying gaussian distribution
    """

    def __init__(self, mu, sigma, scale=1.0):
        """C'tor

        Parameters
        ----------

        mu : float
            Mean value of the function

        sigma : float
            Variance of the underlying gaussian distribution

        scale : float
            Scale factor applied to input values.

        """
        super(LognormPrior, self).__init__('lognorm')
        self._mu = mu
        self._sigma = sigma

    def normalization(self):
        """Normalization, i.e. the integral of the function
           over the normalization_range.
        """
        return 1.

    def mean(self):
        """Mean value of the function.
        """
        return self._mu

    def sigma(self):
        """ The 'width' of the function.
        What this means depend on the function being used.
        """
        return self._sigma

    def __call__(self, x):
        """Return the Prior value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        return stats.lognorm(self._sigma, scale=self._mu).pdf(x*self.scale)


class NormPrior(PriorFunctor):
    """ A wrapper around the normal function.

    Parameters
    ----------
    mu : float
        Mean value of the function

    sigma : float
        Variance of the underlying gaussian distribution
    """

    def __init__(self, mu, sigma, scale=1.0):
        """C'tor

        Parameters
        ----------

        mu : float
            Mean value of the function

        sigma : float
            Variance of the underlying gaussian distribution

        scale : float
            Scale factor applied to input values.

        """
        super(NormPrior, self).__init__('norm', scale)
        self._mu = mu
        self._sigma = sigma

    def normalization(self):
        """Normalization, i.e. the integral of the function
           over the normalization_range.
        """
        return 1.

    def mean(self):
        """Mean value of the function.
        """
        return self._mu

    def sigma(self):
        """ The 'width' of the function.
        What this means depend on the function being used.
        """
        return self._sigma

    def __call__(self, x):
        """Return the Prior value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        return stats.norm(loc=self._mu, scale=self._sigma).pdf(x*self.scale)


class FileFuncPrior(PriorFunctor):
    """A wrapper around the interpolated function.

    Parameters
    ----------

    filename : string
       File with the function parameters

    """
    def __init__(self, filename, scale=1.0):
        """C'tor

        Parameters
        ----------

        filename : string
            File with the function parameters

        """
        super(FileFuncPrior, self).__init__('file')
        self._filename = filename
        d = tools.yaml_load(self._filename)
        self._mu = d['mean']
        self._sigma = d['sigma']
        self._x = d['x']
        self._y = d['y']
        self._kind = d.get('kind', 'linear')
        self._fill_value = d.get('fill_value', 0)
        self._interpfunc = interp1d(self._x, self._y, kind=self._kind,
                                    bounds_error=False, fill_value=self._fill_value)

    def mean(self):
        """Mean value of the function.
        """
        return self._mu

    def sigma(self):
        """ The 'width' of the function.
        What this means depend on the function being used.
        """
        return self._sigma

    def __call__(self, x):
        """Return the Prior value

        Parameters
        ----------

        x :`numpy.ndarray`
            Input values

        Returns
        -------

        y : `numpy.ndarray`
            Output values, same shape as x

        """
        return self._interpfunc(x*self.scale)


def create_prior_functor(d):
    """Build a prior from a dictionary.

    Parameters
    ----------
    d     :  A dictionary, it must contain:
       d['functype'] : a recognized function type
                       and all of the required parameters for the
                       prior_functor of the desired type

    Returns
    ----------
    A sub-class of '~fermipy.stats_utils.prior_functor'

    Recognized types are:

    'lognorm'       : Scipy lognormal distribution
    'norm'          : Scipy normal distribution
    'gauss'         : Gaussian truncated at zero
    'lgauss'        : Gaussian in log-space
    'lgauss_like'   : Gaussian in log-space, with arguments reversed.
    'lgauss_logpdf' : ???
    """
    functype = d.get('functype', 'lgauss_like')
    mu = d['mu']
    sigma = d['sigma']
    scale = d.get('scale', 1.0)
    if functype == 'norm':
        return NormPrior(mu, sigma, scale)
    elif functype == 'lognorm':
        return LognormPrior(mu, sigma, scale)
    elif functype == 'gauss':
        return GaussPrior(mu, sigma, scale)
    elif functype == 'lgauss':
        return LGaussPrior(mu, sigma, scale)
    elif functype in ['lgauss_like', 'lgauss_lik']:
        return LGaussLikePrior(mu, sigma, scale)
    elif functype == 'lgauss_log':
        return LGaussLogPrior(mu, sigma, scale)
    elif functype == 'interp':
        return FileFuncPrior(d['filename'], scale)
    else:
        raise KeyError("Unrecognized prior_functor type %s" % functype)


def factory(ptype, **kwargs):
    """Factor method to create Priors
    Keyword arguments are passed to class c'tor

    Parameters
    ----------

    ptype : str
        Prior type

    Returns
    -------

    prior : `PriorFunctor`
        Newly created object

    """
    import dmsky.factory

    prior_copy = kwargs.copy()
    return dmsky.factory.factory(ptype, module=__name__, **prior_copy)
