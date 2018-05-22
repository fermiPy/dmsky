# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Utilities for statistical operations
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import scipy.stats as stats


def norm(x, mu, sigma=1.0):
    """ Scipy norm function """
    return stats.norm(loc=mu, scale=sigma).pdf(x)


def ln_norm(x, mu, sigma=1.0):
    """ Natural log of scipy norm function truncated at zero """
    return np.log(stats.norm(loc=mu, scale=sigma).pdf(x))


def lognorm(x, mu, sigma=1.0):
    """ Log-normal function from scipy """
    return stats.lognorm(sigma, scale=mu).pdf(x)


def log10norm(x, mu, sigma=1.0):
    """ Scale scipy lognorm from natural log to base 10
    x     : input parameter
    mu    : mean of the underlying log10 gaussian
    sigma : variance of underlying log10 gaussian
    """
    return stats.lognorm(sigma * np.log(10), scale=mu).pdf(x)


def ln_log10norm(x, mu, sigma=1.0):
    """ Natural log of base 10 lognormal """
    return np.log(stats.lognorm(sigma * np.log(10), scale=mu).pdf(x))


def gauss(x, mu, sigma=1.0):
    """Gaussian """
    s2 = sigma * sigma
    return 1. / np.sqrt(2 * s2 * np.pi) * np.exp(-(x - mu) * (x - mu) / (2 * s2))


def lngauss(x, mu, sigma=1.0):
    """Natural log of a Gaussian"""
    s2 = sigma * sigma
    return -0.5 * np.log(2 * s2 * np.pi) - np.power(x - mu, 2) / (2 * s2)


def lgauss(x, mu, sigma=1.0, logpdf=False):
    """ Log10 normal distribution...

    Parameters
    ----------

    x     : `numpy.array` or list
        Parameter of interest for scanning the pdf

    mu    : float
        Peak of the lognormal distribution (mean of the underlying
        normal distribution is log10(mu)

    sigma : float
        Standard deviation of the underlying normal distribution

    logpdf : bool
        Define the PDF in log space

    Returns
    -------

    vals : `numpy.array`
        Output values, same shape as x

    """
    x = np.array(x, ndmin=1)

    lmu = np.log10(mu)
    s2 = sigma * sigma

    lx = np.zeros(x.shape)
    v = np.zeros(x.shape)

    lx[x > 0] = np.log10(x[x > 0])

    v = 1. / np.sqrt(2 * s2 * np.pi) * np.exp(-(lx - lmu)**2 / (2 * s2))

    if not logpdf:
        v /= (x * np.log(10.))

    v[x <= 0] = -np.inf

    return v


def lnlgauss(x, mu, sigma=1.0, logpdf=False):
    """Log-likelihood of the natural log of a Gaussian
    """
    x = np.array(x, ndmin=1)

    lmu = np.log10(mu)
    s2 = sigma * sigma

    lx = np.zeros(x.shape)
    v = np.zeros(x.shape)

    mask = x > 0
    inv_mask = np.invert(mask)

    lx[mask] = np.log10(x[mask])

    v = -0.5 * np.log(2 * s2 * np.pi) - np.power(lx - lmu, 2) / (2 * s2)
    if not logpdf:
        v -= 2.302585 * lx + np.log(np.log(10.))

    if inv_mask.any():
        v[inv_mask] = -np.inf

    return v
