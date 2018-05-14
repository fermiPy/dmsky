#!/usr/bin/env python
"""
Module for plotting stuff in dmsky.
"""

import numpy as np
import matplotlib.pyplot as plt

from dmsky.utils.units import Units


def plot_density(name, targets, xlims=None, nstep=100):
    """Make a plot of the density as a function of distance from
    the target center.

    Parameters
    ----------

    name : str
        Name for the plot

    targets : list
        List of targets to include in the plot

    xlims : tuple
        Range for the x-axis of the plot

    nteps : int
        Number of points to include in the plot

    Returns
    -------

    fig : `matplotlib.Figure`

    ax : `matplotlib.Axes`

    leg : `matplotlib.Legend`

    """
    fig = plt.figure(name)
    ax = fig.add_subplot(111)

    if xlims is None:
        xmin = Units.convert_to(targets[0].density.rmin, 'kpc')
        xmax = Units.convert_to(targets[0].density.rmax, 'kpc')
    else:
        xmin = xlims[0]
        xmax = xlims[1]

    if np.isinf(xmax):
        xmax = 10. * Units.convert_to(targets[0].density.rs, 'kpc')

    fig = plt.figure(name)
    ax = fig.add_subplot(111)

    ax.set_xlim((xmin, xmax))
    ax.set_yscale('log')
    ax.set_xlabel('Distance [kpc]')
    ax.set_ylabel(r'Density [GeV / cm$^3$]')

    xvals = np.linspace(xmin, xmax, nstep)
    for target in targets:
        density = target.density
        yvals = Units.convert_to(density(Units.convert_from(xvals, 'kpc')), 'gev_cm3')
        ax.plot(xvals, yvals, color=target.color, label=target.name)

    leg = ax.legend(loc="upper right", fontsize=10, ncol=2)

    return fig, ax, leg


def plot_j_profile(name, targets, xlims=None, nstep=100, ylims=None):
    """Make a plot of the J-factor as a function of the angle from
    the target center.

    Parameters
    ----------

    name : str
        Name for the plot

    targets : list
        List of targets to include in the plot

    xlims : tuple
        Range for the x-axis of the plot

    nteps : int
        Number of points to include in the plot

    ylims : tuple
        Range for the y-axis of the plot

    Returns
    -------

    fig : `matplotlib.Figure`

    ax : `matplotlib.Axes`

    leg : `matplotlib.Legend`

    """
    fig = plt.figure(name)
    ax = fig.add_subplot(111)

    if xlims is None:
        xmin = 0.
        xmax = targets[0].psi_max
    else:
        xmin = xlims[0]
        xmax = xlims[1]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlim((xmin, xmax))
    if ylims is not None:
        ax.set_ylim((ylims[0], ylims[1]))

    ax.set_yscale('log')
    ax.set_xlabel('Angular Separation [deg]')
    ax.set_ylabel(r'$dJ/d\Omega$ [GeV$^2$ / cm$^5$ sr]')

    xvals = np.linspace(xmin, xmax, nstep)[1:]
    for target in targets:
        j_profile = target.j_profile
        yvals = Units.convert_to(j_profile(xvals, degrees=True), 'gev2_cm5')
        ax.plot(xvals, yvals, color=target.color, label=target.name)

    leg = ax.legend(loc="upper right", fontsize=10, ncol=2)

    return fig, ax, leg


def plot_j_integ(name, targets, xlims=None, nstep=100, ylims=None):
    """Make a plot of the integrated J-factor as a function of the angle from
    the target center.

    Parameters
    ----------

    name : str
        Name for the plot

    targets : list
        List of targets to include in the plot

    xlims : tuple
        Range for the x-axis of the plot

    nteps : int
        Number of points to include in the plot

    ylims : tuple
        Range for the y-axis of the plot

    Returns
    -------

    fig : `matplotlib.Figure`

    ax : `matplotlib.Axes`

    leg : `matplotlib.Legend`

    """
    fig = plt.figure(name)
    ax = fig.add_subplot(111)

    if xlims is None:
        xmin = 0.
        xmax = targets[0].psi_max
    else:
        xmin = xlims[0]
        xmax = xlims[1]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlim((xmin, xmax))
    if ylims is not None:
        ax.set_ylim((ylims[0], ylims[1]))

    ax.set_yscale('log')
    ax.set_xlabel('Angular Separation [deg]')
    ax.set_ylabel(r'$J$ [GeV$^2$ / cm$^5$]')

    xvals = np.linspace(xmin, xmax, nstep)[1:]
    for target in targets:
        j_profile = target.j_profile
        yvals = Units.convert_to(j_profile.angularIntegral(xvals), 'gev2_cm5')
        ax.plot(xvals, yvals, color=target.color, label=target.name)

    leg = ax.legend(loc="lower right", fontsize=10, ncol=2)

    return fig, ax, leg
