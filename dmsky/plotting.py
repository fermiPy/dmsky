#!/usr/bin/env python
"""
Module for plotting stuff in dmsky.
"""
import sys
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from dmsky.utils.units import Units


def plot_density(name,targets,xlims=None,nstep=100):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if xlims is None:
        xmin = Units.convert_to(targets[0].density.rmin,'kpc')
        xmax = Units.convert_to(targets[0].density.rmax,'kpc')
    else:
        xmin = xlims[0]
        xmax = xlims[1]

    if np.isinf(xmax):
        xmax = 10.* Units.convert_to(targets[0].density.rs,'kpc')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlim((xmin,xmax))
    ax.set_yscale('log')
    ax.set_xlabel('Distance [kpc]')
    ax.set_ylabel(r'Density [GeV / cm$^3$]')
    
    xvals = np.linspace(xmin,xmax,nstep)
    for target in targets:
        density = target.density
        yvals = Units.convert_to(density(Units.convert_from(xvals,'kpc')),'gev_cm3')
        ax.plot(xvals,yvals,color=target.color,label=target.name)
    
    leg = ax.legend(loc="upper right",fontsize=10,ncol=2)

    return fig,ax,leg


def plot_j_profile(name,targets,xlims=None,nstep=100,ylims=None):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if xlims is None:
        xmin = 0.
        xmax = targets[0].psi_max
    else:
        xmin = xlims[0]
        xmax = xlims[1]        

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlim((xmin,xmax))
    if ylims is not None:
        ax.set_ylim((ylims[0],ylims[1]))

    ax.set_yscale('log')
    ax.set_xlabel('Angular Separation [deg]')
    ax.set_ylabel(r'$dJ/d\Omega$ [GeV$^2$ / cm$^5$ sr]')
    
    xvals = np.linspace(xmin,xmax,nstep)[1:]
    for target in targets:
        j_profile = target.j_profile
        yvals = Units.convert_to(j_profile(xvals,degrees=True),'gev2_cm5')
        ax.plot(xvals,yvals,color=target.color,label=target.name)
    
    leg = ax.legend(loc="upper right",fontsize=10,ncol=2)

    return fig,ax,leg



def plot_j_integ(name,targets,xlims=None,nstep=100,ylims=None):
    """
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if xlims is None:
        xmin = 0.
        xmax = targets[0].psi_max
    else:
        xmin = xlims[0]
        xmax = xlims[1]        

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlim((xmin,xmax))
    if ylims is not None:
        ax.set_ylim((ylims[0],ylims[1]))

    ax.set_yscale('log')
    ax.set_xlabel('Angular Separation [deg]')
    ax.set_ylabel(r'$J$ [GeV$^2$ / cm$^5$]')
    
    xvals = np.linspace(xmin,xmax,nstep)[1:]
    for target in targets:
        j_profile = target.j_profile
        yvals = Units.convert_to(j_profile.angularIntegral(xvals),'gev2_cm5')
        ax.plot(xvals,yvals,color=target.color,label=target.name)
    
    leg = ax.legend(loc="lower right",fontsize=10,ncol=2)

    return fig,ax,leg





if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    rl = RosterLibrary()
    ro = rl.create_roster('martinez2015_nfw')
  
    plot_dict = {}
    plot_dict['Density'] = plot_density('Density',tlist,xlims=(0.,1.))
    plot_dict['dJ_dOmega'] = plot_j_profile('dJ_dOmega',tlist,xlims=(0.,1.),ylims=(1e19,1e27))
    plot_dict['J_int'] = plot_j_integ('J',tlist,xlims=(0.,0.2),ylims=(1e16,1e20))
        
    FORMATS = ['png']

    for k,v in plot_dict.items():
        for fmt in FORMATS:
            outname = "%s.%s"%(k,fmt)
            print "Writing %s"%(outname)
            v[0].savefig(outname)

    print ""
  
    
