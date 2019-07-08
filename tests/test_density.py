#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Alex Drlica-Wagner"
import os,sys
from collections import OrderedDict as odict
import inspect

import numpy as np
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')
import pylab as plt

import dmsky.density

dp = dmsky.density.DensityProfile()

module = 'dmsky.density'
fn = lambda member: inspect.isclass(member) and member.__module__==module
classes = odict(inspect.getmembers(sys.modules[module], fn))

radius = np.linspace(0.001,5,1000)

# Plot the different profiles
def test_densities():
    fig,ax = plt.subplots()
    y_min = 1e-4
    for i,cls in enumerate(classes):
        if cls == 'DensityProfile': continue
        dp = dmsky.density.factory(cls)
        dp.set_rho_r(1,dp.rs)
        denp = dp(radius)
        y_min = min(y_min, denp.min())
        print(dp)
        # Small offset for visualization
        epsilon = i*0.2e-1
        #epsilon = 0
        plt.plot(radius, denp * (1+epsilon),'-',label=cls)
        print('')
    ax.set_ymin(
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlabel("Radius")
    plt.ylabel("Density")
    plt.legend()

test_densities()

def test_min_max():
    fig,ax = plt.subplots()

    dp = dmsky.density.factory('NFWProfile',rmin=0.1)
    print(dp)
    epsilon = 0*2e-2
    plt.plot(radius,dp(radius) * (1+epsilon),'-',label='rmin=%g'%dp.rmin)
     
    dp = dmsky.density.factory('NFWProfile',rmax=2)
    print(dp)
    epsilon = 1*2e-2
    plt.plot(radius,dp(radius) * (1+epsilon),'-',label='rmax=%g'%dp.rmax)
     
    dp = dmsky.density.factory('NFWProfile',rhomax=100)
    print(dp)
    epsilon = 2*2e-2
    plt.plot(radius,dp(radius) * (1+epsilon),'-',label='rhomax=%g'%dp.rhomax)
     
    plt.xlabel("Radius")
    plt.ylabel("Density")
    plt.legend()
     
    ax.set_yscale('log')
    ax.set_xscale('log')

test_min_max()

plt.ion()
plt.show()
