#!/usr/bin/env python
"""
Tests for J-factor calculations:
dmsky.jcalc.py
"""

from dmsky.jcalc import *
from dmsky.density import NFWProfile
import pylab as plt

rhos = 1
rmax = 300
psi = np.logspace(-6,np.log10(np.pi),500)

rs_list    = [0.01,1.0,20.0] # 10 pc, 1kpc, 20kpc
dhalo_list = [10,100,1000] # 10 kpc, 100kpc, 1000kpc

print("Testing 1D profiles" )
for rs in rs_list:
    fig,ax = plt.subplots(2,1)
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[1].set_xlabel('Offset Angle (deg)')
    ax[1].set_ylabel('Fractional Deviation')

    for dhalo in dhalo_list:
        dp = NFWProfile(rs=rs,rhos=rhos,rmax=rmax)
        slow = LoSIntegral(dp,dhalo)
        fast = LoSIntegralFast(dp,dhalo)
        itrp = LoSIntegralInterp(dp,dhalo)
         
        slow_jval = slow(psi,dhalo)
        fast_jval = fast(psi,dhalo)
        itrp_jval = itrp(psi,dhalo)
        
        zero = slow_jval == 0
        fast_jval[zero] = np.nan
        itrp_jval[zero] = np.nan

        fast_vs_slow = (fast_jval-slow_jval)/slow_jval
        itrp_vs_slow = (itrp_jval-slow_jval)/slow_jval
        itrp_vs_fast = (itrp_jval-fast_jval)/fast_jval

        ax[0].plot(np.degrees(psi),slow_jval,label='slow')
        ax[0].plot(np.degrees(psi),fast_jval,label='fast')
        ax[0].plot(np.degrees(psi),itrp_jval,label='itrp')
     
        ax[1].plot(np.degrees(psi),fast_vs_slow,label='fast vs slow')
        ax[1].plot(np.degrees(psi),itrp_vs_slow,label='itrp vs slow')

    plt.sca(ax[0])
    plt.title('rs=%g, rhos=%g, rmax=%g'%(rs,rhos,rmax))
    plt.sca(ax[0]); plt.legend(fontsize=10,loc='lower left',ncol=3)
    plt.sca(ax[1]); plt.legend(fontsize=10,loc='lower left',ncol=3)
    ax[1].set_ylim(-0.005,0.005)

# Now test 2D integrals
dhalo = np.logspace(1,3,5)
_dhalo,_psi = np.meshgrid(dhalo,psi)
print("Testing 2D profiles...")
for rs in rs_list:
    dp = NFWProfile(rs=rs,rhos=rhos,rmax=rmax)
    slow = LoSIntegral(dp,dhalo)
    fast = LoSIntegralFast(dp,dhalo)
    itrp = LoSIntegralInterp(dp,dhalo)
         
    slow_jval = slow(_psi,_dhalo)
    fast_jval = fast(_psi,_dhalo)
    itrp_jval = itrp(_psi,_dhalo)

    fast_vs_slow = (fast_jval-slow_jval)/slow_jval
    itrp_vs_slow = (itrp_jval-slow_jval)/slow_jval

    fig,ax = plt.subplots(2,2)
    plt.suptitle('rs=%g, rhos=%g, rmax=%g'%(rs,rhos,rmax))
    kwargs = dict(aspect='auto',interpolation='none',
                  vmax=np.log10(np.max(slow_jval)),
                  vmin=np.log10(slow_jval[slow_jval>0].min()))

    ax[0][0].imshow(np.log10(slow_jval),**kwargs)
    ax[0][1].imshow(np.log10(fast_jval),**kwargs)
    ax[1][0].imshow(np.log10(itrp_jval),**kwargs)
    im = ax[1][1].imshow(itrp_vs_slow,**dict(kwargs,vmin=-0.01,vmax=0.01))
    plt.colorbar(im)

plt.ion()
plt.show()
