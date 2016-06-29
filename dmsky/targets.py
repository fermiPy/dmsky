#!/usr/bin/env python
"""
Module for target objects. A target carries the physical information
about a target (ra, dec, distance, density profile, etc.) and
interfaces with the `jcalc` module to calculate the l.o.s. integral at
a given sky position. Classes inherit from the `Target` baseclass and
can be created with the `factory` function.
"""
import sys
import os, os.path
from os.path import abspath, dirname, join

import numpy as np

from collections import OrderedDict as odict

from dmsky.jcalc import LoSIntegral, LoSIntegralFast, LoSIntegralInterp
from dmsky.utils import coords
from dmsky.utils.units import Units
from dmsky.utils.tools import update_dict, merge_dict, yaml_load, get_items, item_version
from dmsky.library import ObjectLibrary
from dmsky.density import factory as density_factory
from dmsky.density import DensityProfile


from pymodeler.model import Model
from pymodeler.parameter import *


class Target(Model):
    _params = odict([('title'      ,Property(dtype=str   ,required=True , comment='Human-readable name')),
                     ('name'       ,Property(dtype=str   ,required=True , comment='Machine-readable name')),
                     ('abbr'       ,Property(dtype=str   ,required=True , comment='Title abbreviation')),
                     ('profile'    ,Property(dtype=dict  ,required=True , comment='Density Profile (see `jcalc`)')),
                     ('version'    ,Property(dtype=str   ,default=None,   comment='Which version of this target?')),
                     ('nickname'   ,Property(dtype=str   ,default=None,   comment='Do we need this?')),
                     ('altnames'   ,Property(dtype=list  ,default=[],     comment='Alternative names')),
                     ('ra'         ,Property(dtype=float ,format='%.3f', default=0.0     ,unit='deg', 
                                             comment='Right Ascension')),
                     ('dec'        ,Property(dtype=float ,format='%.3f', default=0.0     ,unit='deg',
                                             comment='Declination')),
                     ('distance'   ,Property(dtype=float ,format='%.1f', default=0.0     ,unit='kpc',
                                             comment='Distance')),
                     ('dist_err'   ,Property(dtype=float ,format='%.1f', default=0.0     ,unit='kpc',
                                             comment='Distance Uncertainty')),
                     ('major_axis' ,Property(dtype=float ,format='%.3f', default=np.nan  ,unit='kpc',
                                             comment='Major axis')),
                     ('ellipticity',Property(dtype=float ,format='%.3f', default=np.nan  ,unit='kpc',
                                             comment='Major axis')),
                     ('references' ,Property(dtype=list  ,default=[],
                                             comment='Literature references')),
                     ('color'      ,Property(dtype=str   ,default='k',
                                             comment='Plotting color')),
                     ('mode'       ,Property(dtype=str   ,default='interp',
                                             comment='L.o.S. Integration mode')),
                     ('density'    ,Derived(dtype=DensityProfile ,comment='Density profile object')),
                     ('proftype'   ,Derived(dtype=str            ,comment='Profile type (see `jcalc`)')), 
                     ('prof_par'   ,Derived(dtype=np.ndarray     ,comment='Profile parameters')),
                     ('prof_err'   ,Derived(dtype=np.ndarray     ,comment='Profile uncertainties')),
                     ('glat'       ,Derived(dtype=float, format='%.3f', unit='deg', 
                                            comment='Galactic Longitude')),
                     ('glon'       ,Derived(dtype=float, format='%.3f', unit='deg',         
                                            comment='Galactic Latitude')),
                     ('rad_max'    ,Derived(dtype=float, format='%.1f', unit='kpc',         
                                            comment='Maximum integration radius')),
                     ('psi_max'    ,Derived(dtype=float, format='%.3f', unit='deg'
                                            ,comment='Maximum integration angle')),
                     ('j_integ'    ,Derived(dtype=float, format='%.2e', unit='GeV2 cm-5',
                                            comment='Integrated J factor')),
                     ('j_sigma'    ,Derived(dtype=float, format='%.2e', unit='GeV2 cm-5',
                                            comment='Uncertainty on integ. J factor')),
                     ('d_integ'    ,Derived(dtype=float, format='%.2e', unit='GeV cm-2',
                                            comment='Integrated D factor')),
                     ('d_sigma'    ,Derived(dtype=float, format='%.2e', unit='GeV cm-2',
                                            comment='Uncertainty on integ. D factor')),
                     ('j_profile'  ,Derived(dtype=LoSIntegral, comment='J factor profile')),
                     ('j_derivs'   ,Derived(dtype=dict,        comment='J factor profile derivatives')),
                     ('j_map_file' ,Derived(dtype=str,         comment='File with J factor map')),
                     ('d_profile'  ,Derived(dtype=LoSIntegral, comment='D factor profile')),
                     ('d_derivs'   ,Derived(dtype=dict,        comment='D factor profile derivatives')),
                     ('d_map_file' ,Derived(dtype=str,         comment='File with D factor map'))])                    


    def __init__(self, **kwargs):
        super(Target,self).__init__(**kwargs)    

    def _density(self):
        prof_copy = self.profile.copy()
        ptype = prof_copy.pop('type',"NFW")
        return density_factory(ptype,**prof_copy)

    def _proftype(self):
        return self.profile.get('type')

    def _prof_par(self):
        return self.density.param_values()
  
    def _prof_err(self):
        return self.density.param_errors()

    def _glat(self):
        return coords.cel2gal(self.ra,self.dec)[1]

    def _glon(self):
        return coords.cel2gal(self.ra,self.dec)[0]

    def _rad_max(self):
        return self.density.rmax

    def _psi_max(self):
        rmax = self.rad_max
        dist_kpc = self.distance * Units.kpc
        if rmax > dist_kpc:
            return 180.
        else:
            return np.degrees(np.arcsin(rmax/dist_kpc))
      
    def _j_integ(self):
        jprof = self.j_profile        
        return jprof.angularIntegral(self.psi_max)[0]
 
    def _j_sigma(self):
        jd = self.j_derivs
        den = self.density
        dv = np.matrix(np.zeros((len(den.deriv_params))))
        for i,pname in enumerate(den.deriv_params):
            dv[0,i] = jd[pname].angularIntegral(self.psi_max)[0]
        return np.sqrt((dv * self.density.covar * dv.T)[0,0])

    def _j_map_file(self):
        raise Exception('Not implemented')

    def _d_integ(self):
        dprof = self.d_profile
        return dprof.angularIntegral(self.psi_max)[0]
 
    def _d_sigma(self):
        dd = self.d_derivs
        den = self.density
        dv = np.matrix(np.zeros((len(den.deriv_params))))
        for i,pname in enumerate(den.deriv_params):
            dv[i] = dd[pname].angularIntegral(self.psi_max)[0]
        return np.sqrt((dv.T * self.density.covar * dv)[0,0])

    def _d_map_file(self):
        raise Exception('Not implemented')


    def _density_integral(self,ann=True,derivPar=None):
        """ return a functor that calculates the LoS integral for various cases

        ann      : build the functor for annihilation (i.e., integrate density^2 instead of density)
        derivPar : build the functor for the derivative w.r.t. this parameter 
        """
        if self.mode  == 'interp':
            return LoSIntegralInterp(self.density, self.distance*Units.kpc, ann=True, derivPar=derivPar)
        elif self.mode == 'fast':
            return LoSIntegralFast(self.density, self.distance*Units.kpc, ann=True, derivPar=derivPar)
        else:
            return LoSIntegral(self.density, self.distance*Units.kpc, ann=True, derivPar=derivPar)
   
    def _j_profile(self):
        return self._density_integral(ann=True,derivPar=None)

    def _j_derivs(self):
        retDict = {}
        for pname in self.density.deriv_params:            
            retDict[pname] = self._density_integral(ann=True,derivPar=pname)
        return retDict

    def _d_profile(self):
        return self._density_integral(ann=False,derivPar=None)

    def _d_derivs(self):
        retDict = {}
        for pname in self.density.deriv_params:
            retDict[pname] = self._density_integral(ann=False,derivPar=pname)
        return retDict
      
    def __str__(self):
        ret = self.__class__.__name__
        for k in ['name','ra','dec','distance','density']:
            v = getattr(self,k)
            ret += '\n  %-15s: %s'%(k,v)
        return ret

    def jvalue(self,ra,dec):
        sep = coords.angsep(self.ra,self.dec,ra,dec)
        return self.j_profile(np.radians(sep))

    def jsigma(self,ra,dec):
        raise Exception('Not implemented')

    def dvalue(self,ra,dec):
        sep = coords.angsep(self.ra,self.dec,ra,dec)        
        return self.d_profile(np.radians(sep))
    
    def dsigma(self,ra,dec):
        raise Exception('Not implemented')


class Galactic(Target): pass
class Dwarf(Target): pass
class Galaxy(Target): pass
class Cluster(Target): pass
class Isotropic(Target): pass

def factory(type, **kwargs):
    import dmsky.factory
    return dmsky.factory.factory(type, module=__name__,**kwargs)

class TargetLibrary(ObjectLibrary):
    _defaults = (
        ('path', join(dirname(abspath(__file__)),'data/targets')),
    )

    def get_target_dict(self, name, version=None, **kwargs):
        """ Step through the various levels of dependencies to get the
        full dictionary for a target.

        target: version -> ... -> target: default -> default: type
        """
        n,v = item_version(name)
        if version is not None and v is not None:
            msg = "Version specified twice: %s, %s"%(name,version)
            raise ValueError(msg)

        if v is not None:   version = v
        if version is None: version = 'default'
        name = n

        # Start with the target:version requested
        ret = self.library[name][version]

        # Walk down the chain until we either return None or the
        # 'default' version
        while (version is not None) and (version != 'default'):
            version = ret.get('base','default')
            ret = merge_dict(self.library[name][version], ret)
        ret['version'] = version
        kwargs['name'] = name    
        # And finally, overwrite with kwargs
        update_dict(ret,kwargs)
        return ret

    def create_target(self, name, version=None, **kwargs):
        kw = self.get_target_dict(name,version,**kwargs)
        return factory(**kw)

if __name__ == "__main__":
    import argparse
    description = __doc__
    parser = argparse.ArgumentParser(description=description)
    args = parser.parse_args()

    
    prof_dict = dict(type='NFW',
                     units='msun_kpc3_kpc',
                     rs=0.27,
                     rhos=2.64e+08,
                     rmax=1.0)
    targ = Target(distance=50.,profile=prof_dict)
