#!/usr/bin/env python
"""
Module for target objects. A target carries the physical information
about a target (ra, dec, distance, density profile, etc.) and
interfaces with the `jcalc` module to calculate the l.o.s. integral at
a given sky position. Classes inherit from the `Target` baseclass and
can be created with the `factory` function.
"""
import sys
from os.path import abspath, dirname, join
from collections import OrderedDict as odict

import numpy as np


from dmsky.jcalc import LoSIntegral, LoSIntegralFast, LoSIntegralInterp
from dmsky.utils import coords
from dmsky.utils.units import Units
from dmsky.utils.tools import update_dict, merge_dict, item_version
from dmsky.library import ObjectLibrary
from dmsky.density import factory as density_factory
from dmsky.density import DensityProfile
from dmsky.priors import PriorFunctor
from dmsky.priors import factory as prior_factory

from pymodeler.model import Model
from pymodeler.parameter import Property, Derived


class Target(Model):
    """A DM search target """
    _params = odict([('title', Property(dtype=str,
                                        required=True, help='Human-readable name')),
                     ('name', Property(dtype=str,
                                       required=True, help='Machine-readable name')),
                     ('abbr', Property(dtype=str,
                                       required=True, help='Title abbreviation')),
                     ('profile', Property(dtype=dict,
                                          required=True, help='Density Profile (see `jcalc`)')),
                     ('version', Property(dtype=str,
                                          default=None, help='Which version of this target?')),
                     ('ver_key', Property(dtype=str,
                                          default=None, help='Short key for the version?')),
                     ('nickname', Property(dtype=str,
                                           default=None, help='Do we need this?')),
                     ('altnames', Property(dtype=list,
                                           default=[], help='Alternative names')),
                     ('ra', Property(dtype=float, format='%.3f',
                                     default=0.0, unit='deg',
                                     help='Right Ascension')),
                     ('dec', Property(dtype=float, format='%.3f',
                                      default=0.0, unit='deg',
                                      help='Declination')),
                     ('distance', Property(dtype=float, format='%.1f',
                                           default=0.0, unit='kpc',
                                           help='Distance')),
                     ('dist_err', Property(dtype=float, format='%.1f',
                                           default=0.0, unit='kpc',
                                           help='Distance Uncertainty')),
                     ('major_axis', Property(dtype=float, format='%.3f',
                                             default=np.nan, unit='kpc',
                                             help='Major axis')),
                     ('ellipticity', Property(dtype=float, format='%.3f',
                                              default=np.nan, unit='kpc',
                                              help='Major axis')),
                     ('references', Property(dtype=list, default=[],
                                             help='Literature references')),
                     ('color', Property(dtype=str, default='k',
                                        help='Plotting color')),
                     ('mode', Property(dtype=str, default='interp',
                                       help='L.o.S. Integration mode')),
                     ('j_rad_file', Property(dtype=str, default=None,
                                             help='File with J factor radial profile')),
                     ('d_rad_file', Property(dtype=str, default=None,
                                             help='File with D factor radial profile')),
                     ('j_map_file', Property(dtype=str, default=None,
                                             help='File with J factor map')),
                     ('d_map_file', Property(dtype=str, default=None,
                                             help='File with D factor map')),
                     ('j_prior_def', Property(dtype=dict,
                                              help='Details of J-factor prior distribution')),
                     ('d_prior_def', Property(dtype=dict,
                                              help='Details of D-factor prior distribution')),
                     ('density', Derived(dtype=DensityProfile, help='Density profile object')),
                     ('proftype', Derived(dtype=str, help='Profile type (see `jcalc`)')),
                     ('prof_par', Derived(dtype=np.ndarray,
                                          help='Profile parameters')),
                     ('prof_err', Derived(dtype=np.ndarray,
                                          help='Profile uncertainties')),
                     ('glat', Derived(dtype=float, format='%.3f', unit='deg',
                                      help='Galactic Longitude')),
                     ('glon', Derived(dtype=float, format='%.3f', unit='deg',
                                      help='Galactic Latitude')),
                     ('rad_max', Derived(dtype=float, format='%.1f', unit='kpc',
                                         help='Maximum integration radius')),
                     ('psi_max', Derived(dtype=float, format='%.3f',
                                         unit='deg', help='Maximum integration angle')),
                     ('j_integ', Derived(dtype=float, format='%.2e', unit='GeV2 / cm5',
                                         help='Integrated J factor')),
                     ('j_sigma', Derived(dtype=float, format='%.2e', unit='GeV2 / cm5',
                                         help='Uncertainty on integ. J factor')),
                     ('d_integ', Derived(dtype=float, format='%.2e', unit='GeV / cm2',
                                         help='Integrated D factor')),
                     ('d_sigma', Derived(dtype=float, format='%.2e', unit='GeV / cm2',
                                         help='Uncertainty on integ. D factor')),
                     ('j_profile', Derived(dtype=LoSIntegral, help='J factor profile')),
                     ('j_derivs', Derived(dtype=dict, help='J factor profile derivatives')),
                     ('j_prior', Derived(dtype=PriorFunctor,
                                         help='J factor likelihood prior distribution')),
                     ('d_profile', Derived(dtype=LoSIntegral, help='D factor profile')),
                     ('d_derivs', Derived(dtype=dict, help='D factor profile derivatives')),
                     ('d_prior', Derived(dtype=PriorFunctor,
                                         help='D factor likelihood prior distribution'))])

    def _density(self):
        """Create the `DensityProfile` object for this target
        """
        prof_copy = self.profile.copy()
        ptype = prof_copy.pop('type', "NFW")
        return density_factory(ptype, **prof_copy)

    def _proftype(self):
        """Get the type of `DensityProfile`
        """
        return self.profile.get('type')

    def _prof_par(self):
        """Get the profile parameters"""
        return self.density.param_values()

    def _prof_err(self):
        """Get the errors on the profile parameters"""
        return self.density.param_errors()

    def _glat(self):
        """Return the galactic latitude"""
        return coords.cel2gal(self.ra, self.dec)[1]

    def _glon(self):
        """Return the galactic longitude"""
        return coords.cel2gal(self.ra, self.dec)[0]

    def _rad_max(self):
        """Return the maximum integration radius"""
        units = self.getp('rad_max').unit
        return Units.convert_to(self.density.rmax, units)

    def _psi_max(self):
        """Return the maximum integration angle"""
        rmax = self.rad_max
        dist_kpc = self.distance
        if rmax > dist_kpc:
            return 180.
        return np.degrees(np.arcsin(rmax / dist_kpc))

    def _j_integ(self):
        """Compute the integrated J-factor
        """
        jprof = self.j_profile
        units = self.getp('j_integ').unit
        return Units.convert_to(jprof.angularIntegral(self.psi_max)[0], units)

    def _j_sigma(self):
        """Compute the uncertaintiy on the integrated J-factor
        """
        jd = self.j_derivs
        den = self.density
        dv = np.matrix(np.zeros((len(den.deriv_params))))
        for i, pname in enumerate(den.deriv_params):
            dv[0, i] = jd[pname].angularIntegral(self.psi_max)[0]
        return np.sqrt((dv * self.density.covar * dv.T)[0, 0])

    def _d_integ(self):
        """Compute the integrated D-factor
        """
        dprof = self.d_profile
        units = self.getp('j_integ').unit
        return Units.convert_to(dprof.angularIntegral(self.psi_max)[0], units)

    def _d_sigma(self):
        """Compute the uncertaintiy on the integrated D-factor
        """
        dd = self.d_derivs
        den = self.density
        dv = np.matrix(np.zeros((len(den.deriv_params))))
        for i, pname in enumerate(den.deriv_params):
            dv[0, i] = dd[pname].angularIntegral(self.psi_max)[0]
        return np.sqrt((dv * self.density.covar * dv.T)[0, 0])

    def _density_integral(self, ann=True, derivPar=None):
        """Return a functor that calculates the LoS integral for various cases

        Parameters
        ----------

        ann : bool
            build the functor for annihilation (i.e., integrate density^2 instead of density)

        derivPar : str
            build the functor for the derivative w.r.t. this parameter

        Returns
        -------

        los : `LoSIntegral`
           Object that caculates line-of-sight integrals

        """
        if self.mode == 'interp':
            return LoSIntegralInterp(self.density, self.distance *
                                     Units.kpc, ann=True, derivPar=derivPar)
        elif self.mode == 'fast':
            return LoSIntegralFast(self.density, self.distance *
                                   Units.kpc, ann=True, derivPar=derivPar)
        return LoSIntegral(self.density, self.distance * Units.kpc, ann=ann, derivPar=derivPar)

    def _j_profile(self):
        """Return an object that compute the J-factor at any direction

        Returns
        -------

        los : `LoSIntegral`
           Object that caculates line-of-sight integrals for J-factors

        """
        return self._density_integral(ann=True, derivPar=None)

    def _j_derivs(self):
        """Return an dict of objects that compute
        derivatives of the J-factor w.r.t. the profile parameters
        at any direction
        """
        retDict = {}
        for pname in self.density.deriv_params:
            retDict[pname] = self._density_integral(ann=True, derivPar=pname)
        return retDict

    def _d_profile(self):
        """Return an object that compute the D-factor at any direction

        Returns
        -------

        los : `LoSIntegral`
           Object that caculates line-of-sight integrals for D-factors

        """
        return self._density_integral(ann=False, derivPar=None)

    def _d_derivs(self):
        """Return an dict of objects that compute
        derivatives of the D-factor w.r.t. the profile parameters
        at any direction
        """
        retDict = {}
        for pname in self.density.deriv_params:
            retDict[pname] = self._density_integral(ann=False, derivPar=pname)
        return retDict

    def _j_prior(self):
        """Return the Prior on the J-factor

        Returns
        -------

        prior : `PriorFunctor`
           Object that compute Prior value

        """
        prior_copy = self.j_like_def.copy()
        prior_type = prior_copy.pop('type', 'l')
        the_prior = prior_factory(prior_type, **prior_copy)
        return the_prior

    def _d_prior(self):
        """Return the Prior on the D-factor

        Returns
        -------

        prior : `PriorFunctor`
           Object that compute Prior value

        """
        prior_copy = self.d_like_def.copy()
        prior_type = prior_copy.pop('type', 'l')
        the_prior = prior_factory(prior_type, **prior_copy)
        return the_prior

    def __str__(self, indent=0):
        """String represetation,
        includes name, ra, dec, distance and density
        """
        ret = '{0:>{1}}'.format('', indent)
        ret += self.__class__.__name__
        for k in ['name', 'ra', 'dec', 'distance', 'density']:
            v = getattr(self, k)
            ret += '\n  %-15s: %s' % (k, v)
        return ret

    def jvalue(self, ra, dec):
        """Return the J-factor in any direction

        Parameters
        ----------
        ra : `numpy.ndarray`
            Right-accension (in degrees)

        dec : `numpy.ndarray`
            Declination (in degrees)

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input ra, dec

        """
        sep = coords.angsep(self.ra, self.dec, ra, dec)
        return self.j_profile(np.radians(sep))

    def jsigma(self, ra, dec):
        """Return the uncertainty of J in any direction

        Parameters
        ----------
        ra : `numpy.ndarray`
            Right-accension (in degrees)

        dec : `numpy.ndarray`
            Declination (in degrees)

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input ra, dec

        """
        raise RuntimeError('jsigma computation not implemented')

    def dvalue(self, ra, dec):
        """Return the D-factor in any direction

        Parameters
        ----------
        ra : `numpy.ndarray`
            Right-accension (in degrees)

        dec : `numpy.ndarray`
            Declination (in degrees)

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input ra, dec

        """
        sep = coords.angsep(self.ra, self.dec, ra, dec)
        return self.d_profile(np.radians(sep))

    def dsigma(self, ra, dec):
        """Return the uncertainty of J in any direction

        Parameters
        ----------
        ra : `numpy.ndarray`
            Right-accension (in degrees)

        dec : `numpy.ndarray`
            Declination (in degrees)

        Returns
        -------

        values : `numpy.array`
             Return values, same shape as the input ra, dec

        """
        raise RuntimeError('dsigma computation not implemented')

    def _create_map(self, func, npix=150, subsample=4, coordsys='CEL', projection='AIT'):
        """Create a spatial map of a function

        Parameters
        ----------

        func : function
            Must take ra, dec as inputs the function we are mapping,

        npix : int
            Number of pixels along one axis of output map

        subsample : int
            Number of subsamples to take per pixel

        coordsys : str
            Coordniate system: 'GAL' or 'CEL'

        projection : str
            Map projection

        Returns
        -------

        image : `numpy.ndarray`
            Image data

        pix : `numpy.ndarray`
            Pixel coordinatates

        wcs : `WCS.wcs`
            WCS object for map

        """
        from dmsky.utils.wcs import create_image_wcs, get_pixel_skydirs
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        skydir = SkyCoord(self.ra * u.deg, self.dec * u.deg)

        cdelt = 2 * self.psi_max / npix
        subnpix = npix * subsample
        subcdelt = cdelt / subsample

        subimage, wcs = create_image_wcs(skydir, subcdelt, subnpix, coordsys, projection)
        pix = get_pixel_skydirs(subimage.shape, wcs)
        subimage = func(pix.ra, pix.dec).reshape(subimage.shape)

        # Take the mean of the subsampled pixels
        if subsample > 1:
            image, wcs = create_image_wcs(skydir, cdelt, npix, coordsys, projection)
            pix = get_pixel_skydirs(image.shape, wcs)
            image = (subimage.reshape(npix, subsample, npix, subsample)).mean(axis=3).mean(axis=1)
        else:
            image = subimage

        return image, pix, wcs

    def create_jmap(self, npix=150, subsample=4, coordsys='CEL', projection='AIT'):
        """Create a J-factor map

        Parameters
        ----------

        npix : int
            Number of pixels along one axis of output map

        subsample : int
            Number of subsamples to take per pixel

        coordsys : str
            Coordniate system: 'GAL' or 'CEL'

        projection : str
            Map projection

        Returns
        -------

        image : `numpy.ndarray`
            Image data

        pix : `numpy.ndarray`
            Pixel coordinatates

        wcs : `WCS.wcs`
            WCS object for map

        """
        return self._create_map(self.jvalue, npix, subsample, coordsys, projection)

    def create_dmap(self, npix=150, subsample=4, coordsys='CEL', projection='AIT'):
        """Create a D-factor map

        Parameters
        ----------

        npix : int
            Number of pixels along one axis of output map

        subsample : int
            Number of subsamples to take per pixel

        coordsys : str
            Coordniate system: 'GAL' or 'CEL'

        projection : str
            Map projection

        Returns
        -------

        image : `numpy.ndarray`
            Image data

        pix : `numpy.ndarray`
            Pixel coordinatates

        wcs : `WCS.wcs`
            WCS object for map

        """
        return self._create_map(self.dvalue, npix, subsample, coordsys, projection)

    def write_j_rad_file(self, j_rad_file=None, npts=50, minpsi=1e-4):
        """Write a text file with the J-fractor radial profile

        Parameters
        ----------

        j_rad_file : str
            Filename to write to

        npts : int
            Number of angles to write

        minpsi : float
            Value for smallest angle to write

        """
        if j_rad_file is None:
            j_rad_file = self.j_rad_file
        self._write_radial_profile(j_rad_file, self.j_profile, npts, minpsi)

    def write_d_rad_file(self, d_rad_file=None, npts=50, minpsi=1e-4):
        """Write a text file with the D-fractor radial profile

        Parameters
        ----------

        d_rad_file : str
            Filename to write to

        npts : int
            Number of angles to write

        minpsi : float
            Value for smallest angle to write

        """
        if d_rad_file is None:
            d_rad_file = self.d_rad_file
        self._write_radial_profile(d_rad_file, self.d_profile, npts, minpsi)

    def _write_radial_profile(self, filepath, profile, npts=50, minpsi=1e-4):
        """ Write radial profile to a text file

        The Fermi-LAT science tools expects a file two columns:
        Column 1: Angular offset in degrees
        Column 2: Profile value

        Column 2 should be normalized so that the angular integral
        over the sphere is 1.
        """

        psi_vals = np.zeros((npts))
        prof_vals = np.zeros((npts))
        psi_vals[0:-1] = np.logspace(np.log10(minpsi), np.log10(self.psi_max), npts - 1)
        psi_vals[-1] = 180.
        prof_vals[0:-1] = profile(np.radians(psi_vals[0:-1]))

        x_vals = np.radians(psi_vals)
        y_vals = 2. * np.pi * x_vals * prof_vals

        # Trapezoid summation
        norm = ((x_vals[1:] - x_vals[0:-1]) * (y_vals[1:] + y_vals[0:-1]) / 2.).sum()

        prof_vals /= norm

        if filepath is None:
            fout = sys.stdout
        else:
            fout = open(filepath, 'w!')
        for psi, prof in zip(psi_vals, prof_vals):
            fout.write("%0.3e %0.3e\n" % (psi, prof))


    def write_jmap_wcs(self, filename, npix=150, clobber=False,
                       map_kwargs=None, file_kwargs=None):
        """Write the J-factor to a template map.
        """
        from dmsky.utils.wcs import create_image_hdu

        if map_kwargs is None:
            map_kwargs = {}

        if file_kwargs is None:
            file_kwargs = {}

        image, pix, wcs = self.create_jmap(npix=npix, **map_kwargs)

        # This assumes square pixels.
        norm = np.sum(image) * np.radians(wcs.wcs.cdelt[0])**2
        norm_comment = "[%s] Normalization factor." % (self.getp('j_integ').unit)
        try:
            normerr = self.j_sigma
            normerr_comment = "[%s] Normalization uncertainty." % (self.getp('j_sigma').unit)
        except RuntimeError:
            normerr = None

        # Create the HDU
        hdu = create_image_hdu(image / norm, wcs)
        hdu.header.set('NORM', value=norm, comment=norm_comment)
        if normerr is not None:
            hdu.header.set('NORMERR', value=normerr, comment=normerr_comment)

        self.setp('j_map_file', value=filename)
        return hdu.writeto(filename, clobber=clobber, **file_kwargs)

    #def write_jmap_hpx(self, filename):
    #    """Write the J-factor to a template map.
    #    """
    #    raise RuntimeError('write_jmap_hpx not implemented')

    write_jmap = write_jmap_wcs

    def write_dmap_wcs(self, filename, npix=150, clobber=False,
                       map_kwargs=None, file_kwargs=None):
        """Write the D-factor to a template map.
        """
        from dmsky.utils.wcs import create_image_hdu
        if map_kwargs is None:
            map_kwargs = {}

        if file_kwargs is None:
            file_kwargs = {}

        image, pix, wcs = self.create_dmap(npix=npix, **map_kwargs)

        # This assumes square pixels.
        norm = np.sum(image) * np.radians(wcs.wcs.cdelt[0])**2
        norm_comment = "[%s] Normalization factor." % (self.getp('d_integ').unit)
        try:
            normerr = self.d_sigma
            normerr_comment = "[%s] Normalization uncertainty." % (self.getp('d_sigma').unit)
        except RuntimeError:
            normerr = None

        # Create the HDU
        hdu = create_image_hdu(image / norm, wcs)
        hdu.header.set('NORM', value=norm, comment=norm_comment)
        if normerr is not None:
            hdu.header.set('NORMERR', value=normerr, comment=normerr_comment)

        self.setp('d_map_file', value=filename)
        return hdu.writeto(filename, clobber=clobber, **file_kwargs)

    #def write_dmap_hpx(self, filename):
    #    """Write the D-factor to a template map.
    #    """
    #    raise RuntimeError('write_dmap_hpx not implemented')

    write_dmap = write_dmap_wcs


class Galactic(Target):
    """Class to add specifics for Galactic DM targets"""
    pass


class Dwarf(Target):
    """Class to add specifics for Dwarf Galaxy DM targets"""
    pass


class Galaxy(Target):
    """Class to add specifics for Galaxy DM targets"""
    pass


class Cluster(Target):
    """Class to add specifics for Galaxy Cluster DM targets"""
    pass


class Isotropic(Target):
    """Class to add specifics for Isotropic DM targets"""
    pass


def factory(ttype, **kwargs):
    """Factory function to tuild a `Target` objects
    """
    import dmsky.factory
    return dmsky.factory.factory(ttype, module=__name__, **kwargs)


class TargetLibrary(ObjectLibrary):
    """A top-level object, keeping track of all the `Target` objects
    that we have created
    """

    _suffix = 'targets'

    _defaults = (
        ('path', join(dirname(abspath(__file__)), 'data', _suffix)),
    )

    def get_target_dict(self, name, version=None, **kwargs):
        """Step through the various levels of dependencies to get the
        full dictionary for a target.

        target: version -> ... -> target: default -> default: type
        """
        n, v = item_version(name)
        if version is not None and v is not None:
            msg = "Version specified twice: %s, %s" % (name, version)
            raise ValueError(msg)

        if v is not None:
            version = v
        if version is None:
            version = 'default'
        name = n

        # Start with the target:version requested
        ret = self.library[name][version]

        # Walk down the chain until we either return None or the
        # 'default' version
        ret['version'] = version
        while (version is not None) and (version != 'default'):
            version = ret.get('base', 'default')
            ret = merge_dict(self.library[name][version], ret)
        kwargs['name'] = name
        # And finally, overwrite with kwargs
        update_dict(ret, kwargs)
        return ret

    def create_target(self, name, version=None, **kwargs):
        """Create a `Target`

        Parameters
        ----------

        name : str
            A name for the `Target`

        version : str
            Key that species which set of parameters to used
            for this target

        Returns
        -------

        target : `Target`
            The newly created Target

        """
        kw = self.get_target_dict(name, version, **kwargs)
        ttype = kw.get('type')
        return factory(ttype, **kw)
