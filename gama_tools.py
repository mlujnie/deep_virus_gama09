import sys
import os
import os.path
import subprocess
import numpy as np

import tables as tb
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.visualization import ZScaleInterval

from astropy.table import Table, vstack
from astropy.io import ascii
#from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u 
from scipy.ndimage import gaussian_filter
from astropy.stats import biweight_location, biweight_scale
from scipy.optimize import curve_fit

import glob

#from hetdex_api.shot import *

def gaia_query(ra_deg, dec_deg, rad_deg, maxmag=20, 
               maxsources=10000): 
    """
    Query Gaia DR1 @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 
                             'phot_g_mean_mag'], 
                    column_filters={"phot_g_mean_mag": 
                                    ("<%f" % maxmag)}, 
                    row_limit = maxsources) 
 
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, 
                           unit=(u.deg, u.deg), 
                           frame='icrs')
    return vquery.query_region(field, 
                               radius=("%fd" % rad_deg), 
                               catalog="I/337/gaia")[0] 

def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=20,
                    maxsources=10000):
    """
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    """
    vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag',
                             'ymag', 'e_ymag'],
                    column_filters={"gmag":
                                    ("<%f" % maxmag)},
                    row_limit=maxsources)

    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               radius=("%fd" % rad_deg),
                               catalog="II/349/ps1")[0]


def get_table(gama_field, keys=["spectrum", "error"], clean=True, perc=90, test=False):
    
    """Open shot files and join tables to get one super-table."""
    
    config = GamaConfig()
    def_wave = config.wavelength
    
    
    N=len(config.shotlist[gama_field])
    all_tables = []
    print("Getting: "+", ".join(keys))
    if test:
        N = 2
    for idx in range(N):
        shotid = config.shotlist[gama_field][idx]
        fileh = tb.open_file(os.path.join(config.mxhf_base, "gama_recon/",shotid + ".h5"))
        t = Table(fileh.root.Info.read())
        fibers = fileh.root.Fibers
        t["shotid"] = [shotid for k in range(len(t))]
        for key in keys:
            t[key] = fibers[:][key]
          
        if not clean:
            all_tables.append(t)
            print("Done with "+shotid+f" ({idx+1}/{N})")
            continue
            
        mask = np.ones(len(t), dtype="bool")
        spec_here = t["spectrum"]

        # exclude extreme continuum values
        #perc = 93

        wlcont_lo = (def_wave > 4000)&(def_wave <= 4500)
        medians_lo = np.nanmedian(spec_here[:,wlcont_lo], axis=1)
        perc_lo = np.nanpercentile(medians_lo, perc)

        wlcont_hi = (def_wave > 4800)&(def_wave <= 5300)
        medians_hi = np.nanmedian(spec_here[:,wlcont_hi], axis=1)
        perc_hi = np.nanpercentile(medians_hi, perc)
        mask[abs(medians_lo)>perc_lo] = False
        mask[abs(medians_hi)>perc_hi] = False
        
        #mask_idx = np.where(~mask)[0]
        #for fiber_idx in mask_idx:
        #    mask[fiber_idx-1] = False
        #    mask[fiber_idx+1] = False

        t["mask"] = mask
        print("Masked {}/{} fibers.".format(len(mask[~mask]), len(mask)))
        print("Finished "+shotid)
            
        all_tables.append(t)
        print("Done with "+shotid+f" ({idx+1}/{N})")
    print("Joining tables...")
    jt = vstack(all_tables)
    fiber_coords = SkyCoord(ra=jt["ra"]*u.deg, dec=jt["dec"]*u.deg)
    
    GAIA = False
    PANSTARRS = False
    if GAIA:
        print("Masking GAIA sources...")
        gaia_sources = ascii.read(os.path.join(config.basedir, "gaia_sources.tab"))
        for gs in gaia_sources:
            source_coords = SkyCoord(ra=gs["RA_ICRS"]*u.deg, dec=gs["DE_ICRS"]*u.deg)
            fiber_dists = source_coords.separation(fiber_coords).arcsec
            jt["mask"][fiber_dists < 2] = False
    if PANSTARRS:
        print("Masking PANSTARRS sources...")
        gaia_sources = ascii.read(os.path.join(config.basedir, "panstarrs_sources.tab"))
        for gs in gaia_sources:
            source_coords = SkyCoord(ra=gs["RAJ2000"]*u.deg, dec=gs["DEJ2000"]*u.deg)
            fiber_dists = source_coords.separation(fiber_coords).arcsec
            jt["mask"][fiber_dists < 2] = False        
    print("Done.")
    return jt

def fit_moffat(dist, amp, fwhm):
    beta = 3.
    gamma = fwhm/(2*np.sqrt(2**(1/beta) - 1))
    norm = (beta-1)/(np.pi*gamma**2)
    return amp * norm * (1+(dist/gamma)**2)**(-1*beta)

def integrate_moffat(dist, amp, fwhm):
    """integrates the moffat function over the fiber area"""

    INTSTEP = 0.1

    dist_xy = dist/np.sqrt(2)
    gridrange = np.arange(dist_xy-0.75, dist_xy+0.75+INTSTEP, INTSTEP) # diameter of a fiber is 1.5'' -> radius = 0.75''
    xgrid = np.array([gridrange for i in range(len(gridrange))])
    ygrid = xgrid.T

    fiber_r = np.sqrt((xgrid-dist_xy)**2 + (ygrid-dist_xy)**2)
    disthere = fiber_r <= 0.75

    grid_r = np.sqrt(xgrid**2 + ygrid**2)
    grid_r[~disthere] = np.nan

    psf_grid = fit_moffat(grid_r, amp, fwhm)
    mean_psf = np.nansum(psf_grid[disthere])*(INTSTEP**2)
    return mean_psf

def int_moffat(dist, amp, fwhm):
    """returns integrate_moffat() for an array of distances"""
    return [integrate_moffat(x, amp, fwhm) for x in dist]

def star_model(ra_dec, mid_ra, mid_dec, amp, fwhm):
    # ra_dec has to be a numpy array with shape [n_fibers, 2], where each fiber has an [ra, dec] entry.
    
    coords = SkyCoord(ra=ra_dec[:,0] * u.deg, dec=ra_dec[:,1] * u.deg)
    mid_coords = SkyCoord(ra=mid_ra*u.deg, dec=mid_dec*u.deg)
    r = mid_coords.separation(coords).arcsec
    
    return int_moffat(r, amp, fwhm)

def save_stars(gama_field, overwrite=False):
    
    config = GamaConfig()
   
    stars = ascii.read(config.stars_faint) #config.stars[gama_field])
    jt = get_table(gama_field)

    fiber_coords = SkyCoord(ra=jt["ra"]*u.deg, dec=jt["dec"]*u.deg)
    wl = config.wavelength
    wl_here = (wl > 4550) & (wl <= 4650)
    for star in stars:
        
        starID = star["ID"]
        if starID[0] != gama_field[-1]:
            continue
        
        outfile = os.path.join(config.basedir, "star_profiles/{}.tab".format(starID))
        if os.path.exists(outfile):
            if not overwrite:
                print("File exists, skipping "+starID)
                continue
            elif overwrite:
                print("Overwriting "+starID)
            
        mid_ra, mid_dec = star["ra"], star["dec"]

        star_coords = SkyCoord(ra=mid_ra*u.deg, dec=mid_dec*u.deg)
        f_dist = star_coords.separation(fiber_coords).arcsec
        f_here = f_dist < 8

        shotids = jt["shotid"][f_here]
        weights = jt["error"][f_here][:,wl_here] ** (-2)

        flux = jt["spectrum"][f_here][:,wl_here]

        print(flux.shape, weights.shape)

        weight_sum = np.nansum(weights, axis=1)
        flux_mean = np.nansum(weights*flux, axis=1) / weight_sum
        
        flux_mean_error = 1./np.sqrt(weight_sum)
        
        #empty = np.ones(flux.shape)
        #empty[~np.isfinite(flux)] = np.nan
        #N = np.nansum(empty, axis=1)
        #flux_minus_mean_squared = (flux.T - flux_mean.T).T**2
        #flux_mean_error = np.sqrt(np.nansum(weights * flux_minus_mean_squared, axis=1)/((N-1)*weight_sum))

        star_tab = {"flux": flux_mean,
                        "flux_error": flux_mean_error,
                        "ra": jt["ra"][f_here],
                        "dec": jt["dec"][f_here],
                        "shotid": shotids}

        ascii.write(star_tab, outfile, overwrite=overwrite)
    return 
            

class GamaConfig:
    
    def __init__(self):
        self.shotlist = {}
        __ = """20191203_0000024
20191221_0000022
20191203_0000025
20191222_0000023
20191221_0000023
20191224_0000024
20191222_0000024
20191229_0000023
20191231_0000024
20200101_0000019
20191231_0000025
20200101_0000020"""
        self.shotlist["gama09E"] = __.split()


        __ = """20200118_0000017
20200215_0000016
20200119_0000018
20200217_0000014
20200124_0000016
20200217_0000015
20200119_0000019
20200225_0000016
20200125_0000020
20200315_0000012
20200126_0000020"""
        self.shotlist["gama09F"] = __.split()

        self.basedir = "/home/idies/workspace/Temporary/maja/scratch/gama09"
        self.mxhf_base = "/home/idies/workspace/Storage/mxhf/gama09"
        self.gama_EF_catalog = os.path.join(self.mxhf_base, "data/gama09EF.cat")
        self.stars = {}
        self.stars["gama09E"] = os.path.join(self.basedir, "stars_gama09E.tab")
        self.stars["gama09F"] = os.path.join(self.basedir, "stars_gama09F.tab")
        self.stars_faint = os.path.join(self.basedir, "stars_faint.tab")
        self.survey = os.path.join(self.basedir, "survey.tab")    

        self.wavelength = np.arange(3470, 5542, 2)
