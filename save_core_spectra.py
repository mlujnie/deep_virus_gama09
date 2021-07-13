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

from gama_tools import *
config = GamaConfig()

shotid = "20191222_0000024"
filename = os.path.join(config.mxhf_base, "gama_recon/",shotid + ".h5")

lae_info = ascii.read("lae_info.tab")
lae_info["redshift"] = lae_info["wl_com"]/1215.67 - 1

def_wave = config.wavelength

perc=93

for gama_field in ["gama09E", "gama09F"]:

    laes_here = lae_info[lae_info["field"]==gama_field]
    print(len(laes_here), "LAEs in", gama_field)

    N=len(config.shotlist[gama_field])
    all_tables = []
    for idx in range(N):
        shotid = config.shotlist[gama_field][idx]
        fileh = tb.open_file(os.path.join(config.mxhf_base, "gama_recon/",shotid + ".h5"))
        t = Table(fileh.root.Info.read())
        fibers = fileh.root.Fibers
        t["shotid"] = [shotid for k in range(len(t))]
        for key in ["spectrum"]:
            t[key] = fibers[:][key]

        mask = np.ones(len(t), dtype="bool")
        spec_here = t["spectrum"]

        # exclude extreme continuum values

        wlcont_lo = (def_wave > 4000)&(def_wave <= 4500)
        medians_lo = np.nanmedian(spec_here[:,wlcont_lo], axis=1)
        perc_lo = np.nanpercentile(medians_lo, perc)

        wlcont_hi = (def_wave > 4800)&(def_wave <= 5300)
        medians_hi = np.nanmedian(spec_here[:,wlcont_hi], axis=1)
        perc_hi = np.nanpercentile(medians_hi, perc)
        mask[abs(medians_lo)>perc_lo] = False
        mask[abs(medians_hi)>perc_hi] = False

        spec_here = spec_here - np.nanmedian(spec_here[mask], axis=0) # this is how I ran it, but I should change it in the future.

        t["mask"] = mask
        print("Masked {}/{} fibers.".format(len(mask[~mask]), len(mask)))

        fiber_coords = SkyCoord(ra=t["ra"]*u.deg, dec=t["dec"]*u.deg)

        for lae in laes_here:
            lae_coords = SkyCoord(ra=lae["ra_com"]*u.deg, dec=lae["dec_com"]*u.deg)
            fiber_dists = lae_coords.separation(fiber_coords).arcsec

            here = fiber_dists < 2
            spec_0_2 = np.nanmedian(spec_here[here], axis=0)
            spec_0_2_m = np.nanmedian(spec_here[here*mask], axis=0)

            here = (fiber_dists >= 2) * (fiber_dists < 5)
            spec_2_5 = np.nanmedian(spec_here[here], axis=0)
            spec_2_5_m = np.nanmedian(spec_here[here*mask], axis=0)

            here = (fiber_dists >= 5) * (fiber_dists < 10)
            spec_5_10 = np.nanmedian(spec_here[here], axis=0)
            spec_5_10_m = np.nanmedian(spec_here[here*mask], axis=0)

            rest_wave = def_wave / (1+lae["redshift"])

            save_tab = {"rest_wave": rest_wave,
                        "spec_0_2": spec_0_2,
                        "spec_2_5": spec_2_5,
                        "spec_5_10": spec_5_10,
                        "spec_0_2_m": spec_0_2_m,
                        "spec_2_5_m": spec_2_5_m,
                        "spec_5_10_m": spec_5_10_m}
            save_name = "core_spectra_unflagged/{}_{}_{}_{}.tab".format(shotid, lae["field"], lae["ifu"], lae["id"])
            ascii.write(save_tab, save_name)
            print("Wrote to ", save_name)

        fileh.close()           

        print("Done with "+shotid+f" ({idx+1}/{N})")
print("Done.")