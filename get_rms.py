from astropy.io import ascii
import numpy as np
import glob
from scipy.optimize import curve_fit
from astropy.coordinates import SkyCoord
import astropy.units as u
import os
import get_interpolation as gi
import matplotlib.pyplot as plt

basedir = "/home/idies/workspace/Temporary/maja/scratch/gama09"
mxhf_base = "/home/idies/workspace/Storage/mxhf/gama09"

wl = np.arange(3470, 5542, 2)

psf_func = gi.get_psf_func(basedir)
powerlaw_func = gi.get_powerlaw_func(basedir)

fiber_area = np.pi*0.75**2
def lae_profile_powerlaw(rs, fwhm, c1, c2):
	redshift = 2.3
	center_profile = c1 * psf_func([(fwhm, r) for r in rs])
	halo_profile = c2 * powerlaw_func([(fwhm, redshift, r) for r in rs])

	return center_profile + halo_profile


ff = glob.glob("star_profiles/*_faint_*.tab")
ff = np.sort(ff)
print("Found ", len(ff), " faint star files.")

star_info = ascii.read("stars_faint.tab")

fwhm_dict = {}
fwhm_tab = ascii.read(os.path.join(basedir, "survey.tab")) # change this!
for i in range(len(fwhm_tab)):
	fwhm_dict[fwhm_tab[i]["shotid"]] = float(fwhm_tab[i]["fwhm"])
    
r_plot = np.arange(0,8,0.1)

rms_dict = {
    "ID" : [],
    "shotid": [],
    "rms": []
}
i=0
for fin in ff:
    tab = ascii.read(fin)
    tab = tab[np.isfinite(tab["flux"])]
    starID = fin.split("/")[-1][:-4]
    si = star_info[star_info["ID"] == starID]
    ra, dec = si["ra"].data[0], si["dec"].data[0]
    star_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    fiber_coords = SkyCoord(ra=tab["ra"]*u.deg, dec=tab["dec"]*u.deg)
    fiber_dists = star_coord.separation(fiber_coords).arcsec

    for shotid in np.unique(tab["shotid"]):
        if shotid == "20191229_0000023":
            continue
        shot_here = tab["shotid"]==shotid
        tmptab = tab[shot_here]
        tmp_rs = fiber_dists[shot_here]
        fwhm = fwhm_dict[shotid]
        
        popt, pcov = curve_fit(lae_profile_powerlaw, tmp_rs, tmptab["flux"], sigma=tmptab["flux_error"])
        #print(popt)
    
        fwhm_here = tmp_rs < 1.5*fwhm
        rms = np.sqrt(np.nanmean((tmptab["flux"][fwhm_here]-lae_profile_powerlaw(tmp_rs[fwhm_here], *popt))**2))
        
        rms_dict["ID"].append(starID)
        rms_dict["shotid"].append(shotid)
        rms_dict["rms"].append(rms)
        
    i+=1
    if i%10 == 0:
        print(f"Finished {i}/79.")
        
#ascii.write(rms_dict, "RMS_stars.tab")