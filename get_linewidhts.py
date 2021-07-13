from scipy.optimize import curve_fit
from gama_tools import *
from astropy.coordinates import SkyCoord
import os

config = GamaConfig()
wl = config.wavelength

survey = ascii.read(config.survey)
basedir = config.basedir
gama_field = "gama09F"

jt = get_table(gama_field)

gold = ascii.read(os.path.join(basedir, "classified_golden_sample.cat"))
detcat = ascii.read(os.path.join(config.mxhf_base, "data", f"{gama_field}.cat"))
detcat["ifu"] = np.array(detcat["ifu"], dtype=int)
detcat["id"] = np.array(detcat["id"], dtype=int)

fiber_coords = SkyCoord(ra=jt["ra"]*u.deg, dec = jt["dec"]*u.deg)

def gaussian(x, mu, sigma, amp):
    return amp * np.exp(-0.5*(x-mu)**2/sigma**2)
    
all_spectra = {}
linewidths = {}
linewidth = {}
for lae in gold:
    if lae["field"] != gama_field:
        continue
        
    lae_info = detcat[(detcat["field"]==gama_field)&(detcat["ifu"] == int(lae["ifu"]))&(detcat["id"]==lae["id"])][0]
    
    ra, dec = lae_info["ra_com"], lae_info["dec_com"]
    lae_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    wave = lae_info["wl_com"]
    wl_here = abs(wl-wave) < 20
    
    fiber_dists = lae_coords.separation(fiber_coords).arcsec
    here = fiber_dists < 5
    
    jt_here = jt[here]
    fiber_dists = fiber_dists[here]
    
    spectra = []
    linewidths[(lae["field"], lae["ifu"], lae["id"])] = []
    for shotid in config.shotlist[gama_field]:
        if shotid == "20191229_0000023":
            continue
            
        shot_here = jt_here["shotid"]==shotid
        jt_tmp = jt_here[shot_here]
        
        fwhm = survey["fwhm"][survey["shotid"]==shotid][0]
        if fwhm > 2:
            continue
        
        fiber_tmp = fiber_dists[shot_here]
        jt_tmp = jt_tmp[fiber_tmp < 1.2*fwhm]
        fiber_tmp = fiber_tmp[fiber_tmp < 1.2*fwhm]
        
        psf = np.array(int_moffat(fiber_tmp, amp=1, fwhm=fwhm))
        
        spectrum = np.nansum((jt_tmp["spectrum"].T*psf.T).T, axis=0) / np.nansum(psf)
        spectrum_err = np.sqrt(np.nansum(((jt_tmp["error"]**2).T * psf.T**2).T, axis=0)) / np.nansum(psf)
        
        try:
            popt, pcov = curve_fit(gaussian, wl[wl_here], spectrum[wl_here], sigma=spectrum_err[wl_here], p0=[wave, 3.0, 0.5])
            if popt[-1] < 0:
                continue
            if popt[1] > 10:
                continue
            if (popt[0]==wave)&(popt[1]==3.0)&(popt[2]==0.5):
                continue
            
            linewidths[(lae["field"], lae["ifu"], lae["id"])].append(popt[1])
        except:
            pass
        
        spectra.append(spectrum)
        
    try:
        popt, pcov = curve_fit(gaussian, wl[wl_here], np.nanmean(spectra, axis=0)[wl_here], p0=[wave, 3.0, 0.5])
        wave = popt[0]
        linewidth[(lae["field"], lae["ifu"], lae["id"])] = popt[1]
        print(lae["field"], lae["ifu"], lae["id"], popt[1], popt[0])
        # I copied this into the linewidths.tab
    except Exception as e:
        print(e)
        print(lae["field"], lae["ifu"], lae["id"], "failed")
        pass
    
    all_spectra[(lae["field"], lae["ifu"], lae["id"])] = [wave, np.nanmean(spectra, axis=0)]
    
