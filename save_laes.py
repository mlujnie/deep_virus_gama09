from gama_tools import *

config = GamaConfig()
wl = config.wavelength
def_wave = wl
lae_dir = os.path.join(config.basedir, "lae_data/*")

gama_field = "gama09E"
jt = get_table(gama_field, clean=True, perc=90)

fiber_coords = SkyCoord(ra=jt["ra"]*u.deg, dec=jt["dec"]*u.deg)

detcat = ascii.read(os.path.join(config.mxhf_base, "data", f"{gama_field}.cat"))
detcat["ifu"] = np.array(detcat["ifu"], dtype=int)
detcat["id"] = np.array(detcat["id"], dtype=int)
linewidth = ascii.read(os.path.join(config.basedir, "linewidths.tab"))

for lae in linewidth:
    if lae["field"] != gama_field:
        continue
        
    lae_info = detcat[(detcat["field"]==gama_field)&(detcat["ifu"] == int(lae["ifu"]))&(detcat["id"]==lae["id"])][0]
    
    ra, dec = lae_info["ra_com"], lae_info["dec_com"]
    lae_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    wave = lae["wavelength"]
    
    wl_here = abs(wl-wave) <= 1.5*lae["linewidth"]
    
    fiber_dists = lae_coords.separation(fiber_coords).arcsec
    here = fiber_dists < 10
    
    jt_here = jt[here]
    fiber_dists = fiber_dists[here]
    
    ones = np.ones(jt_here['spectrum'].shape)
    ones[~np.isfinite(jt_here['spectrum'])] = 0
       
    cont_wlhere = (abs(def_wave - wave) <= 40) & (abs(def_wave - wave)>2.5*lae['linewidth'])
    trough_continuum = np.nanmedian(jt_here['spectrum'][:,cont_wlhere], axis=1)
    N = np.nansum(ones[:,cont_wlhere], axis=1)
    trough_continuum_error = biweight_scale(jt_here['spectrum'][:,cont_wlhere], axis=1) / np.sqrt(N)
    trough_contsub = (jt_here['spectrum'].T - trough_continuum).T
    trough_contsub_error = np.sqrt(jt_here['error'].T**2 + trough_continuum_error**2).T
    
    ras, decs, fluxes, flux_errors, shots = [], [], [], [], []
    troughsubs, err_troughsubs = [], []
    mask = []
    for shotid in config.shotlist[gama_field]:
        if shotid == "20191229_0000023":
            continue
            
        shot_here = jt_here["shotid"]==shotid
        jt_tmp = jt_here[shot_here]
        trough_contsub_tmp = trough_contsub[shot_here]
        trough_contsub_error_tmp = trough_contsub_error[shot_here]
        
        spec_sum_tmp = np.nansum(jt_tmp["spectrum"][:,wl_here], axis=1)
        err_sum_tmp = np.sqrt(np.nansum(jt_tmp["error"][:,wl_here]**2, axis=1))
        troughsub_sum_tmp = np.nansum(trough_contsub_tmp[:,wl_here], axis=1)
        err_troughsub_sum_tmp = np.sqrt(np.nansum(trough_contsub_tmp[:,wl_here]**2, axis=1))
        
        here = np.isfinite(spec_sum_tmp) & np.isfinite(err_sum_tmp) & (spec_sum_tmp!=0) & (err_sum_tmp!=0) & np.isfinite(troughsub_sum_tmp) & np.isfinite(err_troughsub_sum_tmp) & (troughsub_sum_tmp!=0) & (err_troughsub_sum_tmp!=0)
        # should be there! But was not there when I ran it.
        spec_sum_tmp = spec_sum_tmp[here]
        err_sum_tmp = err_sum_tmp[here]
        jt_tmp = jt_tmp[here]
        troughsub_sum_tmp = troughsub_sum_tmp[here]
        err_troughsub_sum_tmp = err_troughsub_sum_tmp[here]
        
        ras.append(jt_tmp["ra"])
        decs.append(jt_tmp["dec"])
        fluxes.append(spec_sum_tmp)
        flux_errors.append(err_sum_tmp)
        shots.append(jt_tmp["shotid"])
        mask.append(np.array(jt_tmp["mask"], dtype=int))
        troughsubs.append(troughsub_sum_tmp)
        err_troughsubs.append(err_troughsub_sum_tmp)
        
    tab = {"flux": np.concatenate(fluxes), "flux_error": np.concatenate(flux_errors), "ra": np.concatenate(ras), "dec": np.concatenate(decs), 
           "shotid": np.concatenate(shots),
          "mask":np.concatenate(mask),
          "flux_troughsub": np.concatenate(troughsubs),
          "flux_troughsub_error": np.concatenate(err_troughsubs)}
    save_file = os.path.join(config.basedir, "lae_data_troughsub", "{}_{}_{}.dat".format(lae["field"], lae["ifu"], lae["id"]))
    ascii.write(tab, save_file)
    print("Wrote to "+save_file)
