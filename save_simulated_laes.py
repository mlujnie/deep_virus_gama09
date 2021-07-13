from gama_tools import *
import random
import get_interpolation as gi

config = GamaConfig()
wl = config.wavelength
def_wave = wl
lae_dir = os.path.join(config.basedir, "lae_data/*")
basedir = config.basedir

gama_field = "gama09F"
jt = get_table(gama_field, clean=True, perc=90)

fiber_coords = SkyCoord(ra=jt["ra"]*u.deg, dec=jt["dec"]*u.deg)

detcat = ascii.read(os.path.join(config.mxhf_base, "data", f"{gama_field}.cat"))
detcat["ifu"] = np.array(detcat["ifu"], dtype=int)
detcat["id"] = np.array(detcat["id"], dtype=int)
linewidth = ascii.read(os.path.join(config.basedir, "linewidths.tab"))

# get the psf and powerlaw functions
psf_func = gi.get_psf_func(basedir=basedir)
powerlaw_func = gi.get_powerlaw_func(basedir=basedir)

def lae_profile_powerlaw(rs, fwhm, c1, c2, redshift):
	center_profile = c1 * psf_func([(fwhm, r) for r in rs])
	halo_profile = c2 * powerlaw_func([(fwhm, redshift, r) for r in rs])

	return center_profile + halo_profile


fwhm_tab = ascii.read(os.path.join(basedir, "survey.tab"))
fwhm_dict = {}
for i in range(len(fwhm_tab)):
        fwhm_dict[fwhm_tab[i]["shotid"]] = float(fwhm_tab[i]["fwhm"])
        
LINEWIDTH = 3.0
c1 = 2.5 # change this!!! It might not be a good number anymore.

sim_lae_info = {}
sim_lae_info["field"] = []
sim_lae_info["ifu"] = []
sim_lae_info["id"] = []
sim_lae_info["ra_com"] = []
sim_lae_info["dec_com"] = []
sim_lae_info["wl_com"] = []
sim_lae_info["redshift"] = []

fiber_area = np.pi*0.75**2
for c2 in [0, 0.05, 0.1, 0.2, 0.5]:
    for iternum in range(100):

        fibnum = random.randint(0, len(jt))
        while ~jt["mask"][fibnum]:
            print("Recalculating fibnum.")
            fibnum = random.randint(0, len(jt))

        dra, ddec = random.uniform(-1, 1)/3600., random.uniform(-1, 1)/3600.
        ra, dec = jt["ra"][fibnum]+dra, jt[fibnum]["dec"]+ddec
        lae_coords = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        wave = random.uniform(3470, 5540)
        redshift = wave/1215.67 - 1

        wl_here = abs(wl-wave) <= 1.5*LINEWIDTH

        fiber_dists = lae_coords.separation(fiber_coords).arcsec
        here = fiber_dists < 10

        jt_here = jt[here]
        fiber_dists = fiber_dists[here]
        
        ones = np.ones(jt_here['spectrum'].shape)
        ones[~np.isfinite(jt_here['spectrum'])] = 0

        cont_wlhere = (abs(def_wave - wave) <= 40) & (abs(def_wave - wave)>2.5*LINEWIDTH)
        trough_continuum = np.nanmedian(jt_here['spectrum'][:,cont_wlhere], axis=1)
        N = np.nansum(ones[:,cont_wlhere], axis=1)
        trough_continuum_error = biweight_scale(jt_here['spectrum'][:,cont_wlhere], axis=1) / np.sqrt(N)
        trough_contsub = (jt_here['spectrum'].T - trough_continuum).T
        trough_contsub_error = np.sqrt(jt_here['error'].T**2 + trough_continuum_error**2).T

        ras, decs, fluxes, flux_errors, shots = [], [], [], [], []
        troughsubs, err_troughsubs = [], []
        mask = []
        
        for shotid in np.unique(jt_here["shotid"]):#config.shotlist[gama_field]:
            if shotid == "20191229_0000023":
                continue

            shot_here = jt_here["shotid"]==shotid
            jt_tmp = jt_here[shot_here]
            trough_contsub_tmp = trough_contsub[shot_here]
            trough_contsub_error_tmp = trough_contsub_error[shot_here]
            
            fwhm = fwhm_dict[shotid]
            sim_lae = lae_profile_powerlaw(rs=fiber_dists[shot_here], fwhm=fwhm, c1=c1, c2=c2, redshift=redshift)
            sim_lae = sim_lae * fiber_area

            spec_sum_tmp = np.nansum(jt_tmp["spectrum"][:,wl_here], axis=1)
            err_sum_tmp = np.sqrt(np.nansum(jt_tmp["error"][:,wl_here]**2, axis=1))
            troughsub_sum_tmp = np.nansum(trough_contsub_tmp[:,wl_here], axis=1)
            err_troughsub_sum_tmp = np.sqrt(np.nansum(trough_contsub_tmp[:,wl_here]**2, axis=1))

            here = np.isfinite(spec_sum_tmp) & np.isfinite(err_sum_tmp) & (spec_sum_tmp!=0) & (err_sum_tmp!=0) & np.isfinite(troughsub_sum_tmp) & np.isfinite(err_troughsub_sum_tmp) & (troughsub_sum_tmp!=0) & (err_troughsub_sum_tmp!=0)
            spec_sum_tmp = spec_sum_tmp[here]
            err_sum_tmp = err_sum_tmp[here]
            jt_tmp = jt_tmp[here]
            sim_lae = sim_lae[here]
            troughsub_sum_tmp = troughsub_sum_tmp[here]
            err_troughsub_sum_tmp = err_troughsub_sum_tmp[here]

            ras.append(jt_tmp["ra"])
            decs.append(jt_tmp["dec"])
            fluxes.append(spec_sum_tmp + sim_lae)
            flux_errors.append(err_sum_tmp)
            shots.append(jt_tmp["shotid"])
            mask.append(np.array(jt_tmp["mask"], dtype=int))
            troughsubs.append(troughsub_sum_tmp + sim_lae)
            err_troughsubs.append(err_troughsub_sum_tmp)

        tab = {"flux": np.concatenate(fluxes), "flux_error": np.concatenate(flux_errors), "ra": np.concatenate(ras), "dec": np.concatenate(decs), "shotid": np.concatenate(shots),
              "mask":np.concatenate(mask),
              "flux_troughsub": np.concatenate(troughsubs),
              "flux_troughsub_error": np.concatenate(err_troughsubs)}
        if len(tab["flux"]) == 0:
            continue
        save_file = os.path.join(config.basedir, "lae_simulated_troughsub", "{}_{}_{}.dat".format(gama_field, int(100*c2), iternum))
        sim_lae_info["field"].append(gama_field)
        sim_lae_info["ifu"].append(int(100*c2))
        sim_lae_info["id"].append(iternum)
        sim_lae_info["ra_com"].append(ra)
        sim_lae_info["dec_com"].append(dec)
        sim_lae_info["wl_com"].append(wave)
        sim_lae_info["redshift"].append(redshift)
        ascii.write(tab, save_file)
        print("Wrote to "+save_file)
        
ascii.write(sim_lae_info, os.path.join(config.basedir, f"sim_lae_info_troughsub_{gama_field}.tab"))
print("Wrote info to ", os.path.join(config.basedir, f"sim_lae_info_troughsub_{gama_field}.tab"))
