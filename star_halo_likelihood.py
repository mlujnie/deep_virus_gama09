import gama_tools as gt
import get_interpolation as gi

from astropy.io import ascii
import numpy as np
import astropy.units as u
import glob
import os

from astropy.coordinates import SkyCoord
from cobaya import likelihood

fit_moffat = gt.fit_moffat # dist, amp, fwhm
config = gt.GamaConfig()
wl = config.wavelength

basedir = "/home/idies/workspace/Temporary/maja/scratch/gama09"
#"/work/05865/maja_n/stampede2/master/gama09"
mxhf_base = "/home/idies/workspace/Storage/mxhf/gama09"

# get the psf and powerlaw functions
psf_func = gi.get_psf_func(basedir=basedir)
powerlaw_func = gi.get_powerlaw_func(basedir=basedir)

def lae_profile_powerlaw(rs, fwhm, c1, c2):
	redshift = 2.3 # so that r_cutoff = 1''.
	center_profile = c1 * psf_func([(fwhm, r) for r in rs])
	halo_profile = c2 * powerlaw_func([(fwhm, redshift, r) for r in rs])

	return center_profile + halo_profile

class Star:
	def __init__(self, starID, mid_ra, mid_dec, ra, dec, flux, flux_error, syst_error, shotid, fwhm):
		self.ID = starID
		self.shotid = shotid
		self.fwhm = fwhm
		self.mid_ra = mid_ra
		self.mid_dec = mid_dec
		self.ra = ra
		self.dec = dec
		self.flux = flux
		self.error = np.sqrt(flux_error**2 + syst_error**2)

		star_coords = SkyCoord(ra = self.mid_ra * u.deg, dec = self.mid_dec * u.deg)
		fiber_coords = SkyCoord(ra = self.ra * u.deg, dec = self.dec * u.deg)
		self.dist = star_coords.separation(fiber_coords).arcsec

class StarHaloLike(likelihood.Likelihood):
	def initialize(self):
		shotlist = {}
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
		shotlist["gama09E"] = __.split()


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
		shotlist["gama09F"] = __.split()

		print("Analyzing star "+self.starID)
		self.gama_field = "gama09"+self.starID[0]
		print("Analyzing field "+self.gama_field)
		self.shotlist = shotlist[self.gama_field]
		
		# get starting values for the star positions
		#if self.gama_field == "gama09E":
		#    stars = ascii.read(os.path.join(basedir, "stars_gama09E.tab"))
		#elif self.gama_field == "gama09F":
		#    stars = ascii.read(os.path.join(basedir, "stars_gama09F.tab"))
		#else:
		#    print("Error: no valid field specified: ", gama_field)
		stars = ascii.read(os.path.join(basedir, "stars_faint.tab"))
		print(f"Found {len(stars)} stars.")

		if self.fwhm_file == "survey":
			fwhm_file = "survey.tab"
		elif self.fwhm_file == "own":
			fwhm_file = "fwhm_posterior_mean.tab"
		else:
			print("Unknown option for fwhm_file: ", self.fwhm_file)
		fwhm_file = os.path.join(basedir, fwhm_file)
		fwhm_tab = ascii.read(fwhm_file)
		print("Using ", fwhm_file, " for the FWHM.")
		fwhm_dict = {}
		for i in range(len(fwhm_tab)):
			fwhm_dict[fwhm_tab[i]["shotid"]] = float(fwhm_tab[i]["fwhm"])

		print("Halo profile: "+self.halo_profile)

		self.starIDs = np.unique(stars["ID"])

		start_values = stars# ascii.read(os.path.join(basedir, "star_ra_dec_amp.tab"))
		rms_values = ascii.read(os.path.join(basedir, "RMS_stars.tab"))

		# initiate stars

		self.stars = []
		star = stars[stars["ID"] == self.starID]
		#for star in stars:
		starID = star["ID"].data[0]
		mid_ra, mid_dec = start_values["ra"][start_values["ID"]==starID].data[0], start_values["dec"][start_values["ID"]==starID].data[0]

		path = os.path.join(basedir, f"star_profiles/{starID}.tab")
		print("Reading path ", path)
		star_tab = ascii.read(os.path.join(basedir, f"star_profiles/{starID}.tab"))

		fiber_area = np.pi*0.75**2

		FWHM_cutoff = 1.5
		for shotid in np.unique(star_tab["shotid"]):
			#if fwhm_dict[shotid] > FWHM_cutoff:
			#	print("Skipping "+shotid+" because FWHM>{}''.".format(FWHM_cutoff))
			#	continue
			if shotid == "20191229_0000023":
				print("Skipping star in "+shotid)
				continue

			shot_here = star_tab["shotid"] == shotid
			star_here = star_tab[shot_here]
			fiber_coords = SkyCoord(ra=star_here["ra"]*u.deg, dec=star_here["dec"]*u.deg)
			star_coord = SkyCoord(ra=mid_ra*u.deg, dec=mid_dec*u.deg)
			fiber_dists = star_coord.separation(fiber_coords).arcsec
			mask = fiber_dists < 5
			mask *= np.isfinite(star_here["flux"])&np.isfinite(star_here["flux_error"])
			mask *= (star_here["flux"]!=0)&(star_here["flux_error"]!=0)
			fiber_dists = fiber_dists[mask]
			star_here = star_here[mask]
			
			if len(fiber_dists[fiber_dists < 2.0]) < 3:
				print("Skipping star "+starID+" in "+shotid+" because of too few points within 2''.")
				continue

			rms = rms_values[(rms_values["ID"]==starID)&(rms_values["shotid"]==shotid)]["rms"][0]
			rms = float(rms)
			self.stars.append(Star( starID = starID,
					 shotid = shotid,
					 fwhm = fwhm_dict[shotid],
					 mid_ra = mid_ra,
					 mid_dec = mid_dec,
					 ra = star_here["ra"],
					 dec = star_here["dec"], 
					 flux = star_here["flux"]/fiber_area,
					 flux_error = star_here["flux_error"]/fiber_area,
					 syst_error = rms/fiber_area)) 
				
	def logp_exp(self, **kwargs):

		logp = 0
		for star in self.stars:
			model = lae_profile(rs = star.dist, fwhm = star.fwhm, c1 = kwargs["c1"], c2 = kwargs["c2"], rn = kwargs["rn"]) 
			logp += np.nansum(- 0.5*(star.flux - model)**2/star.error**2 - 0.5*np.log(2*np.pi*star.error**2))

		return logp
	
	def logp_powerlaw(self, **kwargs):
		
		logp = 0
		for star in self.stars:
			model = lae_profile_powerlaw(rs = star.dist, fwhm = star.fwhm, c1 = kwargs["c1"], c2 = kwargs["c2"]) 
			logp += np.nansum(- 0.5*(star.flux - model)**2/star.error**2 - 0.5*np.log(2*np.pi*star.error**2))

		return logp
	

	def logp(self, **kwargs):
		if self.halo_profile == "exponential":
			return self.logp_exp(**kwargs)

		elif self.halo_profile == "powerlaw":
			return self.logp_powerlaw(**kwargs)

		else:
			print("ERROR: invalid halo_profile name. It must be 'exponential' or 'powerlaw'.")

