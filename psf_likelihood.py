import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d
import glob
import os

from astropy.table import Table, vstack
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u 
from astropy.coordinates import SkyCoord
from scipy.ndimage import gaussian_filter

from cobaya import likelihood

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
    mean_psf = np.nansum(psf_grid[disthere]) * (INTSTEP**2) # changed to a proper grid sum integration.
    return mean_psf

def int_moffat(dist, amp, fwhm):
    """returns integrate_moffat() for an array of distances"""
    return [integrate_moffat(x, amp, fwhm) for x in dist]

def star_model(ras, decs, mid_ra, mid_dec, amp, fwhm):
    
    coords = SkyCoord(ra=ras * u.deg, dec=decs * u.deg)
    mid_coords = SkyCoord(ra=mid_ra*u.deg, dec=mid_dec*u.deg)
    r = mid_coords.separation(coords).arcsec
    
    return int_moffat(r, amp, fwhm)

class Star:
	def __init__(self, starID, mid_ra, mid_dec, ra, dec, flux, flux_error, syst_error, shotid):
		self.ID = starID
		self.shotid = shotid
		self.mid_ra = mid_ra
		self.mid_dec = mid_dec
		self.ra = ra
		self.dec = dec
		self.flux = flux
		self.error = np.sqrt(flux_error**2 + syst_error**2)

		star_coords = SkyCoord(ra = self.mid_ra * u.deg, dec = self.mid_dec * u.deg)
		fiber_coords = SkyCoord(ra = self.ra * u.deg, dec = self.dec * u.deg)
		self.dist = star_coords.separation(fiber_coords).arcsec

basedir = "/work/05865/maja_n/stampede2/master/gama09"
mxhf_base = "/home/idies/workspace/Storage/mxhf/gama09"

wl = np.arange(3470, 5542, 2)

class PSFLike(likelihood.Likelihood):
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

		print("Analyzing field "+self.gama_field)
		self.shotlist = shotlist[self.gama_field]
		
		# get starting values for the star positions
		if self.gama_field == "gama09E":
		    stars = ascii.read(os.path.join(basedir, "stars_gama09E.tab"))
		elif self.gama_field == "gama09F":
		    stars = ascii.read(os.path.join(basedir, "stars_gama09F.tab"))
		else:
		    print("Error: no valid field specified: ", gama_field)	
		print(f"Found {len(stars)} stars.")

		self.starIDs = np.unique(stars["ID"])

		start_values = ascii.read(os.path.join(basedir, "star_ra_dec_amp.tab"))
		rms_values = ascii.read(os.path.join(basedir, "RMS_stars.tab"))

		# initiate stars

		self.stars = []
		for star in stars:
			#mid_ra, mid_dec = star["ra"], star["dec"]
			starID = star["ID"]
			mid_ra, mid_dec = start_values["ra"][start_values["ID"]==starID].data[0], start_values["dec"][start_values["ID"]==starID].data[0]

			star_tab = ascii.read(os.path.join(basedir, f"star_profiles/{starID}.tab"))

			for shotid in np.unique(star_tab["shotid"]):
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
				fiber_dists = fiber_dists[mask]
				star_here = star_here[mask]
				
				if len(fiber_dists[fiber_dists < 2.0]) < 3:
					print("Skipping star "+starID+" in "+shotid+" because of too few points within 2''.")
					continue

				rms = rms_values[(rms_values["ID"]==starID)&(rms_values["shotid"]==shotid)]["rms"][0]
				self.stars.append(Star( starID = starID,
						 shotid = shotid,
						 mid_ra = mid_ra,
						 mid_dec = mid_dec,
						 ra = star_here["ra"], #[shot_here],
						 dec = star_here["dec"], #[shot_here],
						 flux = star_here["flux"], #[shot_here],
						 flux_error = star_here["flux_error"],
						 syst_error = rms)) #[shot_here]))

	def logp(self, **kwargs):

		logp = 0
		for star in self.stars:
			#PSF = star_model(ras = star.ra,
			#		 decs = star.dec,
			#		 mid_ra = kwargs["ra_{}".format(star.ID)],
			#		 mid_dec = kwargs["dec_{}".format(star.ID)],
			#		 amp = kwargs["A_{}".format(star.ID)],
			#		 fwhm = kwargs["FWHM_{}".format(star.shotid)])

			PSF = int_moffat(dist = star.dist, amp = kwargs["A_{}".format(star.ID)], fwhm = kwargs["FWHM_{}".format(star.shotid)])

			logp += np.nansum(- 0.5*(star.flux - PSF)**2/star.error**2 - 0.5*np.log(2*np.pi*star.error**2))

		return logp

class PSFsingleLike(likelihood.Likelihood):
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

		self.shotid = str(self.shotid)
		self.shotid = self.shotid[:8]+"_"+self.shotid[8:]

		if str(self.shotid) not in shotlist[self.gama_field]:
			print("\nNo valid combination of shotid and gama_field.\n")
			return

		print("Analyzing shot "+self.shotid+" in field "+self.gama_field)
		self.shotlist = [self.shotid]
		
		# get starting values for the star positions
		if self.gama_field == "gama09E":
		    stars = ascii.read(os.path.join(basedir, "stars_gama09E.tab"))
		elif self.gama_field == "gama09F":
		    stars = ascii.read(os.path.join(basedir, "stars_gama09F.tab"))
		else:
		    print("Error: no valid field specified: ", gama_field)	
		print(f"Found {len(stars)} stars.")

		self.starIDs = np.unique(stars["ID"])

		start_values = ascii.read(os.path.join(basedir, "star_ra_dec_amp.tab"))

		# initiate stars

		self.stars = []
		for star in stars:
			#mid_ra, mid_dec = star["ra"], star["dec"]
			starID = star["ID"]
			mid_ra, mid_dec = start_values["ra"][start_values["ID"]==starID].data[0], start_values["dec"][start_values["ID"]==starID].data[0]

			star_tab = ascii.read(os.path.join(basedir, f"star_profiles/{starID}.tab"))

			if self.shotid in np.unique(star_tab["shotid"]):
				shotid = self.shotid
				shot_here = star_tab["shotid"] == shotid
				star_here = star_tab[shot_here]
				fiber_coords = SkyCoord(ra=star_here["ra"]*u.deg, dec=star_here["dec"]*u.deg)
				star_coord = SkyCoord(ra=mid_ra*u.deg, dec=mid_dec*u.deg)
				fiber_dists = star_coord.separation(fiber_coords).arcsec
				mask = fiber_dists < 5
				star_here = star_here[mask]
				self.stars.append(Star( starID = starID,
						 shotid = shotid,
						 mid_ra = mid_ra,
						 mid_dec = mid_dec,
						 ra = star_here["ra"], #[shot_here],
						 dec = star_here["dec"], #[shot_here],
						 flux = star_here["flux"], #[shot_here],
						 error = star_here["flux_error"])) #[shot_here]))

	def logp(self, **kwargs):

		logp = 0
		for star in self.stars:
			#PSF = star_model(ras = star.ra,
			#		 decs = star.dec,
			#		 mid_ra = kwargs["ra_{}".format(star.ID)],
			#		 mid_dec = kwargs["dec_{}".format(star.ID)],
			#		 amp = kwargs["A_{}".format(star.ID)],
			#		 fwhm = kwargs["FWHM_{}".format(star.shotid)])

			PSF = int_moffat(dist = star.dist, amp = kwargs["A_{}".format(star.ID)], fwhm = kwargs["FWHM_{}".format(star.shotid)])

			logp += np.nansum(- 0.5*(star.flux - PSF)**2/star.error**2 - 0.5*np.log(2*np.pi*star.error**2))

		return logp
