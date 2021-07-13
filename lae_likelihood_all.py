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

class LAE:
	def __init__(self, laeID, mid_ra, mid_dec, ra, dec, flux, flux_error, shotid, fwhm, z):
		self.ID = laeID
		self.shotid = shotid
		self.fwhm = fwhm
		self.mid_ra = mid_ra
		self.mid_dec = mid_dec
		self.ra = ra
		self.dec = dec
		self.z = z
		self.flux = flux
		self.error = flux_error 

		star_coords = SkyCoord(ra = self.mid_ra * u.deg, dec = self.mid_dec * u.deg)
		fiber_coords = SkyCoord(ra = self.ra * u.deg, dec = self.dec * u.deg)
		self.dist = star_coords.separation(fiber_coords).arcsec

basedir = "/home/idies/workspace/Temporary/maja/scratch/gama09"
#"/work/05865/maja_n/stampede2/master/gama09"
mxhf_base = "/home/idies/workspace/Storage/mxhf/gama09"

# get the psf and powerlaw functions
psf_func = gi.get_psf_func(basedir=basedir)
powerlaw_func = gi.get_powerlaw_func(basedir=basedir)

def lae_profile_powerlaw(rs, fwhm, c1, c2, redshift):
	center_profile = c1 * psf_func([(fwhm, r) for r in rs])
	halo_profile = c2 * powerlaw_func([(fwhm, redshift, r) for r in rs])

	return center_profile + halo_profile

#############
# Detectids #
#############
lae_ids = """gama09E_104_45
gama09E_106_34
gama09E_16_143
gama09E_23_195
gama09E_24_60
gama09E_24_66
gama09E_28_67
gama09E_33_41
gama09E_33_7
gama09E_34_51
gama09E_35_24
gama09E_35_52
gama09E_36_199
gama09E_36_50
gama09E_43_23
gama09E_44_28
gama09E_46_40
gama09E_47_33
gama09E_47_72
gama09E_68_27
gama09E_68_35
gama09E_72_131
gama09E_75_62
gama09E_77_48
gama09E_83_62
gama09E_84_135
gama09E_84_34
gama09E_85_41
gama09E_86_217
gama09E_88_28
gama09E_94_94
gama09E_95_49
gama09E_97_101
gama09F_105_203
gama09F_106_125
gama09F_106_31
gama09F_13_47
gama09F_14_113
gama09F_16_34
gama09F_24_98
gama09F_27_87
gama09F_34_170
gama09F_37_121
gama09F_37_332
gama09F_38_88
gama09F_42_72
gama09F_43_197
gama09F_44_62
gama09F_44_63
gama09F_45_107
gama09F_53_192
gama09F_67_98
gama09F_76_199
gama09F_77_136
gama09F_81_68
gama09F_82_84
gama09F_87_131
gama09F_87_384
gama09F_87_418
gama09F_87_85
gama09F_91_52
gama09F_97_318""".split("\n")

class LAEallLike(likelihood.Likelihood):
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
            
		self.shotlist = shotlist

		print("Halo profile: "+self.halo_profile)
		print("Type self.simulated: ", type(self.simulated))

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
            
		if self.troughsub:
			print('Using locally continuum-subtracted fluxes.')
			DIR_APX = '_troughsub'
		else:
			print('No local continuum subtraction.')
			DIR_APX = ''

		if self.simulated:
			lae_info = ascii.read(os.path.join(basedir, f"sim_lae_info{DIR_APX}.tab"))
			fin = os.path.join(basedir, f"lae_simulated{DIR_APX}/{self.laeID}.dat")
		else:
			lae_info = ascii.read(os.path.join(basedir, "lae_info.tab"))
			#fin = os.path.join(basedir, f"lae_data/{self.laeID}.dat")
		
		print("Redshift cut? ", self.redshift_cut)
		if self.redshift_cut:
			print("All LAEs: ", len(lae_info))
			lae_info["mask"] = (lae_info["wl_com"]/1215.67 - 1) <= 2.6
			print("LAEs left after redshift cut: ", len(lae_info[lae_info["mask"]]))
        
		self.laes = []
		fiber_area = np.pi*0.75**2
		self.NUMBER_OF_DATA_POINTS = 0
		for laeID in lae_ids:
				fin = os.path.join(basedir, f"lae_data{DIR_APX}/{laeID}.dat")
				print("Reading "+laeID)
                
				field, ifu, id = laeID.split("_")
				linf = lae_info[(lae_info["field"]==field)&(np.array(lae_info["ifu"], dtype=int)==int(ifu))&(np.array(lae_info["id"], dtype=int)==int(id))][0]
				mid_ra, mid_dec = linf["ra_com"], linf["dec_com"]
				mid_wave = linf["wl_com"]
				redshift = mid_wave/1215.67 - 1
				if self.redshift_cut:
					if redshift > 2.6:
						print("LAE has redshift z>2.6: ", redshift)
						continue
                        
				lae_tab = ascii.read(fin)
				lae_tab["mask"] = np.array(lae_tab["mask"], dtype=bool)
				lae_tab["mask"] *= lae_tab["flux"]!= 0.0
                
                
				if self.troughsub:
					print("Not masking any fibers.")
					print("Using extensions "+f"flux{DIR_APX}"+" and "+f"flux{DIR_APX}_error")
				else:
					print("Masking continuum fibers.")
					lae_tab = lae_tab[lae_tab["mask"]]
					print("Using extensions "+f"flux{DIR_APX}"+" and "+f"flux{DIR_APX}_error")


				#FWHM_cutoff = 1.8
				for shotid in np.unique(lae_tab["shotid"]):
				#		if fwhm_dict[shotid] > FWHM_cutoff:
				#			print("Skipping "+shotid+" because FWHM>{}''.".format(FWHM_cutoff))
				#			continue
						if shotid == "20191229_0000023":
								print("Skipping LAE in 20191229_0000023")
								continue
                                      
						shot_here = lae_tab["shotid"]==shotid
						lae_here = lae_tab[shot_here]
						self.laes.append(LAE( laeID = laeID,
                                                 shotid = shotid,
                                                 fwhm = fwhm_dict[shotid],
                                                 mid_ra = mid_ra,
                                                 mid_dec = mid_dec,
                                                 ra = lae_here["ra"],
                                                 dec = lae_here["dec"],
						 z = redshift,
                                                 flux = lae_here[f"flux{DIR_APX}"]/fiber_area, 
                                                 flux_error = lae_here[f"flux{DIR_APX}_error"]/fiber_area))
						self.NUMBER_OF_DATA_POINTS += len(lae_here)
		print("Number of data points: ", self.NUMBER_OF_DATA_POINTS)
            
	def logp_exp(self, **kwargs):

		logp = 0
		for lae in self.laes:
			model = lae_profile(rs = lae.dist, fwhm = lae.fwhm, c1 = kwargs[f"c1_{lae.ID}"], c2 = kwargs[f"c2_{lae.ID}"], rn = kwargs[f"rn_{lae.ID}"]) 
			logp += np.nansum(- 0.5*(lae.flux - model)**2/lae.error**2 - 0.5*np.log(2*np.pi*lae.error**2))

		return logp
	
	def logp_powerlaw(self, **kwargs):
		logp = 0
		for lae in self.laes:
			model = lae_profile_powerlaw(rs = lae.dist, fwhm = lae.fwhm, c1 = kwargs["c1"], c2 = kwargs["c2"], redshift=lae.z) 
			logp += np.nansum(- 0.5*(lae.flux - model)**2/(lae.error**2+kwargs["sigma_int"]**2) - 0.5*np.log(2*np.pi*lae.error**2))
		return logp
    
	
	def logp(self, **kwargs):
		if self.halo_profile == "exponential":
			return self.logp_exp(**kwargs)

		elif self.halo_profile == "powerlaw":
			return self.logp_powerlaw(**kwargs)

		else:
			print("ERROR: invalid halo_profile name. It must be 'exponential' or 'powerlaw'.")



