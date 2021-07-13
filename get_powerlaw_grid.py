import numpy as np
import sys
import os
import time
from scipy.interpolate import interpn
from astropy.convolution import convolve_fft
from astropy.io import ascii
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import multiprocessing
import gama_tools

# set up cosmology
cosmo = FlatLambdaCDM(H0=67.37, Om0=0.3147)

fit_moffat = gama_tools.fit_moffat

def lae_profile(r, r_cutoff):
        return np.where(r>r_cutoff, (r/r_cutoff)**(-2.4), 1)

max_radius = 50
diff = 0.1
INTSTEP = diff
tmp = np.arange(-max_radius, max_radius+diff, diff)
print(len(tmp))
INDEX = np.argmin(abs(tmp))
print("Middle row:", INDEX)
x = np.array([tmp for i in range(len(tmp))])
y = x.T
dist_map = np.sqrt(x**2+y**2)

def get_convolved_profile(fwhm, redshift):
	kernel = fit_moffat(dist_map, amp=1, fwhm=fwhm)

	kpc_per_arcsec_here = cosmo.kpc_proper_per_arcmin(redshift)/60*u.arcmin/u.kpc
	r_cutoff = kpc_proper_arcsec_mid / kpc_per_arcsec_here

	lae_img = lae_profile(dist_map, r_cutoff)
	result = convolve_fft(lae_img, kernel)

	kernel_fiber = np.where(dist_map <= 0.75, 1, 0)
	kernel_fiber = kernel_fiber / (np.nansum(kernel_fiber)*INTSTEP**2)
	result_fiber = convolve_fft(result, kernel_fiber)


	result_dict = {"r": dist_map[INDEX, INDEX:],
			"profile_raw": result[INDEX, INDEX:],
			"profile_fiber": result_fiber[INDEX, INDEX:]}

	ascii.write(result_dict, "convolution_files/powerlaw_conv_int_{:.3f}_{:.1f}.tab".format(fwhm, redshift), overwrite=True)
	return 1

z_min = 1.9
z_max = 3.5
z_mid = 2.3
kpc_proper_arcsec_mid = cosmo.kpc_proper_per_arcmin(z_mid)/60*u.arcmin/u.kpc
print('r_cutoff [kpc] = ', kpc_proper_arcsec_mid)
sys.exit()

fwhm_input = np.arange(0.6, 3.5, 0.1)
z_input = np.arange(z_min, z_max+0.1, 0.1)

done = 0
N = len(fwhm_input)*len(z_input)
for fwhm in fwhm_input:
	for redshift in z_input:
		tmp = get_convolved_profile(fwhm, redshift)
		done += 1
		print(f"Done with {done}/{N}.")
