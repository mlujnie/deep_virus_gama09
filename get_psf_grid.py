import numpy as np
import sys
import os
import time
from scipy.interpolate import interpn
from astropy.convolution import convolve_fft
from astropy.io import ascii
import multiprocessing
import gama_tools

fit_moffat = gama_tools.fit_moffat

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

def get_convolved_profile(fwhm):
	kernel_fwhm = fit_moffat(dist_map, amp=1, fwhm=fwhm)

	kernel_fiber = np.where(dist_map <= 0.75, 1, 0)
	kernel_fiber = kernel_fiber / (np.nansum(kernel_fiber)*INTSTEP**2)
	result_fiber = convolve_fft(kernel_fwhm, kernel_fiber)

	result_dict = {"r": dist_map[INDEX, INDEX:],
			"profile_raw": kernel_fwhm[INDEX, INDEX:],
			"profile_fiber": result_fiber[INDEX, INDEX:]}

	ascii.write(result_dict, "convolution_files/psf_int_{:.1f}.tab".format(fwhm), overwrite=True)
	return 1

fwhm_input = np.arange(0.6, 3.5, 0.1)

done = 0
N = len(fwhm_input)
for fwhm in fwhm_input:
	tmp = get_convolved_profile(fwhm)
	done += 1
	print(f"Done with {done}/{N}.")
