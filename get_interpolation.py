from astropy.io import ascii
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import glob
import os

def get_psf_func(basedir):
	# read all PSF files
	ff = glob.glob(os.path.join(basedir, "convolution_files/psf_int_???.tab"))
	fwhm_psf = np.sort(np.unique([float(fin.split("_")[-1][:-4]) for fin in ff]))
	files = []
	if len(ff)==0:
		print("Wrong base directory: ", basedir)
		return
	for fwhm in fwhm_psf:
		tmp = ascii.read(os.path.join(basedir,"convolution_files/psf_int_{:.1f}.tab".format(fwhm)))
		files.append(tmp["profile_fiber"].data)

	files = np.array(files)
	points = (fwhm_psf, tmp["r"])
	psf_func = RegularGridInterpolator(points, files, bounds_error=False) # will be extrapolated
	return psf_func

def get_powerlaw_func(basedir):
	# read all powerlaw files
	ff = glob.glob(os.path.join(basedir, "convolution_files/powerlaw_conv_int_?????_???.tab"))
	fwhm_psf = np.sort(np.unique([float(fin.split("_")[-2]) for fin in ff]))
	redshift = np.sort(np.unique([float(fin.split("_")[-1][:-4]) for fin in ff]))
	files = []
	if len(ff)==0:
		print("Wrong base directory: ", basedir)
		return
	for fwhm in fwhm_psf:
		tmp_files = []
		for z in redshift:
			tmp = ascii.read(os.path.join(basedir, "convolution_files/powerlaw_conv_int_{:.3f}_{:.1f}.tab".format(fwhm, z)))
			tmp_files.append(tmp["profile_fiber"].data)
		files.append(tmp_files)

	files = np.array(files)
	points = (fwhm_psf, redshift, tmp["r"])
	powerlaw_func = RegularGridInterpolator(points, files, bounds_error=False) # will be extrapolated
	return powerlaw_func

