import numpy as np
import glob
from astropy.io import ascii
from astropy.table import vstack
from astropy.stats import biweight_location, biweight_scale

bogus = glob.glob("bogus_laes_troughsub/*.dat")
all_tabs = []
i = 0
N = len(bogus)
for fin in bogus:
    tmp = ascii.read(fin)
    tmp["mask"] = np.array(tmp["mask"], dtype=bool)
    tmp["shotid"] = np.array(tmp["shotid"], dtype=str)
    #tmp = tmp[tmp["mask"]]
    tmp = tmp[tmp["flux"]!=0]
    tmp = tmp[tmp["flux_error"]!=0]
    all_tabs.append(tmp)
    i+=1
    if i%100==0:
        print(f"Done with {i}/{N}.")
        
big_tab = vstack(all_tabs)
all_values = big_tab["flux"][big_tab['mask']]
flux_sky_level = np.nanmean(all_values)
flux_sky_level_median = np.nanmedian(all_values)
flux_sky_level_err = biweight_scale(all_values) / np.sqrt(len(all_values))
print("Sky level (flux units): ", flux_sky_level)
fiber_area = np.pi*0.75**2
print("Sky level (surface brightness units): ", flux_sky_level/fiber_area)
print("\nSky level median (flux units): ", flux_sky_level_median)
print("Sky level median (surface brightness units): ", flux_sky_level_median/fiber_area)
print('\n sky level error (surface brightness units): ', flux_sky_level_err/fiber_area)

print('\n.........................................................................\n')

all_values = big_tab["flux_troughsub"]
print('After trough subtraction and without masking:')
flux_sky_level = np.nanmean(all_values)
flux_sky_level_median = np.nanmedian(all_values)
flux_sky_level_err = biweight_scale(all_values) / np.sqrt(len(all_values))
print("Sky level (flux units): ", flux_sky_level)
fiber_area = np.pi*0.75**2
print("Sky level (surface brightness units): ", flux_sky_level/fiber_area)
print("\nSky level median (flux units): ", flux_sky_level_median)
print("Sky level median (surface brightness units): ", flux_sky_level_median/fiber_area)
print('\n sky level error (surface brightness units): ', flux_sky_level_err/fiber_area)