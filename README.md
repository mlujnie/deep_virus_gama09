# deep_virus_gama09
Python scripts for the analysis of deep VIRUS data in the GAMA09 field for my master thesis.

## description of the scripts
* gama_tools.py - general setup and loading the data from the h5 files
* save radial profiles of LAEs, simulated LAEs, stars, and random locations: save_laes.py, save_simulated_laes.py, save_stars.py, save_bogus.py
* save average rest-frame spectra of the LAEs: save_core_spectra.py, LAE_sample_plot_and_average_spectra.ipynb
* prepare convolved PSF and power-law profiles and interpolate between them: get_psf_grid.py, get_powerlaw_grid.py, get_interpolation.py 
* get Lya line widths and centers of the LAEs: get_linewidhts.py
* individual LAE analysis: lae_likelihood.py, lae_yaml_generator.py
* joint LAE analysis and getting $\sigma_\mathrm{int}$: lae_likelihood_all.py, lae_all_yaml_generator.py, get_simga_int.py
* simulated LAE analysis: lae_likelihood.py, sim_lae_yaml_generator.py
* get empirical errors for star profiles, likelihood to get the PSF from star profiles: get_rms.py, psf_likelihood.py
* star halo profile analysis: star_halo_likelihood.py, star_halo_yaml_generator.py
* get background surface brightness from random pointings: get_sky_level.py
* master thesis plots: master_thesis_plots.ipynb
