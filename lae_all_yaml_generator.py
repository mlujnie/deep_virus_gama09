import os
import sys

troughsub = 'True'
fwhm_file = "survey"
redshift_cut = False

template = """likelihood:
  lae_likelihood_all.LAEallLike:
    python_path: /home/idies/workspace/Temporary/maja/scratch/gama09
    stop_at_error: True
    halo_profile: powerlaw
    simulated: False
    troughsub: {}
    fwhm_file: {}
    redshift_cut: {}
    input_params: [c1, c2, sigma_int]
params:
  c1:
    prior:
      min: 0
      max: 100
    ref:
      dist: norm
      loc: 5.0
      scale: 1.0
    proposal: 0.05
  c2:
    prior:
      min: -10.0
      max: 10.0
    ref:
      dist: norm
      loc: 0.0
      scale: 1.0
    proposal: 0.05
  sigma_int: 0.0
sampler:
  mcmc:
    Rminus1_stop: 0.01
    Rminus1_cl_stop: 0.1
    max_tries: 1000000
output: {}/{}
"""

if redshift_cut:
    path = "/home/idies/workspace/Temporary/maja/scratch/gama09/cobaya-chains/LAEs_all/low_z"
else:
    path = "/home/idies/workspace/Temporary/maja/scratch/gama09/cobaya-chains/LAEs_all/all_z"
try: 
    os.mkdir(path) 
except OSError as error: 
    print(error)     

with open(os.path.join(path,"all.yaml"), "w") as f:
    f.write(template.format(troughsub, fwhm_file, redshift_cut, path, "all"))

f.close()
print(f"Finished {1}/{1}.")
