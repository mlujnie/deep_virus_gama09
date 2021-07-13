import os
import sys
import glob

TROUGHSUB = True
if TROUGHSUB:
    DIR_APX = '_troughsub'
    troughsub = 'True'
else:
    DIR_APX = ''
    troughsub = 'False'

lae_ids = glob.glob(f"lae_simulated{DIR_APX}/*")
lae_ids = [x.split("/")[-1][:-4] for x in lae_ids]

template = """likelihood:
    lae_likelihood.LAELike:
      python_path: /home/idies/workspace/Temporary/maja/scratch/gama09
      stop_at_error: True
      laeID: {}
      simulated: True
      halo_profile: powerlaw
      troughsub: {}
      fwhm_file: {}
      input_params: [c1, c2]
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
sampler:
  mcmc:
    Rminus1_stop: 0.01
    Rminus1_cl_stop: 0.1
    max_tries: 1000000
output: {}/{}"""

fwhm_file = "survey"

inum = 0
N = len(lae_ids)
for laeID in lae_ids:
	path = "/home/idies/workspace/Temporary/maja/scratch/gama09/cobaya-chains/simulated_LAEs/{}".format(laeID)
	try: 
		os.mkdir(path) 
	except OSError as error: 
		print(error)     

	with open(os.path.join(path,"lae_profile.yaml"), "w") as f:
		f.write(template.format(laeID, troughsub, fwhm_file, path, laeID))

	f.close()
	inum += 1
	print(f"Finished {inum}/{N}.")
