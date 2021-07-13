import os
import sys

lae_ids = """E_faint_1
E_faint_10
E_faint_11
E_faint_12
E_faint_13
E_faint_14
E_faint_15
E_faint_16
E_faint_17
E_faint_18
E_faint_19
E_faint_2
E_faint_20
E_faint_21
E_faint_22
E_faint_23
E_faint_24
E_faint_25
E_faint_26
E_faint_27
E_faint_28
E_faint_29
E_faint_3
E_faint_30
E_faint_31
E_faint_32
E_faint_33
E_faint_34
E_faint_35
E_faint_36
E_faint_37
E_faint_4
E_faint_5
E_faint_6
E_faint_7
E_faint_8
E_faint_9
F_faint_1
F_faint_10
F_faint_11
F_faint_12
F_faint_13
F_faint_14
F_faint_15
F_faint_16
F_faint_17
F_faint_18
F_faint_19
F_faint_2
F_faint_20
F_faint_21
F_faint_22
F_faint_23
F_faint_24
F_faint_25
F_faint_26
F_faint_27
F_faint_28
F_faint_29
F_faint_3
F_faint_30
F_faint_31
F_faint_32
F_faint_33
F_faint_34
F_faint_35
F_faint_36
F_faint_37
F_faint_38
F_faint_39
F_faint_4
F_faint_40
F_faint_41
F_faint_42
F_faint_5
F_faint_6
F_faint_7
F_faint_8
F_faint_9""".split("\n")

flagged = """E_faint_10
E_faint_12
E_faint_13
E_faint_17
E_faint_19
E_faint_21
E_faint_25
E_faint_26
E_faint_3
E_faint_33
F_faint_10
F_faint_11
F_faint_13
F_faint_21
F_faint_27
F_faint_28
F_faint_29
F_faint_32
F_faint_38
F_faint_39
F_faint_8
F_faint_9""".split("\n")

fwhm_file = "survey"

template = """likelihood:
    star_halo_likelihood.StarHaloLike:
      python_path: /home/idies/workspace/Temporary/maja/scratch/gama09
      stop_at_error: True
      starID: {}
      halo_profile: powerlaw
      fwhm_file: {}
      input_params: [c1, c2]
params:
  c1:
    prior:
      min: 0
      max: 10000
    ref:
      dist: norm
      loc: 100
      scale: 10.0
    proposal: 1
  c2:
    prior:
      min: -100.0
      max: 100.0
    ref:
      dist: norm
      loc: 0.0
      scale: 1.0
    proposal: 0.1
sampler:
  mcmc:
    Rminus1_stop: 0.01
    Rminus1_cl_stop: 0.1
    max_tries: 1000000
output: {}/{}"""

inum = 0
N = len(lae_ids)
for laeID in lae_ids:
	if laeID in flagged:
		continue
	path = "/home/idies/workspace/Temporary/maja/scratch/gama09/cobaya-chains/star_halo/{}".format(laeID)
	try: 
		os.mkdir(path) 
	except OSError as error: 
		print(error)     

	with open(os.path.join(path,"star_profile.yaml"), "w") as f:
		f.write(template.format(laeID, fwhm_file, path, laeID))

	f.close()
	inum += 1
	print(f"Finished {inum}/{N}.")
