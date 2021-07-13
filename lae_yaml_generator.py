import os
import sys

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
gama09F_97_318
gama09F_98_436""".split("\n")

template = """likelihood:
    lae_likelihood.LAELike:
      python_path: /home/idies/workspace/Temporary/maja/scratch/gama09
      stop_at_error: True
      laeID: {}
      simulated: False
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

troughsub = 'True'
fwhm_file = "survey"

inum = 0
N = len(lae_ids)
for laeID in lae_ids:
	path = "/home/idies/workspace/Temporary/maja/scratch/gama09/cobaya-chains/LAEs_individual/{}".format(laeID)
	try: 
		os.mkdir(path) 
	except OSError as error: 
		print(error)     

	with open(os.path.join(path,"lae_profile.yaml"), "w") as f:
		f.write(template.format(laeID, troughsub, fwhm_file, path, laeID))

	f.close()
	inum += 1
	print(f"Finished {inum}/{N}.")
