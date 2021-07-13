from lae_likelihood_all import * 
from scipy.interpolate import interp1d

# before you can run this, you have to comment out the non-chisq part in logp_powerlaw().

lae_like = LAEallLike()

if lae_like.redshift_cut:
	print("This is only for the LAEs with z<2.6.")
	mid_c1, mid_c2 = 2.33810445, 0.12454964
else:
	print("This is for all LAEs.")
	mid_c1, mid_c2 = 1.8320348, -0.060899096#2.4851914, 0.12999488

kwargs = {"sigma_int":0, "c1":mid_c1, "c2":mid_c2}

sigma_ints = np.arange(0.02, 0.5, 0.01)

chisq_wanted = lae_like.NUMBER_OF_DATA_POINTS - 2

chisqs = []                                                                                                                                                                                       
for sigma_int in sigma_ints: 
	kwargs["sigma_int"] = sigma_int 
	logp = lae_like.logp(**kwargs) 
	chisqs.append(-2*logp)
chisqs = np.array(chisqs)
    
print('min, max of chisq/chisq_wanted: ', np.nanmin(chisqs/chisq_wanted), np.nanmax(chisqs/chisq_wanted))

interp = interp1d(chisqs/chisq_wanted, sigma_ints)

print("Wanted sigma_int: ", interp(1))