## Some preamble... after reading the DESI manual, it seems like pairs of points below 0.05 degrees may produce biased statistics due to close pairs. So, we should avoid pairs closer than 0.05 * pi / 180 * 1825 h^-1 Mpc (distance is estimated from the mean survey redshift) = 1.6 h^-1 Mpc. Maybe to be safe, we will make a scale cut above 5 Mpc h^-1.  
from astropy.cosmology import Planck18
import astropy.cosmology.units as cu
from astropy.io import fits
import astropy.units as u
import numpy as np
import pycorr
print('done loading packages')

rp_bins = np.logspace(np.log10(5.), np.log10(80.), num = 21) # 5 to 80 h^-1 Mpc 
pi_bins = np.linspace(-100., 100., 2000) # Ensured that pi_max > largest rp bin 
njk = 100 # Checked that this is a safe number. LRG area = 5740 sq degrees
nside_jk = 4096

gals = fits.open('DESI_DR1_LRG.fits')[1].data
rands = fits.open('DESI_DR1_randoms.fits')[1].data

ra_gals = gals['ra']
dec_gals = gals['dec']
z_gals = gals['z'] * cu.redshift # No tomography for this fiducial run 
d_gals = z_gals.to(u.Mpc, cu.redshift_distance(Planck18, kind = "comoving")).value * Planck18.h 
w_gals = gals['weight']

n_rands = len(ra_gals) * 40 
ra_rands = rands['ra']

if len(ra_rands) > n_rands:
    indices = np.random.choice(len(ra_rands), size = n_rands, replace = False) # Randomly downsample randoms if there are a lot
    ra_rands = ra_rands[indices]
    dec_rands = rands['dec'][indices]
    z_rands = rands['z'][indices] * cu.redshift   
    w_rands = rands['weight'][indices]

else:
    dec_rands = rands['dec']
    z_rands = rands['z'] * cu.redshift  
    w_rands = rands['weight']

d_rands = z_rands.to(u.Mpc, cu.redshift_distance(Planck18, kind = "comoving")).value * Planck18.h
combined_ra = np.concatenate((ra_gals, ra_rands))
combined_dec = np.concatenate((dec_gals, dec_rands))
print('creating randoms..', flush = True)
subsampler = pycorr.twopoint_jackknife.KMeansSubsampler("angular", [combined_ra, combined_dec], nsamples = njk, position_type = "rd", random_state = 1, nside = nside_jk)
labels = subsampler.label([combined_ra, combined_dec], position_type = "rd")
labels_gals = labels[:len(ra_gals)]
labels_rands = labels[len(ra_gals):]

print('computing wp..', flush = True)
eta = pycorr.TwoPointCorrelationFunction("rppi", (rp_bins, pi_bins), [ra_gals, dec_gals, d_gals], randoms_positions1 = [ra_rands, dec_rands, d_rands], data_weights1 = w_gals, randoms_weights1 = w_rands, data_samples1 = labels_gals, randoms_samples1 = labels_rands, estimator = "landyszalay", bin_type = "custom", position_type = "rdd", los = "midpoint", compute_sepsavg = False, nthreads = 8)
np.save('wp_desi.npy', eta, allow_pickle = True)