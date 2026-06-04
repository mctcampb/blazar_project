## This script is not generalizable
from astropy.io import fits 
import numpy as np 
import pycorr 

theta_bins = np.linspace(5, 40, num = 11) # Max bin estimated from rough size of patches 

rands = fits.open('randoms.fits')[1].data
ra_rands = rands['ra']
dec_rands = rands['dec']

def compute_wtheta(sample1, njk, output_fname, sample2 = None):
    if sample2 == None:
        objs1 = fits.open(sample1)[1].data
        ra_objs1 = objs1['ra']
        dec_objs1 = objs1['dec']

        combined_ra = np.concatenate((ra_objs1, ra_rands))
        combined_dec = np.concatenate((dec_objs1, dec_rands))

        subsampler = pycorr.twopoint_jackknife.KMeansSubsampler('angular', [combined_ra, combined_dec], nsamples = njk, position_type = 'rd', random_state = 1, nside = 64)
        labels = subsampler.label([combined_ra, combined_dec], position_type = 'rd')
        labels_objs1 = labels[:len(ra_objs1)]
        labels_rands = labels[len(ra_objs1):]
        res = pycorr.TwoPointCorrelationFunction('theta', theta_bins, [ra_objs1, dec_objs1], randoms_positions1 = [ra_rands, dec_rands], estimator = 'landyszalay', bin_type = 'custom', position_type = 'rd', compute_sepsavg = False, nthreads = 1, data_samples1 = labels_objs1, randoms_samples1 = labels_rands)
        np.save(output_fname + '.npy', res, allow_pickle = True)
        
    else:
        objs1 = fits.open(sample1)[1].data
        ra_objs1 = objs1['ra']
        dec_objs1 = objs1['dec']

        objs2 = fits.open(sample2)[1].data
        ra_objs2 = objs2['ra']
        dec_objs2 = objs2['dec']

        combined_ra = np.concatenate((ra_objs1, ra_objs2, ra_rands))
        combined_dec = np.concatenate((dec_objs1, dec_objs2, dec_rands))

        subsampler = pycorr.twopoint_jackknife.KMeansSubsampler('angular', [combined_ra, combined_dec], nsamples = njk, position_type = 'rd', random_state = 1, nside = 64)
        labels = subsampler.label([combined_ra, combined_dec], position_type = 'rd')
        labels_objs1 = labels[:len(ra_objs1)]
        labels_objs2 = labels[len(ra_objs1):len(ra_objs1) + len(ra_objs2)]
        labels_rands = labels[len(ra_objs1) + len(ra_objs2):]
        res = pycorr.TwoPointCorrelationFunction('theta', theta_bins, [ra_objs1, dec_objs1], data_positions2 = [ra_objs2, dec_objs2], randoms_positions1 = [ra_rands, dec_rands], randoms_positions2 = [ra_rands, dec_rands],  estimator = 'landyszalay', bin_type = 'custom', position_type = 'rd', compute_sepsavg = False, nthreads = 1, data_samples1 = labels_objs1, data_samples2 = labels_objs2, randoms_samples1 = labels_rands, randoms_samples2 = labels_rands)
        np.save(output_fname + '.npy', res, allow_pickle = True)
