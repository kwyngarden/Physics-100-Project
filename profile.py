from collections import Counter
from constants import *
import numpy as np
from parse import read_data
import matplotlib.pyplot as plt


def convert_to_cartesian(ra, dec, dist=DISTANCE_TO_PERSEUS):
    # http://www.neoprogrammics.com/distance_between_two_stars/
    # http://www.projectrho.com/public_html/starmaps/trigonometry.php
    alpha = np.radians(ra * (360./24.))
    delta = np.radians(dec)
    x = dist * np.cos(alpha) * np.cos(delta)
    y = dist * np.sin(alpha) * np.cos(delta)
    z = dist * np.sin(delta)
    return x, y, z

# Input: right ascensions and declinations in decimal hours and degrees
# Output: apparent distance between points to observer inside sphere
def get_apparent_sphere_distance(ra1, dec1, ra2, dec2):
    x1, y1, z1 = convert_to_cartesian(ra1, dec1)
    x2, y2, z2 = convert_to_cartesian(ra2, dec2)
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

def plot_radial_velocities(data_list):
    dist_list = [get_apparent_sphere_distance(data[RA], data[DEC], CENTER_RA, CENTER_DEC) for data in data_list]
    hrvs = [data[HRV] for data in data_list]
    hrv_errs = [data[HRV_ERR] for data in data_list]
    xs, ys, xerrs, yerrs = [], [], [], []
    for dist, hrv, hrv_err in zip(dist_list, hrvs, hrv_errs):
        if hrv and hrv_err:
            xs.append(dist)
            ys.append(hrv)
            xerrs.append(0)
            yerrs.append(hrv_err)
    plt.errorbar(xs, ys, xerr=xerrs, yerr=yerrs, fmt='o')
    # plt.hist(dist_list, bins=40)
    plt.show()


# Input: list of distances of galaxies from center to cluster
# Output: map from bin index (0 -> len(BINS)-1) to counts of galaxies observed in that bin
def get_observed_densities_per_bin(dist_list):
    counts = Counter()
    for dist in dist_list:
        for i in range(len(BINS)):
            start, end = BINS[i]
            if dist >= start:
                counts[i] += 1
    return counts

# Input: list of galaxy data dicts
# Output: list of galaxy counts per unit volume, one entry per bin
def get_bin_densities(data_list):
    # Image in fb chat:
    # C_i,i+1 are the deprojected galaxy counts
    # S_m,m+1 are projected (observed) galaxy counts
    dist_list = [get_apparent_sphere_distance(data[RA], data[DEC], CENTER_RA, CENTER_DEC) for data in data_list]
    observed_densities = get_observed_densities_per_bin(dist_list)
    print 'Observed densities: %s' % (observed_densities)
    num_bins = len(BINS)
    C_list = list()
    A_list = [get_A(i, i + 1) for i in range(num_bins - 1)]
    b = ((np.pi * mpc_to_cm(DISTANCE_TO_PERSEUS)) / 10800)**2
    for i in range(num_bins - 1):
        the_sum = observed_densities[i] / (b / A_list[i])
        partial_sum = sum([sub_sum(j, i, C_list) for j in range(i)])
        C_i = (the_sum - partial_sum) / get_V(i, i + 1)
        C_list.append(C_i)
    return C_list

def get_A(i, j):
    return np.pi * (get_r(i)**2 - get_r(j)**2)

def sub_sum(i, m, C_list):
    return C_list[i] * (get_V(i, m + 1) - get_V(i + 1, m + 1) - get_V(i, m) + get_V(i + 1, m))
    

def get_V(i, j):
    if i >= j:
        return 0.0
    return (4 / 3) * np.pi * (get_r(i)**2 - get_r(j)**2)**1.5

def get_r(i):
    return BINS[i][1]    

# Alp ^
# Keith v


# Input: list of bin densities; index of bin for which to get probabilities
# Output: list of index+1 probabilities, where the i'th probability is the
#         probability that a galaxy in bin #bin_index is actually in bin i
def get_bin_probabilities(bin_densities, bin_index):
    bin_densities = bin_densities[:bin_index+1]
    total_density = sum(bin_densities)
    return [bin_density / total_density for bin_density in bin_densities]

# Input: list of bin densities
# Output: map of bin index -> list of source bin probabilities
def get_bin_probability_map(bin_densities):
    bin_prob_map = {}
    for i in range(len(BINS)):
        bin_prob_map[i] = get_bin_probabilities(bin_densities, i)
    return bin_prob_map

# Input: list of galaxy data
# Output: map from bin index to list of galaxy data in that bin
def get_bin_map(data_list):
    bin_map = {}
    for data in data_list:
        dist = get_apparent_sphere_distance(data[RA], data[DEC], CENTER_RA, CENTER_DEC)
        for i in range(len(BINS)):
            start, end = BINS[i]
            if dist >= start:
                if i not in bin_map:
                    bin_map[i] = []
                bin_map[i].append(data)
    return bin_map

# Input: list of galaxy data dicts (ones with RV only), list of densities (one entry per bin)
# Output: list of velocity dispersions, one entry per bin
def get_bin_velocity_dispersions(data_list, bin_densities):
    bin_map = get_bin_map(data_list)
    bin_prob_map = get_bin_probability_map(bin_densities)
    means = []
    for bin_index in range(len(BINS)):
        observed_mean_rv = np.mean([data[HRV] for data in bin_map[bin_index]])
        own_bin_prob = bin_prob_map[bin_index][bin_index]
        prev_bin_rv_sum = 0.
        for prev_bin_index in range(bin_index):
            prev_bin_rv_sum += bin_prob_map[bin_index][prev_bin_index] * means[prev_bin_index]
        means.append((observed_mean_rv - prev_bin_rv_sum) / own_bin_prob)

    variances = []
    for bin_index in range(len(BINS)):
        observed_var_rv = np.var([data[HRV] for data in bin_map[bin_index]], ddof=1)
        own_bin_prob = bin_prob_map[bin_index][bin_index]
        prev_bin_weighted_vars = 0.
        weighted_sq_means = 0.
        weighted_means_sq = 0.

        for prev_bin_index in range(bin_index+1):
            prev_bin_prob = bin_prob_map[bin_index][prev_bin_index]
            if prev_bin_index < bin_index:
                prev_bin_weighted_vars += prev_bin_prob * variances[prev_bin_index]
            weighted_sq_means += prev_bin_prob * (means[prev_bin_index] ** 2)
            weighted_means_sq += prev_bin_prob * means[prev_bin_index]
        weighted_means_sq = weighted_means_sq ** 2

        variances.append((observed_var_rv - prev_bin_weighted_vars - weighted_sq_means + weighted_means_sq) / own_bin_prob)

    dispersions = [np.sqrt(var) for var in variances]
    return dispersions

# Input: list of galaxy data dicts (ones with RV only), list of densities (one entry per bin)
# Output: list of velocity dispersions (one per bin), list of errors on those dispersions
#         obtained via Monte Carlo simulation
def get_bin_dispersions_and_errors(data_list, bin_densities, iters=2000, print_progress=True):
    bin_dispersion_lists = {i: [] for i in range(len(BINS))}
    for iter_num in xrange(iters):
        if print_progress and iter_num % 500 == 0 and iter_num > 0:
            print 'Iteration %s' % (iter_num)
        fake_data_list = [
            {
                HRV: np.random.normal(data[HRV], data[HRV_ERR]),
                RA: data[RA],
                DEC: data[DEC],
            } for data in data_list
        ]
        dispersions = get_bin_velocity_dispersions(fake_data_list, bin_densities)
        for i in range(len(BINS)):
            bin_dispersion_lists[i].append(dispersions[i])

    bin_dispersions = []
    bin_dispersion_errs = []
    for bin_index in range(len(BINS)):
        bin_dispersion_list = bin_dispersion_lists[bin_index]
        bin_dispersions.append(np.mean(bin_dispersion_list))
        bin_dispersion_errs.append(np.std(bin_dispersion_list) / np.sqrt(len(bin_dispersion_list)))
    return bin_dispersions, bin_dispersion_errs

if __name__=='__main__':
    data_list = read_data()

    bin_densities = get_bin_densities(data_list)
    # bin_densities = [0.005, 0.01, 0.025, 0.06, 0.11, 0.175, 0.29]
    print 'Using bin densities: %s' % (bin_densities)

    data_with_rv = [data for data in data_list if data[HRV] and data[HRV]>MIN_RV and data[HRV]<MAX_RV]
    bin_dispersions, bin_dispersion_errs = get_bin_dispersions_and_errors(data_with_rv, bin_densities)
    for i in range(len(BINS)):
        print 'Bin %s: %.2f km/s (+/- %.2f)' % (i, bin_dispersions[i], bin_dispersion_errs[i])

    print 'Done'
