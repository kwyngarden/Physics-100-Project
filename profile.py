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


# C_i,i+1 are the deprojected galaxy counts
# S_m,m+1 are projected (observed) galaxy counts
# TODO: figure out how to actually use this in the bin densities

# Input: list of galaxy data dicts
# Output: list of galaxy counts per unit volume, one entry per bin
def get_bin_densities(data_list):
    pass # TODO


# Alp ^
# Keith v


# Input: list of bin densities; index of bin for which to get probabilities
# Output: list of index+1 probabilities, where the i'th probability is the
#         probability that a galaxy in bin #bin_index is actually in bin i
def get_bin_probabilities(bin_densities, bin_index):
    bin_densities = bin_densities[:bin_index+1]
    total_density = sum(bin_densities)
    return [bin_density / total_density for bin_density in bin_densities]

# Input: list of galaxy data dicts, list of densities (one entry per bin)
# Output: list of velocity dispersions, one entry per bin
def get_bin_velocity_dispersions(data_list, bin_densities):
    pass

if __name__=='__main__':
    data_list = read_data()

    bin_densities = get_bin_densities(data_list)
    print bin_densities
    bin_velocity_dispersions = get_bin_velocity_dispersions(data_list, bin_densities)
    print bin_velocity_dispersions

    print 'Done'