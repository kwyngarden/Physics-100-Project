import numpy as np

def mpc_to_m(mpc):
    return 3.08567758128e+22 * mpc
def km_to_m(km):
    return km * 1000.

G = 6.67384e-11 # Gravitational constant

# These are the (start, end) radii of each shell/bin into which we discretized
# the galaxies. All radii are in Megaparsecs. The bin distributions are
# a little misshaped due to some quirks in our data (hopefully going away after
# some bug squashing). 5.2 Mpc is the max radius in our data.
# BINS = [(4.5, 5.2), (3.5, 4.5), (2.75, 3.5), (2., 2.75), (1.5, 2.), (0., 1.5)]
# Alternate, more "normal-sized" bins for experimenting
BINS = [(4., 5.2), (3., 4.), (2., 3.), (1.5, 2.), (1., 1.5), (0.5, 1.), (0., 0.5)]

# Galaxies per cubic meter, calculated over all our data.
AVG_DENSITY = 4.76186343229e-68

def get_bin_densities():
    # Make up a fairly gentle density decline, since our calculated densities
    # seem a little off at the moment.
    slope = AVG_DENSITY / (len(BINS))
    return [AVG_DENSITY + (i-len(BINS)/2)*slope for i in range(len(BINS))]

def get_bin_dispersions():
    # Again, just make up some dispersion numbers for the purpose of this demo,
    # since ours are currently thrown off by the densities. Judging from
    # literature and our data, the velocity dispersion average across the whole
    # cluster is something like ~1200-1300 km/s, but we notice a slight decrease in
    # dispersion moving away from our data, so that is replicated here.
    AVG_DISPERSION = 1300.
    slope = AVG_DISPERSION / (10*len(BINS))
    return [km_to_m(AVG_DISPERSION + (i-len(BINS)/2)*slope) for i in range(len(BINS))]

def jeans_eq_mass_profile(bin_densities, bin_dispersions):
    masses = {}

    for right_bin_index in range(len(BINS)-1, 0, -1):
        # Note that we need a way to get a difference in densities and dispersion
        # across a change in radius. To accomplish this, we consider the radii
        # marking the boundaries between bins, and calculate the changes in
        # these quantities between adjacent bins over the distance between
        # their midpoints.
        r = BINS[right_bin_index][1] # Note r is in Mpc, but converted to meters where used.

        # Note that radius increases as moving from right bins to left bins, so
        # the change in radius is positive.
        left_midpoint = np.mean([BINS[right_bin_index-1][0], BINS[right_bin_index-1][1]])
        right_midpoint = np.mean([BINS[right_bin_index][0], BINS[right_bin_index][1]])
        log_r = np.log(mpc_to_m(left_midpoint)) - np.log(mpc_to_m(right_midpoint))

        # Calculate differences between bins.
        diff_density = np.log(bin_densities[right_bin_index-1]) - np.log(bin_densities[right_bin_index])
        diff_dispersion = np.log(bin_dispersions[right_bin_index-1]**2) - np.log(bin_dispersions[right_bin_index]**2)

        # This is tricky - what's supposed to go here (the velocity dispersion
        # term outside the parentheses in the equation written on that
        # slide; assumed to be velocity dispersion at r)? We take the average
        # of the dispersions in the two adjacent bins.
        dispersion = np.mean([bin_dispersions[right_bin_index-1]**2, bin_dispersions[right_bin_index]**2])

        # Note the derivatives in the Jeans equations are replaced with the differences.
        # This equation gives masses that are too small! Order of 5e39 kg at r = 4.5 Mpc
        masses[r] = -dispersion * mpc_to_m(r) * ((diff_density + diff_dispersion) / log_r) / G

        # Alternative form, using a little calculus to rewrite the derivatives first.
        # But this is far too high of an estimate, ~2e62 at r = 4.5 Mpc
        # masses[r] = -dispersion * mpc_to_m(r) * ((diff_density + diff_dispersion) * mpc_to_m(r) / (left_midpoint - right_midpoint)) / G
    
    return masses

if __name__=='__main__':
    bin_densities = get_bin_densities()
    print 'Using bin densities: %s' % (bin_densities)
    bin_dispersions = get_bin_dispersions()
    print 'Using bin dispersions: %s' % (bin_dispersions)
    masses = jeans_eq_mass_profile(bin_densities, bin_dispersions)
    for r in sorted(masses):
        print '%r Mpc: %s kg' % (r, masses[r])