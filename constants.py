DISTANCE_TO_PERSEUS = 73.6 # Megaparsecs
def mpc_to_cm(mpc):
    return 3.08567758128e+24 * mpc
def mpc_to_m(mpc):
    return 3.08567758128e+22 * mpc
def km_to_m(km):
    return km * 1000.
def cm3_to_m3_density(density):
    # Convert a density in cm^-3 to a density in m^-3
    return density * (100. ** 3)

CENTER_RA = 3.3300361111111108
CENTER_DEC = 41.511722222222225

# BINS = [(4., 5.2), (3., 4.), (2., 3.), (1.5, 2.), (1., 1.5), (0.5, 1.), (0., 0.5)]
#BINS = [(4.5, 5.2), (3.5, 4.5), (2.75, 3.5), (2., 2.75), (1.5, 2.), (0., 1.5)]
BINS = [(0., 5.2)]
BIN_WIDTHS = [DISTANCE_TO_PERSEUS / (end-start) for (start, end) in BINS]

# Used for graphing purposes
BIN_CENTERS = list(reversed([(a_bin[0] + a_bin[1]) / 2 for a_bin in BINS]))
BIN_ERRORS = list(reversed([(a_bin[1] - a_bin[0]) / 2 for a_bin in BINS]))


MIN_RV = 500
MAX_RV = 12500

G = 6.67384e-11
SOLAR_MASS = 1.9891e30

DATAFILE = 'with_errors.tsv'

NAME = 'Name'
ALT_NAME_1 = 'AltName1'
ALT_NAME_2 = 'AltName2'
RA = 'RA'
DEC = 'Dec'
SUBR = 'SuBr'
LUM = 'Lum'
LUM_ERR = 'Lum_err'
HRV = 'HRV'
HRV_ERR = 'HRV_err'
GTYPE = 'Galaxy_type'
REMARKS = 'Remarks'

HEADERS = [
    NAME,
    ALT_NAME_1,
    ALT_NAME_2,
    RA,
    DEC,
    SUBR,
    LUM,
    LUM_ERR,
    HRV,
    HRV_ERR,
    GTYPE,
    REMARKS,
]
