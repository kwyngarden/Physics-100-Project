DISTANCE_TO_PERSEUS = 73.6 # Megaparsecs
def mpc_to_cm(mpc):
    return 3.085677581e+24 * mpc

CENTER_RA = 3.3300361111111108
CENTER_DEC = 41.511722222222225

BINS = [(4., 5.2), (3., 4.), (2., 3.), (1.5, 2.), (1., 1.5), (0.5, 1.), (0., 0.5)]
BIN_WIDTHS = [DISTANCE_TO_PERSEUS / (end-start) for (start, end) in BINS]

MIN_RV = 500
MAX_RV = 12500

DATAFILE = 'with_errors.tsv'

NAME = 'Name'
ALT_NAME_1 = 'AltName1'
ALT_NAME_2 = 'AltName2'
RA = 'RA'
DEC = 'Dec'
MAG = 'Mag'
HRV = 'HRV'
HRV_ERR = 'HRV_err'
REMARKS = 'Remarks'

HEADERS = [
    NAME,
    ALT_NAME_1,
    ALT_NAME_2,
    RA,
    DEC,
    MAG,
    HRV,
    HRV_ERR,
    REMARKS
]