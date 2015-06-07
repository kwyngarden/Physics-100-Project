DISTANCE_TO_PERSEUS = 73.6

CENTER_RA = 3.3300361111111108
CENTER_DEC = 41.511722222222225

BINS = [(4., 5.5), (3., 4.), (2., 3.), (1., 2.), (0., 1.)]
BIN_WIDTH = DISTANCE_TO_PERSEUS / 1.

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