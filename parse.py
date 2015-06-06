from constants import *
import numpy as np
import requests
from time import sleep


RV_ESTIMATE_STRS = set([
    'V (Kinematic LSR)',
    'V (Galactocentric GSR)',
    'V (Local Group)',
    'V (3K CMB)',
    'V (Virgo Infall only)',
    'V (Virgo + GA only)',
    'V (Virgo + GA + Shapley)',
])


# Given two position strings, e.g. 02 59 32.19, +41 22 33.4, return the NED
# search url.
def get_url(ra_str, dec_str):
    url_ra = '+'.join(ra_str.split())
    url_dec = '+'.join(dec_str.replace('+', '').split())
    return 'https://ned.ipac.caltech.edu/cgi-bin/objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon=%s&lat=%%2B%s&radius=0.1&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ALL&in_objtypes1=Galaxies&search_type=Near+Position+Search&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES' % (
        url_ra, url_dec,
    )

def get_error_or_none(line):
    tokens = line.split()
    if '+/-' not in tokens:
        return None
    error_index = tokens.index('+/-')
    return float(tokens[error_index+1])

def get_max_rv_error(lines):
    matching_lines_errors = [
        get_error_or_none(line) for line in lines
        if any([est_str in line for est_str in RV_ESTIMATE_STRS])
    ]
    return max([rv_err for rv_err in matching_lines_errors if rv_err])

# Given two position strings, e.g. 02 59 32.19, +41 22 33.4, lookup the HRV
# error for the galaxy at that position in the NED database.
def get_hrv_error(data):
    ra_str = data[RA]
    dec_str = data[DEC]
    # Don't look up HRV error, since we already have it or don't need it
    if data[HRV_ERR] != '':
        return float(data[HRV_ERR])
    if data[HRV] == '':
        return ''

    print 'Looking up HRV error for (%s, %s)...' % (ra_str, dec_str)
    sleep(1) # Rate-limit
    r = requests.get(get_url(ra_str, dec_str))
    lines = r.content.split('\n')
    hrv_error_lines = [line for line in lines if ('V (Heliocentric)' in line)]
    if len(hrv_error_lines) == 1:
        ned_hrv_err = get_error_or_none(hrv_error_lines[0])
        if ned_hrv_err:
            return ned_hrv_err
        max_other_err = get_max_rv_error(lines)
        if max_other_err:
            print '\tEstimating error by max error: %s' % (max_other_err)
            return max_other_err
        else:
            print '\tERROR: No error found for (%s, %s)' % (ra_str, dec_str)
            return ''
    print '\tERROR: for (%s, %s), found %s HRV error entries.' % (ra_str, dec_str, len(hrv_error_lines))
    return ''

# Convert measurements like +41 22 33.4 or 02 59 32.19 to decimal numbers.
# Note that the first unit can be degrees or hours, but the other two
# must be minutes and seconds.
def convert_to_fraction(measurement_str):
    whole, minutes, seconds = [float(n) for n in measurement_str.split()]
    return whole + (minutes / 60.) + (seconds / (60.*60.))

def parse_tsv(tsv_file):
    tokenized_lines = []
    with open(tsv_file, 'r') as f:
        for line in f.readlines():
            tokenized_lines.append(line.split('\t'))
    return tokenized_lines[1:]

# Returns list of galaxy data, where each galaxy is a dict of the keys
# listed in constants.py.
def read_data(filename=DATAFILE, headers=HEADERS):
    tokenized_lines = parse_tsv(DATAFILE)
    mapped_data = []
    for line_tokens in tokenized_lines:
        line_data = {}
        for value, header in zip(line_tokens, headers):
            line_data[header] = value

        if line_data[HRV] != '':
            line_data[HRV] = float(line_data[HRV])
        if line_data[MAG] != '':
            line_data[MAG] = float(line_data[MAG])
        line_data[HRV_ERR] = get_hrv_error(line_data)
        
        line_data[RA] = convert_to_fraction(line_data[RA])
        line_data[DEC] = convert_to_fraction(line_data[DEC])

        mapped_data.append(line_data)
    return mapped_data

if __name__=='__main__':
    data = read_data()
    with open('deleteme.txt', 'w') as f:
        f.write(str(data))
    with open('hrv__and_errors.txt', 'w') as f:
        for data_map in data:
            f.write(str(data_map[HRV]) + '\t' + str(data_map[HRV_ERR]) + '\n')
    print 'Done'