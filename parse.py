from constants import *
import numpy as np
import re
import requests
import sys
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
    if not data[HRV]:
        return None
    if data[HRV_ERR]:
        return float(data[HRV_ERR])

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
            return None
    print '\tERROR: for (%s, %s), found %s HRV error entries.' % (ra_str, dec_str, len(hrv_error_lines))
    return None

def get_lum_and_error(lines, ra_str, dec_str):
    lum_lines = [line for line in lines if '<td>Visual</td>' in line]
    if len(lum_lines) == 1:
        vline = lum_lines[0]
        tokens = [token.strip() for token in re.split('</?t(?:r|d)>', vline) if token != '']
        if len(tokens) != 6:
            print '\tERROR: incorrect formatted luminosity line:\n\t\t%s\n\t\t' % (vline, tokens)
        else:
            lum_token = tokens[-1]
            try:
                if '+/-' in lum_token:
                    lum, lum_err = [float(t.strip()) for t in lum_token.split('+/-')]
                    return lum, lum_err
                else:
                    return float(lum_token.strip()), None
            except:
                e = sys.exc_info()[0]
                print '\tERROR: luminosity extraction hit: %s' % (e)
                return None, None
    else:
        print '\tERROR: for (%s, %s), found %s luminosity entries.' % (ra_str, dec_str, len(lum_lines))
    return None, None

def get_gtype(lines, ra_str, dec_str):
    type_lines = [line for line in lines if 'Galaxy Morphology' in line]
    if len(type_lines) == 1:
        tline = type_lines[0]
        tokens = [token.strip() for token in re.split('</?T(?:R|D)>', tline) if token]
        try:
            return tokens[2]
        except:
            e = sys.exc_info()[0]
            print '\tERROR: galaxy type extraction hit: %s' % (e)
            return None
    else:
        print '\tERROR: for (%s, %s), found %s galaxy type entries.' % (ra_str, dec_str, len(type_lines))
    return None

def get_lum_and_type(data):
    ra_str = data[RA]
    dec_str = data[DEC]
    print 'Looking up lum+gtype for (%s, %s)...' % (ra_str, dec_str)
    sleep(0.65) # Rate-limit
    r = requests.get(get_url(ra_str, dec_str))
    lines = r.content.split('\n')
    lum, lum_err = get_lum_and_error(lines, ra_str, dec_str)
    gtype = get_gtype(lines, ra_str, dec_str)
    return lum, lum_err, gtype

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
            line_data[header] = value if value else None

        line_data[HRV] = float(line_data[HRV]) if line_data[HRV] else None
        line_data[HRV_ERR] = get_hrv_error(line_data)
        
        # lum, lum_err, gtype = get_lum_and_type(line_data)
        # line_data[LUM] = lum
        # line_data[LUM_ERR] = lum_err
        # line_data[GTYPE] = gtype

        line_data[LUM] = float(line_data[LUM]) if line_data[LUM] else None
        line_data[LUM_ERR] = float(line_data[LUM_ERR]) if line_data[LUM_ERR] else None

        line_data[RA] = convert_to_fraction(line_data[RA])
        line_data[DEC] = convert_to_fraction(line_data[DEC])

        mapped_data.append(line_data)
    return mapped_data

if __name__=='__main__':
    data = read_data()
    with open('deleteme.txt', 'w') as f:
        f.write(str(data))
    with open('lum_and_errors.txt', 'w') as f:
        for data_map in data:
            f.write(str(data_map[LUM]) + '\t' + str(data_map[LUM_ERR]) + '\n')
    with open('gtypes.txt', 'w') as f:
        for data_map in data:
            f.write(str(data_map[GTYPE]) + '\n')
    print 'Done'