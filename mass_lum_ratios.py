from constants import *
from parse import read_data
import re

DEFAULT_RATIO = (10., 2.1)

M_LUM_RATIOS = {
    ('s', 's0', 's00', 'e', 'e0', 'e1', 'e2', 'e3', 'e4', 'e5'): (10., 2.1),
    ('sa', 'sa0', 'sa00', 's0a'): (6.2, 1.1),
    ('sab', 'sabc', 'sabb', 'sabbc', 'sabd', 'sam', 'sb', 'sbab', 'sbb', 'sba', 'sb0', 'sb0a'): (6.5, 0.5),
    ('sbc', 'sbbc', 'sb0', 'sc'): (4.7, 0.4),
    ('scd', 'sd', 'cd'): (3.9, 0.6),
    ('sdm', 'irr'): (1.7, 0.6)
}

def clean_type(gtype):
    cleaned = gtype.lower().replace('?', '').replace('(r)', '').replace('(s)', '').replace('(r\')', '').replace('(rs)', '')
    cleaned = cleaned.replace('edge-on', '').replace('pec', '')
    cleaned = re.sub('[/^+-]', '', cleaned)
    return cleaned.strip()

def get_mass_estimate(luminosity, gtype):
    cleaned_type = clean_type(gtype)
    for typeset in M_LUM_RATIOS:
        if cleaned_type in typeset:
            lum_ratio, err = M_LUM_RATIOS[typeset]
            return SOLAR_MASS * lum_ratio * luminosity, SOLAR_MASS * err * luminosity
    return SOLAR_MASS * DEFAULT_RATIO[0] * luminosity, SOLAR_MASS * DEFAULT_RATIO[1] * luminosity

if __name__=='__main__':
    data_list = read_data()
    gtypes = set()
    for data in data_list:
        if data[LUM] and data[GTYPE]:
            gtypes.add(clean_type(data[GTYPE]))
    for gtype in gtypes:
        if not any([gtype in k for k in M_LUM_RATIOS]):
            print gtype