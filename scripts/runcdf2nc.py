#!/usr/bin/env python

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import sys
import os
sys.path.append('/Users/dnowacki/Documents/aqdlib')
import aqdlib
import shutil
import argparse
import yaml

parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf format to processed .nc files')
parser.add_argument('cdfname', help='raw .CDF filename')
parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
parser.add_argument('metadata', help='path to ancillary metadata file (YAML formatted)')
parser.add_argument('--p_1ac', help='path to cdf file containing atmospherically corrected pressure data')

args = parser.parse_args()

# initialize metadata from the globalatts file
metadata = aqdlib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
config = yaml.safe_load(open(args.metadata))

for k in config:
    metadata[k] = config[k]

# add a few extra metadata variables
metadata['nominal_sensor_depth_note']: 'WATER_DEPTH-initial_instrument_height'
metadata['nominal_sensor_depth'] = metadata['WATER_DEPTH'] - metadata['initial_instrument_height']
metadata['transducer_offset_from_bottom'] = metadata['initial_instrument_height']

if args.p_1ac:
    press_ac = aqdlib.load_press_ac('press_ac.cdf', ['p_1ac'])
    VEL = aqdlib.cdf_to_nc(args.cdfname, metadata, p_1ac=press_ac['p_1ac'])
else:
    VEL = aqdlib.cdf_to_nc(args.cdfname, metadata)
