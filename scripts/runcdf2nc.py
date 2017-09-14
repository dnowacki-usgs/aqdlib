#!/usr/bin/env python

from __future__ import division, print_function

import sys
sys.path.insert(0, '/Users/dnowacki/Documents/aqdlib')
import aqdlib
import argparse
import yaml

parser = argparse.ArgumentParser(description='Convert raw Aquadopp .cdf format to processed .nc files')
parser.add_argument('cdfname', help='raw .CDF filename')
parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
parser.add_argument('--atmpres', help='path to cdf file containing atmopsheric pressure data')

args = parser.parse_args()

# initialize metadata from the globalatts file
metadata = aqdlib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
config = yaml.safe_load(open(args.config))

for k in config:
    metadata[k] = config[k]

if args.atmpres:
    # press_ac = aqdlib.load_press_ac('press_ac.cdf', ['p_1ac'])
    VEL = aqdlib.cdf_to_nc(args.cdfname, metadata, atmpres=args.atmpres)
else:
    VEL = aqdlib.cdf_to_nc(args.cdfname, metadata)
