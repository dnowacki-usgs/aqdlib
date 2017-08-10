#!/usr/bin/env python

from __future__ import division, print_function

import aqdlib
import argparse
import yaml

parser = argparse.ArgumentParser(description='Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files')
parser.add_argument('basename', help='base name (without extension) of the Aquadopp text files')
parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
parser.add_argument('config', help='path to ancillary config file (YAML formatted)')
# parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')

args = parser.parse_args()

# initialize metadata from the globalatts file
metadata = aqdlib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
config = yaml.safe_load(open(args.config))

for k in config:
    metadata[k] = config[k]

RAW = aqdlib.prf_to_cdf(args.basename, metadata)
