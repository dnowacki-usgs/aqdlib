#!/usr/bin/env python

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime as dt
import sys
import os
sys.path.append('/Users/dnowacki/Documents/python')
sys.path.append('/Users/dnowacki/Documents/aqdlib')
sys.path.append('/Users/dnowacki/Documents/Grand Bay/py')
import aqdlib
import gbts
import shutil
import argparse
import yaml

parser = argparse.ArgumentParser(description='Convert Aquadopp text files to raw .cdf format. Run this script from the directory containing Aquadopp files')
parser.add_argument('basename', help='base name (without extension) of the Aquadopp text files')
parser.add_argument('gatts', help='path to global attributes file (gatts formatted)')
parser.add_argument('metadata', help='path to ancillary metadata file (YAML formatted)')
# parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')

args = parser.parse_args()

# initialize metadata from the globalatts file
metadata = aqdlib.read_globalatts(args.gatts)

# Add additional metadata from metadata config file
config = yaml.safe_load(open(args.metadata))

for k in config:
    metadata[k] = config[k]

RAW = aqdlib.prf_to_cdf(args.basename, metadata)
