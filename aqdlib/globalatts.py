from __future__ import division, print_function
import numpy as np

def read_globalatts(fname):
    a = np.genfromtxt(fname, delimiter=';', dtype=None, autostrip=True)
    metadata = {}
    for row in a:
        if row[0] == 'MOORING':
            metadata[row[0]] = row[1]
        else:
            metadata[row[0]] = str2num(row[1])
    return metadata

def str2num(s):
    try:
        float(s)
        return float(s)
    except ValueError:
        return s
