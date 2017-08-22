from __future__ import division, print_function
import numpy as np
import datetime as dt
import pytz
from .aqdhdr2cdf import compute_time, read_aqd_hdr, check_metadata, check_orientation, define_aqd_cdf_file, write_aqd_cdf_data
import pandas as pd

def wad_to_cdf(basefile, metadata):
    """Main waves load file"""

    instmeta = read_aqd_hdr(basefile)

    RAW = {}
    RAW['instmeta'] = instmeta

    RAW = load_wad(RAW, basefile, metadata)

    RAW = load_whd(RAW, basefile, metadata)

    RAW = compute_time(RAW, instmeta)

    # Deal with metadata peculiarities
    metadata = check_metadata(metadata, instmeta, waves=True)
    metadata['center_first_bin'] = RAW['cellpos'][0]

    print('BIN SIZE:', metadata['bin_size'])

    RAW = check_orientation(RAW, metadata, waves=True)

    # configure file
    cdf_filename = metadata['filename'] + 'wvs-raw.cdf' # TODO: fix the path
    print('Opening %s' % cdf_filename)

    define_aqd_cdf_file(cdf_filename, RAW, metadata, waves=True)
    print('Variables created')

    write_aqd_cdf_data(cdf_filename, RAW, metadata, waves=True)
    print('Variables written')

    return RAW, metadata

def load_whd(RAW, basefile, metadata):
    whdfile = basefile + '.whd'
    print('Loading ' + whdfile)
    WHD = np.loadtxt(whdfile)

    RAW['burst'] = WHD[:,6]
    RAW['nrecs'] = WHD[:,7]
    RAW['cellpos'] = WHD[:,8]
    RAW['battery'] = WHD[:,9]
    RAW['soundspeed'] = WHD[:,10]
    RAW['heading'] = WHD[:,11]
    RAW['pitch'] = WHD[:,12]
    RAW['roll'] = WHD[:,13]
    RAW['minpressure'] = WHD[:,14]
    RAW['temperature'] = WHD[:,16]
    RAW['cellsize'] = WHD[:,17]
    RAW['avgamp1'] = WHD[:,18]
    RAW['avgamp2'] = WHD[:,19]
    RAW['avgamp3'] = WHD[:,20]

    RAW['datetime'] = []
    for year, month, day, hour, minute, second in zip(WHD[:,2], WHD[:,0], WHD[:,1], WHD[:,3], WHD[:,4], WHD[:,5]):
        RAW['datetime'].append(dt.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), tzinfo=pytz.utc))
    RAW['datetime'] = np.array(RAW['datetime'])

    print('Done loading ' + whdfile)

    return RAW

def load_wad(RAW, basefile, metadata):

    wadfile = basefile + '.wad'
    print('Loading wave data from ' + wadfile + '; this may take some time')
    # pd.read_csv is ~10x faster than np.loadtxt or np.genfromtxt
    WAD = pd.read_csv(wadfile, header=None, delim_whitespace=True).values
    # WAD = np.loadtxt(wadfile)

    r, c = np.shape(WAD)
    print(r, c)
    nburst = int(np.floor(r/RAW['instmeta']['WaveNumberOfSamples']))
    nsamps = int(nburst * RAW['instmeta']['WaveNumberOfSamples'])
    wavensamps = int(RAW['instmeta']['WaveNumberOfSamples'])
    print(nburst, nsamps, wavensamps)

    RAW['pressure'] = np.reshape(WAD[0:nsamps, 2], (nburst, wavensamps)).T
    RAW['VEL1'] = np.reshape(WAD[0:nsamps, 5], (nburst, wavensamps)).T
    RAW['VEL2'] = np.reshape(WAD[0:nsamps, 6], (nburst, wavensamps)).T
    RAW['VEL3'] = np.reshape(WAD[0:nsamps, 7], (nburst, wavensamps)).T
    RAW['AMP1'] = np.reshape(WAD[0:nsamps, 9], (nburst, wavensamps)).T
    RAW['AMP2'] = np.reshape(WAD[0:nsamps, 10], (nburst, wavensamps)).T
    RAW['AMP3'] = np.reshape(WAD[0:nsamps, 11], (nburst, wavensamps)).T

    # convert to cm/s
    for n in [1, 2, 3]:
        RAW['VEL' + str(n)] = RAW['VEL' + str(n)] * 100

    print('Done loading ' + wadfile)

    return RAW

# %%
