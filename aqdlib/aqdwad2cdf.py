from __future__ import division, print_function
import numpy as np
import aqdlib

def wad_to_cdf(basefile, metadata):
    """Main waves load file"""

    instmeta = aqdlib.read_aqd_hdr(basefile)

    RAW = {}
    RAW['instmeta'] = instmeta

    RAW = load_whd(RAW, basefile, metadata)

    RAW = load_wad(RAW, basefile, metadata)



    return RAW

def load_whd(RAW, basefile, metadata):
    whdfile = basefile + '.whd'
    print('Loading ' + whdfile)
    WHD = np.genfromtxt(whdfile)

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

    print('Done loading ' + whdfile)

    return RAW

def load_wad(RAW, basefile, metadata):

    wadfile = basefile + '.wad'
    print('Loading wave data from ' + wadfile + '; this may take some time')
    WAD = np.genfromtxt(wadfile)

    r, c = np.shape(WAD)
    print(r, c)
    nburst = np.int(np.floor(r/RAW['instmeta']['WaveNumberOfSamples']))
    nsamps = np.int(nburst * RAW['instmeta']['WaveNumberOfSamples'])
    print(nburst, nsamps)

    RAW['pressure'] = np.reshape(WAD[0:nsamps, 2], (RAW['instmeta']['WaveNumberOfSamples'], nburst))
    RAW['VEL1'] = np.reshape(WAD[0:nsamps, 5], (RAW['instmeta']['WaveNumberOfSamples'], nburst))
    RAW['VEL2'] = np.reshape(WAD[0:nsamps, 6], (RAW['instmeta']['WaveNumberOfSamples'], nburst))
    RAW['VEL3'] = np.reshape(WAD[0:nsamps, 7], (RAW['instmeta']['WaveNumberOfSamples'], nburst))
    RAW['AMP1'] = np.reshape(WAD[0:nsamps, 9], (RAW['instmeta']['WaveNumberOfSamples'], nburst))
    RAW['AMP2'] = np.reshape(WAD[0:nsamps, 10], (RAW['instmeta']['WaveNumberOfSamples'], nburst))
    RAW['AMP3'] = np.reshape(WAD[0:nsamps, 11], (RAW['instmeta']['WaveNumberOfSamples'], nburst))

    print('Done loading ' + wadfile)

    return RAW
# %%
%cd /Volumes/Backstaff/field/gb_proc/1076a/1076a1aqd
metadata = aqdlib.read_globalatts('../glob_att1076a.txt')

# Add additional metadata from metadata config file
import yaml
config = yaml.safe_load(open('config.yaml'))

for k in config:
    metadata[k] = config[k]

RAW = wad_to_cdf('AQ107602', metadata)
# %%
from .aqdhdr2cdf import compute_time
compute_time(asdf)

RAW = compute_time()

# import pickle
# pickle.dump( RAW, open( "RAW.p", "wb" ) )

RAW = pickle.load( open( "RAW.p", "rb" ) )

print(RAW)
