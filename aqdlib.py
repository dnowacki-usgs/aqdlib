import numpy as np

def str2num(s):
    try:
        float(s)
        return float(s)
    except ValueError:
        return s

def read_globalatts(fname):
    a = np.genfromtxt(fname, delimiter=';', dtype=None, autostrip=True)
    metadata = {}
    for row in a:
        if row[0] == 'MOORING':
            metadata[row[0]] = row[1]
        else:
            metadata[row[0]] = str2num(row[1])
    return metadata
    # return dict(zip(a[:,0], str2num(a[:,1])))


import aqdhdr2cdf
reload(aqdhdr2cdf)

gatts = '/Volumes/Backstaff/field/gb_proc/1076a/glob_att1076a.txt'
metadata = read_globalatts(gatts)

baseFile = '/Volumes/Backstaff/field/gb_proc/1076a/1076a1aqd/AQ107602';

aqdnum = '11819';

metadata['instrument_number'] = '1076a1'
metadata['filename='] = '1076a1aqd'   # name of output file, -raw['cdf or ['nc will be appended to this
metadata['LatLonDatum'] = 'NAD83'
metadata['ClockError'] = 0 # sec, negative is slow
metadata['orientation'] = 'UP'          # use this to identify orientation of profiler
metadata['head_rotation'] = 'horizontal'
metadata['initial_instrument_height'] = 0.15  # meters - estimated!!!
metadata['initial_instrument_height_note'] = ''
metadata['nominal_sensor_depth_note'] = 'WATER_DEPTH-initial_instrument_height'
metadata['nominal_sensor_depth'] = metadata['WATER_DEPTH']-metadata['initial_instrument_height']
metadata['transducer_offset_from_bottom'] = metadata['initial_instrument_height']
#metadata['pred_accuracy'] = 4['8; # horizontal predicted accuracy cm/s                         # degrees between magnetic and true north at data location
metadata['zeroed_pressure'] = 'Yes' # was pressure zeroed before deployment
metadata['cutoff_ampl'] = 0   # set to 0, otherwise it automatically specifies cutoff based on amplitude which is a weird thing to do
metadata['trim_method'] = 'Water Level SL'  # Water Level SL trims bin if any part of bin or side lobe is out of water - works best when pressure is corrected for atmospheric

start_time = '2016-08-03 17:00:00' # first ping
inwater_time = '2016-08-04 15:22:00'
outwater_time = '2016-10-19 20:10:00'
metadata['Deployment_date'] = inwater_time # This needs to be actual deployment date and time on the Bottom
metadata['Recovery_date'] = outwater_time  # This needs to be actual recovery date and time on the Bottom
# %%
import aqdhdr2cdf
reload(aqdhdr2cdf)
RAW = aqdhdr2cdf.load(baseFile,metadata)

# aqdhdr2cdf.readAQDprfHeader(baseFile)

# %%
for k, v in metadata.iteritems():
    print k, ':', v
