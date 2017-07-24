from __future__ import division
from netCDF4 import Dataset
import numpy as np
import qaqc

def cdf_to_nc(cdf_filename, metadata, p1_ac=False):

    nc_filename = metadata['filename'] + '.nc'

    rg = Dataset(cdf_filename, 'r')

    VEL = {}

    VEL = load_cdf_amp_vel(rg, VEL)

    define_aqd_nc_file(nc_filename, VEL, metadata)

    write_aqd_nc_file(nc_filename, VEL, metadata)

    rg.close()

    return VEL

def load_cdf_amp_vel(rg, VEL):
    # TODO: need to clip start and end properly

    vel1 = rg['VEL1'][:]
    vel2 = rg['VEL2'][:]
    vel3 = rg['VEL3'][:]

    amp1 = rg['AMP1'][:]
    amp2 = rg['AMP2'][:]
    amp3 = rg['AMP3'][:]

    # initialize arrays
    VEL['U'] = np.zeros(np.shape(vel1))
    VEL['V'] = np.zeros(np.shape(vel2))
    VEL['W'] = np.zeros(np.shape(vel3))

    heading = rg['Heading'][:]
    pitch = rg['Pitch'][:]
    roll = rg['Roll'][:]

    T = rg['TransMatrix'][:]

    N, M = np.shape(vel1)

    VEL['U'], VEL['V'], VEL['W'] = qaqc.coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, rg.AQDCoordinateSystem)

    VEL['pressure'] = rg['Pressure'][:]
    VEL['temp'] = rg['Temperature'][:]
    VEL['time'] = rg['time'][:]
    VEL['time2'] = rg['time2'][:]
    VEL['AGC'] = (amp1 + amp2 + amp3) / 3

    return VEL

def define_aqd_nc_file(nc_filename, VEL, metadata):
    doublefill = 1e35

    N, M = np.shape(VEL['U'])

    double_fill_val = 1e35;

    rg = Dataset('/Volumes/Backstaff/field/gb_proc/1076a/1076a1aqd/' + nc_filename, 'w', format='NETCDF4', clobber=True)

    time = rg.createDimension('time', 0)
    depth = rg.createDimension('depth', M)
    lat = rg.createDimension('lat', 1)
    lon = rg.createDimension('lon', 1)
    Tmatrix = rg.createDimension('Tmatrix', 3)

    timeid = rg.createVariable('time', 'i', ('time',), zlib=True) # 'i' == NC_INT
    timeid.units = 'True Julian Day'
    timeid.type = 'EVEN' # TODO: this is "UNEVEN" in the CDF version
    timeid.epic_code = 624

    time2id = rg.createVariable('time2', 'i', ('time',), zlib=True) # 'i' == NC_INT
    time2id.units = 'msec since 0:00 GMT'
    time2id.type ='EVEN'
    time2id.epic_code = 624

    latid = rg.createVariable('lat', 'f', ('lat',), zlib=True)
    latid.units = 'degree_north'
    latid.epic_code = 500

    lonid = rg.createVariable('lon', 'f', ('lon',), zlib=True)
    lonid.FORTRAN_format = 'F10.4' # TODO: Why is there a FORTRAN_format on here but not the others?????
    lonid.units = 'degree_east'
    lonid.epic_code = 502

    depthid = rg.createVariable('depth', 'f', ('depth',), zlib=True)
    depthid.units = 'm'
    depthid.epic_code = 3
    depthid.long_name = 'mean water depth'
    depthid.initial_instrument_height = metadata['initial_instrument_height']
    # depthid.nominal_instrument_depth = metadata['nominal_instrument_depth'] # FIXME

    bindistid = rg.createVariable('bindist', 'f', ('depth',))
    # TODO: Loop through the attribute names and values as done in Matlab
    bindistid.blanking_distance = metadata['blanking_distance'] # TODO: the Matlab one uses INFO
    bindistid.initial_instrument_height = metadata['initial_instrument_height']
    bindistid.note = 'distance is along profile from instrument head to center of bin'

    u_1205 = rg.createVariable('u_1205', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
    # u_1205.name = 'u'
    u_1205.setncattr('name', 'u')
    u_1205.long_name = 'Eastward Velocity'
    u_1205.generic_name = 'u'
    u_1205.units = 'cm/s'
    u_1205.epic_code = 1205
    u_1205.sensor_type = metadata['INST_TYPE']
    u_1205.initial_instrument_height = metadata['initial_instrument_height']
    # u_1205.nominal_instrument_depth = metadata['nominal_instrument_depth'] # FIXME
    u_1205.height_depth_units = 'm'
    # u_1205.serial_number = metadata['serial_number'] # FIXME
    u_1205.maximum = 0
    u_1205.minimum = 0
    # TODO: trim_method

    v_1206 = rg.createVariable('v_1206', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
    v_1206.setncattr('name', 'v')
    v_1206.long_name = 'Northward Velocity'
    v_1206.generic_name = 'v'
    v_1206.units = 'cm/s'
    v_1206.epic_code = 1206
    v_1206.sensor_type = metadata['INST_TYPE']
    v_1206.initial_instrument_height = metadata['initial_instrument_height']
    # v_1206.nominal_instrument_depth = metadata['nominal_instrument_depth'] # FIXME
    v_1206.height_depth_units = 'm'
    # v_1206.serial_number = metadata['serial_number'] # FIXME
    v_1206.maximum = 0
    v_1206.minimum = 0
    # TODO: trim_method

    w_1204 = rg.createVariable('w_1204', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
    w_1204.setncattr('name', 'w')
    w_1204.long_name = 'Vertical Velocity'
    w_1204.generic_name = 'w'
    w_1204.units = 'cm/s'
    w_1204.epic_code = 1204
    w_1204.sensor_type = metadata['INST_TYPE']
    w_1204.initial_instrument_height = metadata['initial_instrument_height']
    # w_1204.nominal_instrument_depth = metadata['nominal_instrument_depth'] # FIXME
    w_1204.height_depth_units = 'm'
    # w_1204.serial_number = metadata['serial_number'] # FIXME
    w_1204.maximum = 0
    w_1204.minimum = 0
    # TODO: trim_method

    # if isfield(metadata,'trim_method') & strcmp(metadata.trim_method,'Water Level SL')
    #     netcdf.putAtt(ncid,Uid,'note','Velocity bins trimmed if out of water or if side lobes intersect sea surface');
    # end

    Pressid = rg.createVariable('P_1', 'f', ('lat', 'lon', 'time',), zlib=True)
    Pressid.units = 'dbar'
    Pressid.epic_code = 1
    Pressid.setncattr('name', 'P')
    Pressid.long_name = 'PRESSURE (DB)          ' # TODO: why so many extra spaces?
    Pressid.generic_name = 'depth'
    Pressid.minimum = doublefill
    Pressid.maximum = doublefill
    # netcdf.putAtt(ncid,Pressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); #FIXME
    Pressid.initial_instrument_height = metadata['initial_instrument_height']
    # netcdf.putAtt(ncid,Pressid,'nominal_instrument_depth',metadata.nominal_instrument_depth); #FIXME

    Tempid = rg.createVariable('Tx_1211', 'f', ('lat', 'lon', 'time',), zlib=True)
    Tempid.units = 'C'
    Tempid.epic_code = 1211
    Tempid.setncattr('name', 'Tx')
    Tempid.long_name = 'Instrument Transducer Temperature'
    Tempid.generic_name = 'temp'
    Tempid.minimum = doublefill
    Tempid.maximum = doublefill
    # netcdf.putAtt(ncid,Pressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); #FIXME
    Tempid.initial_instrument_height = metadata['initial_instrument_height']
    # netcdf.putAtt(ncid,Pressid,'nominal_instrument_depth',metadata.nominal_instrument_depth); #FIXME

    AGCid = rg.createVariable('AGC_1202', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
    AGCid.units = 'counts'
    AGCid.epic_code = 1202
    AGCid.setncattr('name', 'AGC')
    AGCid.long_name = 'Average Echo Intensity (AGC)'
    AGCid.generic_name = 'AGC'
    # AGCid.sensor_type = INST_TYPE # FIXME
    AGCid.minimum = 0
    AGCid.maximum = 0 # TODO: why are min/max different (0 vs doublefill for others?)
    # netcdf.putAtt(ncid,AGCid,'sensor_type',metadata.cdfmeta.INST_TYPE); #FIXME
    # netcdf.putAtt(ncid,Pressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); #FIXME
    AGCid.initial_instrument_height = metadata['initial_instrument_height']
    AGCid.height_depth_units = 'm'
    # netcdf.putAtt(ncid,Pressid,'nominal_instrument_depth',metadata.nominal_instrument_depth); #FIXME

    rg.close()

def write_aqd_nc_file(nc_filename, VEL, metadata):
    N, M = np.shape(VEL['U'])

    rg = Dataset('/Volumes/Backstaff/field/gb_proc/1076a/1076a1aqd/' + nc_filename, 'r+')

    rg['lat'][:] = metadata['latitude']
    rg['lon'][:] = metadata['longitude']
    # rg['time'][:] = VEL['time']
    # rg['time2'][:] = VEL['time2']

    # rg['bindist'] = VEL['bindist']
    rg['u_1205'][:] = np.reshape(VEL['U'].T, (M, 1, 1, N))
    rg['v_1206'][:] = np.reshape(VEL['V'].T, (M, 1, 1, N))
    rg['w_1204'][:] = np.reshape(VEL['W'].T, (M, 1, 1, N))

    rg['P_1'][:] = VEL['pressure'][np.newaxis, np.newaxis, :]
    rg['Tx_1211'][:] = VEL['temp'][np.newaxis, np.newaxis, :]
    rg['AGC_1202'][:] = np.reshape(VEL['AGC'].T, (M, 1, 1, N))

    rg.close()

    print 'Done writing NetCDF file'
