from __future__ import division, print_function
from netCDF4 import Dataset
import numpy as np
import qaqc
from aqdlib import DOUBLE_FILL
import dateutil
import pytz
import jdcal
import datetime as dt

def cdf_to_nc(cdf_filename, metadata, p_1ac=False):

    nc_filename = metadata['filename'] + '.nc'

    VEL = {}

    VEL = load_cdf_amp_vel(cdf_filename, VEL, metadata, p_1ac=p_1ac)

    define_aqd_nc_file(nc_filename, VEL, metadata)

    write_aqd_nc_file(nc_filename, VEL, metadata)

    return VEL

def load_cdf_amp_vel(cdf_filename, VEL, metadata, p_1ac=False):
    try:
        rg = Dataset(cdf_filename, 'r')

        # clip either by ensemble indices or by the deployment and recovery date specified in metadata
        if 'good_ens' in metadata:
            # we have good ensemble indices in the metadata
            print('Using good_ens')
            S = metadata['good_ens'][0]
            E = metadata['good_ens'][1]
        else:
            # we clip by the times in/out of water as specified in the metadata
            print('Using Deployment_date and Recovery_date')
            time = rg['time'][:]
            time2 = rg['time2'][:]

            times = []
            for t, t2 in zip(time, time2):
                year, mon, day, frac = jdcal.jd2gcal(t - 0.5, t2/86400000)
                hour, minute, second = qaqc.day2hms(frac)
                times.append(dt.datetime(year, mon, day, hour, minute, second, tzinfo=pytz.utc))
            times = np.array(times)

            S = np.argwhere(times >= pytz.utc.localize(dateutil.parser.parse(rg.Deployment_date)))[0]
            E = np.argwhere(times <= pytz.utc.localize(dateutil.parser.parse(rg.Recovery_date)))[-1]

        print('Indices of starting and ending bursts: S:', S, 'E:', E)

        # load data from CDF file, specifying start/end bursts
        # also transpose data so dims are TIME x DEPTH
        vel1 = rg['VEL1'][:, S:E].T
        vel2 = rg['VEL2'][:, S:E].T
        vel3 = rg['VEL3'][:, S:E].T

        amp1 = rg['AMP1'][:, S:E].T
        amp2 = rg['AMP2'][:, S:E].T
        amp3 = rg['AMP3'][:, S:E].T

        heading = rg['Heading'][S:E]
        pitch = rg['Pitch'][S:E]
        roll = rg['Roll'][S:E]

        VEL['bindist'] = rg['bindist'][:]

        T = rg['TransMatrix'][:]

        VEL['pressure'] = rg['Pressure'][S:E]
        VEL['temp'] = rg['Temperature'][S:E]
        VEL['time'] = rg['time'][S:E]
        VEL['time2'] = rg['time2'][S:E]

        if p_1ac is not False:
            VEL['press_ac'] = p_1ac[S:E]


        VEL, metadata = qaqc.create_water_depth(VEL, metadata)

        # initialize arrays
        VEL['U'] = np.zeros(np.shape(vel1))
        VEL['V'] = np.zeros(np.shape(vel2))
        VEL['W'] = np.zeros(np.shape(vel3))

        VEL, T = qaqc.set_orientation(VEL, T, metadata)

        VEL['U'], VEL['V'], VEL['W'] = qaqc.coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, rg.AQDCoordinateSystem)

        VEL['heading'] = heading
        VEL['pitch'] = pitch
        VEL['roll'] = roll

        VEL = qaqc.magvar_correct(VEL, metadata)

        VEL['AGC'] = (amp1 + amp2 + amp3) / 3

        VEL = qaqc.trim_vel(VEL, metadata)

        VEL = qaqc.make_bin_depth(VEL, metadata)

        return VEL

    finally:
        rg.close()



def define_aqd_nc_file(nc_filename, VEL, metadata):
    try:
        N, M = np.shape(VEL['U'])
        print('N:', N, 'M:', M, 'in define_aqd_nc_file')

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
        # TODO: why is the setup of these variables different from the setup in the CDF routines?

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

        # put in attributes common to all velocity components
        def add_vel_attributes(vel, metadata):
            vel.units = 'cm/s'
            vel.sensor_type = metadata['INST_TYPE']
            vel.initial_instrument_height = metadata['initial_instrument_height']
            vel.nominal_instrument_depth = metadata['nominal_instrument_depth']
            vel.serial_number = metadata['serial_number']
            vel.maximum = 0
            vel.minimum = 0
            vel.height_depth_units = 'm'
            # TODO: why do we only do trim_method for Water Level SL?
            if 'trim_method' in metadata and metadata['trim_method'].lower() == 'water level sl':
                vel.note ='Velocity bins trimmed if out of water or if side lobes intersect sea surface'

        u_1205 = rg.createVariable('u_1205', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
        u_1205.setncattr('name', 'u')
        u_1205.long_name = 'Eastward Velocity'
        u_1205.generic_name = 'u'
        u_1205.epic_code = 1205
        add_vel_attributes(u_1205, metadata)

        v_1206 = rg.createVariable('v_1206', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
        v_1206.setncattr('name', 'v')
        v_1206.long_name = 'Northward Velocity'
        v_1206.generic_name = 'v'
        v_1206.epic_code = 1206
        add_vel_attributes(v_1206, metadata)

        w_1204 = rg.createVariable('w_1204', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
        w_1204.setncattr('name', 'w')
        w_1204.long_name = 'Vertical Velocity'
        w_1204.generic_name = 'w'
        w_1204.epic_code = 1204
        add_vel_attributes(w_1204, metadata)

        Pressid = rg.createVariable('P_1', 'f', ('lat', 'lon', 'time',), zlib=True)
        Pressid.units = 'dbar'
        Pressid.epic_code = 1
        Pressid.setncattr('name', 'P')
        Pressid.long_name = 'PRESSURE (DB)'
        Pressid.generic_name = 'depth'
        Pressid.minimum = DOUBLE_FILL
        Pressid.maximum = DOUBLE_FILL
        # netcdf.putAtt(ncid,Pressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); #FIXME
        Pressid.initial_instrument_height = metadata['initial_instrument_height']
        Pressid.nominal_instrument_depth = metadata['nominal_instrument_depth']

        Tempid = rg.createVariable('Tx_1211', 'f', ('lat', 'lon', 'time',), zlib=True)
        Tempid.units = 'C'
        Tempid.epic_code = 1211
        Tempid.setncattr('name', 'Tx')
        Tempid.long_name = 'Instrument Transducer Temperature'
        Tempid.generic_name = 'temp'
        Tempid.minimum = DOUBLE_FILL
        Tempid.maximum = DOUBLE_FILL
        # netcdf.putAtt(ncid,Pressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); #FIXME
        Tempid.initial_instrument_height = metadata['initial_instrument_height']
        Tempid.nominal_instrument_depth = metadata['nominal_instrument_depth']

        AGCid = rg.createVariable('AGC_1202', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
        AGCid.units = 'counts'
        AGCid.epic_code = 1202
        AGCid.setncattr('name', 'AGC')
        AGCid.long_name = 'Average Echo Intensity (AGC)'
        AGCid.generic_name = 'AGC'
        # AGCid.sensor_type = INST_TYPE # FIXME
        AGCid.minimum = 0
        AGCid.maximum = 0 # TODO: why are min/max different (0 vs DOUBLE_FILL for others?)
        # netcdf.putAtt(ncid,AGCid,'sensor_type',metadata.cdfmeta.INST_TYPE); #FIXME
        # netcdf.putAtt(ncid,Pressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); #FIXME
        AGCid.initial_instrument_height = metadata['initial_instrument_height']
        AGCid.height_depth_units = 'm'
        AGCid.nominal_instrument_depth = metadata['nominal_instrument_depth']

        # TODO: add heading/pitch/roll variables

        # TODO: add analog input variables (OBS, NTU, etc)

        bdepid = rg.createVariable('bin_depth', 'f', ('depth', 'lat', 'lon', 'time',), zlib=True)
        bdepid.setncattr('name', 'bin depth')
        bdepid.units = 'm'
        bdepid.initial_instrument_height = metadata['initial_instrument_height']
        bdepid.minimum = 0 # why 0 instead of double_fill??
        bdepid.maximum = 0
        if 'press_ac' in VEL:
            bdepid.note = 'Actual depth time series of velocity bins. Calculated as corrected pressure(P_1ac) - bindist.'
        else:
            bdepid.note = 'Actual depth time series of velocity bins. Calculated as pressure(P_1) - bindist.'

        if 'press_ac' in VEL:
            ACPressid = rg.createVariable('P_1ac', 'f', ('lat', 'lon', 'time',), zlib=True)
            ACPressid.units = 'dbar'
            ACPressid.setncattr('name', 'Pac') # TODO: is Pac correct?
            ACPressid.long_name = 'CORRECTED PRESSURE (DB)'
            ACPressid.minimum = DOUBLE_FILL # why double_fill??
            ACPressid.maximum = DOUBLE_FILL
            # netcdf.putAtt(ncid,CPressid,'serial_number',metadata.cdfmeta.AQDSerial_Number); # FIXME
            ACPressid.initial_instrument_height = metadata['initial_instrument_height']
            ACPressid.nominal_instrument_depth = metadata['nominal_instrument_depth']
            ACPressid.note = 'Corrected for variations in atmospheric pressure using nearby MET station'

    finally:
        rg.close()

def write_aqd_nc_file(nc_filename, VEL, metadata):

    try:
        N, M = np.shape(VEL['U'])
        print('N:', N, 'M:', M, 'in write_aqd_nc_file')
        rg = Dataset('/Volumes/Backstaff/field/gb_proc/1076a/1076a1aqd/' + nc_filename, 'r+')

        rg['lat'][:] = metadata['latitude']
        rg['lon'][:] = metadata['longitude']
        rg['time'][:] = VEL['time']
        rg['time2'][:] = VEL['time2']

        rg['bindist'][:] = VEL['bindist']
        rg['bin_depth'][:] = np.reshape(VEL['bin_depth'].T, (M, 1, 1, N))
        rg['u_1205'][:] = np.reshape(VEL['U'].T, (M, 1, 1, N))
        rg['v_1206'][:] = np.reshape(VEL['V'].T, (M, 1, 1, N))
        rg['w_1204'][:] = np.reshape(VEL['W'].T, (M, 1, 1, N))
        rg['AGC_1202'][:] = np.reshape(VEL['AGC'].T, (M, 1, 1, N))

        rg['P_1'][:] = VEL['pressure'][np.newaxis, np.newaxis, :]

        if 'press_ac' in VEL:
            rg['P_1ac'][:] = VEL['press_ac'][np.newaxis, np.newaxis, :]

        rg['Tx_1211'][:] = VEL['temp'][np.newaxis, np.newaxis, :]

    finally:
        print('Done writing NetCDF file')
        rg.close()
