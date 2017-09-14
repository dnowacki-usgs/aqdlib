from __future__ import division, print_function
from netCDF4 import Dataset
import numpy as np
import qaqc
import dateutil
import pytz
import jdcal
import datetime as dt
import aqdlib
import xarray as xr
import pandas as pd

def cdf_to_nc(cdf_filename, metadata, p_1ac=False):

    nc_filename = metadata['filename'] + '.nc'

    VEL = {}

    VEL, INFO = load_cdf_amp_vel(cdf_filename, VEL, metadata, p_1ac=p_1ac)

    # define_aqd_nc_file(nc_filename, VEL, metadata, INFO)

    # write_aqd_nc_file(nc_filename, VEL, metadata)

    # TODO: Need to add all global attributes from CDF to NC file (or similar)
    qaqc.add_min_max(VEL)
    print('Added min/max values')

    qaqc.add_final_metadata(VEL)
    print('Added final metadata')

    # VEL = VEL.expand_dims('lat')
    # need to associate dimensions with lat/lon
    for var in ['U', 'V', 'W', 'AGC', 'Pressure', 'Temperature', 'Heading', 'Pitch', 'Roll']:
        VEL[var] = reshape(VEL[var])

    VEL = rename(VEL, metadata, INFO)

    VEL.to_netcdf(nc_filename)
    print('Done writing NetCDF file', nc_filename)

    return VEL

def reshape(da):
    """
    Add lon and lat dimensions to data and reshape to conform to our standard order
    """

    da = xr.concat([da], dim='lon')
    da = xr.concat([da], dim='lat')
    if len(da.shape) == 4:
        da = da.transpose('time', 'lon', 'lat', 'bindist')
    elif len(da.shape) == 3:
        da = da.transpose('time', 'lon', 'lat')

    return da

def rename(ds, metadata, INFO):
    """
    Rename DataArrays, drop unused DataArrays, and add attrs
    """

    # put in attributes common to all velocity components
    def add_vel_attributes(vel, metadata, INFO):
        vel.attrs.update({'units': 'cm/s',
            'data_cmnt': 'Velocity in shallowest bin is often suspect and should be used with caution'})

        add_attributes(vel, metadata, INFO)
        # TODO: why do we only do trim_method for Water Level SL?
        if 'trim_method' in metadata and metadata['trim_method'].lower() == 'water level sl':
            vel.attrs.update({'note': 'Velocity bins trimmed if out of water or if side lobes intersect sea surface'})

    def add_attributes(var, metadata, INFO):
        var.attrs.update({'serial_number': INFO['AQDSerial_Number'],
            'initial_instrument_height': metadata['initial_instrument_height'],
            'nominal_instrument_depth': metadata['nominal_instrument_depth'],
            'height_depth_units': 'm', 'sensor_type': INFO['INST_TYPE']})

    varnames = {'U': 'u_1205',
            'V': 'v_1206',
            'W': 'w_1204',
            'AGC': 'AGC_1202',
            'Pressure': 'P_1',
            'Temperature': 'Tx_1211',
            'Heading': 'Hdg_1215',
            'Pitch': 'Ptch_1216',
            'Roll': 'Roll_1217'}

    ds.rename(varnames, inplace=True)

    todrop = ['VEL1', 'VEL2', 'VEL3', 'AMP1', 'AMP2', 'AMP3', 'Battery', 'TransMatrix', 'AnalogInput1', 'AnalogInput2']
    ds = ds.drop(todrop)

    ds['u_1205'].attrs.update({'name': 'u', 'long_name': 'Eastward Velocity', 'generic_name': 'u', 'epic_code': 1205})
    ds['v_1206'].attrs.update({'name': 'v', 'long_name': 'Northward Velocity', 'generic_name': 'v', 'epic_code': 1206})
    ds['w_1204'].attrs.update({'name': 'w', 'long_name': 'Vertical Velocity', 'generic_name': 'w', 'epic_code': 1204})
    ds['P_1'].attrs.update({'units': 'dbar', 'name': 'P', 'long_name': 'Pressure', 'generic_name': 'depth', 'epic_code': 1}) # TODO: is this generic name correct?
    ds['Tx_1211'].attrs.update({'units': 'C', 'name': 'Tx', 'long_name': 'Instrument Transducer Temperature', 'generic_name': 'temp', 'epic_code': 1211})
    ds['AGC_1202'].attrs.update({'units': 'counts', 'name': 'AGC', 'long_name': 'Average Echo Intensity', 'generic_name': 'AGC', 'epic_code': 1202})
    # AGCid.sensor_type = INFO['INST_TYPE']

    ds['Hdg_1215'].attrs.update({'units': 'degrees', 'name': 'Hdg', 'long_name': 'Instrument Heading', 'generic_name': 'hdg', 'epic_code': 1215})
    if 'magnetic_variation_at_site' in metadata:
        ds['Hdg_1215'].attrs.update({'note': 'Heading is degrees true. Converted from magnetic with magnetic variation of ' + str(metadata['magnetic_variation_at_site'])})
    elif 'magnetic_variation' in metadata:
        ds['Hdg_1215'].attrs.update({'note': 'Heading is degrees true. Converted from magnetic with magnetic variation of ' + str(metadata['magnetic_variation'])})

    ds['Ptch_1216'].attrs.update({'units': 'degrees', 'name': 'Ptch', 'long_name': 'Instrument Pitch', 'generic_name': 'ptch', 'epic_code': 1216})
    ds['Roll_1217'].attrs.update({'units': 'degrees', 'name': 'Roll', 'long_name': 'Instrument Roll', 'generic_name': 'roll', 'epic_code': 1217})

    for v in ['P_1', 'Tx_1211', 'AGC_1202', 'Hdg_1215', 'Ptch_1216', 'Roll_1217']:
        add_attributes(ds[v], metadata, INFO)

    for v in ['u_1205', 'v_1206', 'w_1204']:
        add_vel_attributes(ds[v], metadata, INFO)


    return ds

def load_cdf_amp_vel(cdf_filename, VEL, metadata, p_1ac=False):

    ds = xr.open_dataset(cdf_filename, autoclose=True)

    # Create INFO dict with global attributes from dataset
    INFO = ds.attrs

    # clip either by ensemble indices or by the deployment and recovery date specified in metadata
    if 'good_ens' in metadata:
        # we have good ensemble indices in the metadata
        print('Using good_ens')
        S = metadata['good_ens'][0]
        E = metadata['good_ens'][1]

        print('first burst in full file:', ds['time'].min().values)
        print('last burst in full file:', ds['time'].max().values)

        ds = ds.isel(time=slice(S,E))

        print('first burst in trimmed file:', ds['time'].min().values)
        print('last burst in trimmed file:', ds['time'].max().values)

    else:
        # we clip by the times in/out of water as specified in the metadata
        print('Using Deployment_date and Recovery_date')

        print('first burst in full file:', ds['time'].min().values)
        print('last burst in full file:', ds['time'].max().values)

        ds = ds.sel(time=slice(metadata['Deployment_date'], metadata['Recovery_date']))

        print('first burst in trimmed file:', ds['time'].min().values)
        print('last burst in trimmed file:', ds['time'].max().values)

    print(ds.info)

    VEL, metadata = qaqc.create_water_depth(ds, metadata)

    T = ds['TransMatrix'].values

    VEL, T = qaqc.set_orientation(ds, T, metadata, INFO)

    u, v, w = qaqc.coord_transform(VEL['VEL1'].values, VEL['VEL2'].values, VEL['VEL3'].values,
        VEL['Heading'].values, VEL['Pitch'].values, VEL['Roll'].values, T, VEL.attrs['AQDCoordinateSystem'])

    VEL['U'] = xr.DataArray(u, dims=('time', 'bindist'))
    VEL['V'] = xr.DataArray(v, dims=('time', 'bindist'))
    VEL['W'] = xr.DataArray(w, dims=('time', 'bindist'))

    VEL = qaqc.magvar_correct(VEL, metadata)

    VEL['AGC'] = (VEL['AMP1'] + VEL['AMP2'] + VEL['AMP3']) / 3

    VEL = qaqc.trim_vel(VEL, metadata, INFO)

    VEL = qaqc.make_bin_depth(VEL, metadata)


    return VEL, INFO

    #     if p_1ac is not False:
    #         VEL['press_ac'] = p_1ac[S:E]
    #
    return VEL, INFO

def define_aqd_nc_file(nc_filename, VEL, metadata, INFO):

    print(VEL['bindist'])

    # Assign COMPOSITE global attribute (formerly assigned at end)
    VEL.attrs.update({'COMPOSITE': 0})

    # VEL['depth'].attrs.update({'epic_code': 3})

    with Dataset(nc_filename, 'w', format='NETCDF4', clobber=True) as rg:

        # time = rg.createDimension('time', 0)
        # depth = rg.createDimension('depth', M)
        # lat = rg.createDimension('lat', 1)
        # lon = rg.createDimension('lon', 1)
        # Tmatrix = rg.createDimension('Tmatrix', 3)
        #
        # timeid = rg.createVariable('time', 'i', ('time',), zlib=True) # 'i' == NC_INT
        # timeid.units = 'True Julian Day'
        # timeid.type = 'EVEN' # TODO: this is "UNEVEN" in the CDF version
        # timeid.epic_code = 624
        #
        # time2id = rg.createVariable('time2', 'i', ('time',), zlib=True) # 'i' == NC_INT
        # time2id.units = 'msec since 0:00 GMT'
        # time2id.type ='EVEN'
        # time2id.epic_code = 624

        # TODO: depth seems different from the cdf variable. this is a mean water depth
        # depthid = rg.createVariable('depth', 'd', ('depth',), zlib=True)
        # depthid.units = 'm'
        # depthid.epic_code = 3
        # depthid.long_name = 'mean water depth'
        # depthid.initial_instrument_height = metadata['initial_instrument_height']
        # depthid.nominal_instrument_depth = metadata['nominal_instrument_depth']

        # TODO: figure out how to cast DOUBLE_FILL to float
        # TODO: why are these created as double precision, when Ellyn says they should be singles?
        bindistid = rg.createVariable('bindist', 'd', ('depth',), zlib=True, fill_value=aqdlib.DOUBLE_FILL)
        # TODO: Loop through the attribute names and values as done in Matlab
        bindistid.blanking_distance = INFO['AQDBlankingDistance']
        bindistid.initial_instrument_height = metadata['initial_instrument_height']
        bindistid.note = 'distance is along profile from instrument head to center of bin'




        add_vel_attributes(u_1205, metadata, INFO)

        add_vel_attributes(v_1206, metadata, INFO)


        add_vel_attributes(w_1204, metadata, INFO)


        add_attributes(Pressid, metadata, INFO)
        # TODO: why no sensor_type in these vars?



        # TODO: add analog input variables (OBS, NTU, etc)

        bdepid = rg.createVariable('bin_depth', 'd', ('depth', 'lat', 'lon', 'time',), zlib=True, fill_value=aqdlib.DOUBLE_FILL)
        bdepid.setncattr('name', 'bin depth')
        bdepid.units = 'm'
        bdepid.initial_instrument_height = metadata['initial_instrument_height']
        add_attributes(bdepid, metadata, INFO)
        if 'press_ac' in VEL:
            bdepid.note = 'Actual depth time series of velocity bins. Calculated as corrected pressure(P_1ac) - bindist.'
        else:
            bdepid.note = 'Actual depth time series of velocity bins. Calculated as pressure(P_1) - bindist.'

        if 'press_ac' in VEL:
            ACPressid = rg.createVariable('P_1ac', 'd', ('lat', 'lon', 'time',), zlib=True, fill_value=aqdlib.DOUBLE_FILL)
            ACPressid.units = 'dbar'
            ACPressid.setncattr('name', 'Pac') # TODO: is Pac correct?
            ACPressid.long_name = 'CORRECTED PRESSURE (DB)'
            add_attributes(ACPressid, metadata, INFO)
            ACPressid.note = 'Corrected for variations in atmospheric pressure using nearby MET station'

def write_aqd_nc_file(nc_filename, VEL, metadata):

    N, M = np.shape(VEL['U'])
    print('N:', N, 'M:', M, 'in write_aqd_nc_file')

    with Dataset(nc_filename, 'r+') as rg:

        rg['lat'][:] = metadata['latitude']
        rg['lon'][:] = metadata['longitude']
        rg['time'][:] = VEL['time']
        rg['time2'][:] = VEL['time2']
        rg['depth'][:] = VEL['depths']

        rg['bindist'][:] = VEL['bindist']
        rg['bin_depth'][:] = np.reshape(VEL['bin_depth'].T, (M, 1, 1, N))
        rg['u_1205'][:] = np.reshape(VEL['U'].T, (M, 1, 1, N))
        rg['v_1206'][:] = np.reshape(VEL['V'].T, (M, 1, 1, N))
        rg['w_1204'][:] = np.reshape(VEL['W'].T, (M, 1, 1, N))
        rg['AGC_1202'][:] = np.reshape(VEL['AGC'].T, (M, 1, 1, N))

        rg['Tx_1211'][:] = VEL['temp'][np.newaxis, np.newaxis, :]
        rg['P_1'][:] = VEL['pressure'][np.newaxis, np.newaxis, :]

        if 'press_ac' in VEL:
            rg['P_1ac'][:] = VEL['press_ac'][np.newaxis, np.newaxis, :]

        rg['Ptch_1216'][:] = VEL['pitch'][np.newaxis, np.newaxis, :]
        rg['Roll_1217'][:] = VEL['roll'][np.newaxis, np.newaxis, :]
        rg['Hdg_1215'][:] = VEL['heading'][np.newaxis, np.newaxis, :]
