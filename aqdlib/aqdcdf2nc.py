from __future__ import division, print_function
import xarray as xr
import pandas as pd
import numpy as np
import qaqc
import aqdlib

def cdf_to_nc(cdf_filename, metadata, atmpres=False):
    """
    Load a "raw" .cdf file and generate a processed .nc file
    """

    VEL = load_cdf(cdf_filename, metadata, atmpres=atmpres)

    # Create INFO dict with global attributes from dataset
    INFO = VEL.attrs

    # Clip data to in/out water times or via good_ens
    VEL = clip_ds(VEL, metadata)

    # Create water_depth variables
    VEL, metadata = qaqc.create_water_depth(VEL, metadata)

    VEL, T = qaqc.set_orientation(VEL, VEL['TransMatrix'].values, metadata, INFO)

    u, v, w = qaqc.coord_transform(VEL['VEL1'].values, VEL['VEL2'].values, VEL['VEL3'].values,
        VEL['Heading'].values, VEL['Pitch'].values, VEL['Roll'].values, T, VEL.attrs['AQDCoordinateSystem'])

    VEL['U'] = xr.DataArray(u, dims=('time', 'bindist'))
    VEL['V'] = xr.DataArray(v, dims=('time', 'bindist'))
    VEL['W'] = xr.DataArray(w, dims=('time', 'bindist'))

    VEL = qaqc.magvar_correct(VEL, metadata)

    VEL['AGC'] = (VEL['AMP1'] + VEL['AMP2'] + VEL['AMP3']) / 3

    VEL = qaqc.trim_vel(VEL, metadata, INFO)

    VEL = qaqc.make_bin_depth(VEL, metadata)

    # TODO: Need to add all global attributes from CDF to NC file (or similar)
    VEL = qaqc.add_min_max(VEL)

    VEL = qaqc.add_final_metadata(VEL)

    # Reshape and associate dimensions with lat/lon
    for var in ['U', 'V', 'W', 'AGC', 'Pressure', 'Temperature', 'Heading', 'Pitch', 'Roll']:
        VEL = da_reshape(VEL, var)

    # swap_dims from bindist to depth
    VEL = ds_swap_dims(VEL)

    VEL = ds_rename(VEL, metadata, INFO)

    VEL = ds_drop(VEL, metadata, INFO)

    VEL = ds_add_attrs(VEL, metadata, INFO)

    nc_filename = metadata['filename'] + '.nc'

    VEL.to_netcdf(nc_filename)
    print('Done writing NetCDF file', nc_filename)

    print(VEL['u_1205'].attrs)

    return VEL

def load_cdf(cdf_filename, metadata, atmpres=False):

    ds = xr.open_dataset(cdf_filename, autoclose=True)

    if atmpres is not False:
        p = xr.open_dataset(atmpres)
        # TODO: check to make sure this data looks OK
        ds['Pressure_ac'] = xr.DataArray(ds['Pressure'] - (p['atmpres'] - p['atmpres'].offset))

    return ds

def clip_ds(ds, metadata):

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

    return ds

def define_aqd_nc_file(nc_filename, VEL, metadata, INFO):

    # Assign COMPOSITE global attribute (formerly assigned at end)
    VEL.attrs.update({'COMPOSITE': 0})

    # VEL['depth'].attrs.update({'epic_code': 3})

    with Dataset(nc_filename, 'w', format='NETCDF4', clobber=True) as rg:

        # TODO: depth seems different from the cdf variable. this is a mean water depth
        # depthid = rg.createVariable('depth', 'd', ('depth',), zlib=True)


        # TODO: figure out how to cast DOUBLE_FILL to float
        # TODO: why are these created as double precision, when Ellyn says they should be singles?
        bindistid = rg.createVariable('bindist', 'd', ('depth',), zlib=True, fill_value=aqdlib.DOUBLE_FILL)
        # TODO: Loop through the attribute names and values as done in Matlab
        bindistid.blanking_distance = INFO['AQDBlankingDistance']
        bindistid.initial_instrument_height = metadata['initial_instrument_height']
        bindistid.note = 'distance is along profile from instrument head to center of bin'

        # TODO: add analog input variables (OBS, NTU, etc)

def ds_swap_dims(ds):

    # ds['depth'] = xr.DataArray(ds['depth'].values, dims=('bindist'))
    ds = ds.swap_dims({'bindist': 'depth'})
    # ds['bindist'] = ds['bindist'].swap_dims({'bindist': 'depth'})

    return ds

def da_reshape(ds, var):
    """
    Add lon and lat dimensions to DataArrays and reshape to conform to our
    standard order
    """

    # Add the dimensions using concat
    ds[var] = xr.concat([ds[var]], dim=ds['lon'])
    ds[var] = xr.concat([ds[var]], dim=ds['lat'])

    # Reshape using transpose depending on shape
    if len(ds[var].shape) == 4:
        ds[var] = ds[var].transpose('time', 'lon', 'lat', 'bindist')
    elif len(ds[var].shape) == 3:
        ds[var] = ds[var].transpose('time', 'lon', 'lat')

    return ds

def ds_rename(ds, metadata, INFO):
    """
    Rename DataArrays within Dataset for EPIC compliance
    """

    varnames = {'U': 'u_1205',
        'V': 'v_1206',
        'W': 'w_1204',
        'AGC': 'AGC_1202',
        'Pressure': 'P_1',
        'Temperature': 'Tx_1211',
        'Heading': 'Hdg_1215',
        'Pitch': 'Ptch_1216',
        'Roll': 'Roll_1217'}

    if 'Pressure_ac' in ds:
        varnames['Pressure_ac'] = 'P_1ac'

    ds.rename(varnames, inplace=True)

    return ds

def ds_drop(ds, metadata, INFO):
    """
    Drop old DataArrays from Dataset that won't make it into the final .nc file
    """

    todrop = ['VEL1',
        'VEL2',
        'VEL3',
        'AMP1',
        'AMP2',
        'AMP3',
        'Battery',
        'TransMatrix',
        'AnalogInput1',
        'AnalogInput2']

    ds = ds.drop(todrop)

    return ds

def ds_add_attrs(ds, metadata, INFO):
    """
    add attrs
    """

    def add_vel_attributes(vel, metadata, INFO):
        vel.attrs.update({'units': 'cm/s',
            'data_cmnt': 'Velocity in shallowest bin is often suspect and should be used with caution'})

        # TODO: why do we only do trim_method for Water Level SL?
        if 'trim_method' in metadata and metadata['trim_method'].lower() == 'water level sl':
            vel.attrs.update({'note': 'Velocity bins trimmed if out of water or if side lobes intersect sea surface'})

    def add_attributes(var, metadata, INFO):
        var.attrs.update({'serial_number': INFO['AQDSerial_Number'],
            'initial_instrument_height': metadata['initial_instrument_height'],
            'nominal_instrument_depth': metadata['nominal_instrument_depth'],
            'height_depth_units': 'm', 'sensor_type': INFO['INST_TYPE'],
            '_FillValue': 1e35})

    # Update attributes for EPIC and STG compliance
    ds.lat.encoding['_FillValue'] = 1e35

    ds.lon.encoding['_FillValue'] = 1e35

    ds['epic_time'].attrs.update({'units': 'True Julian Day',
        'type': 'EVEN',
        'epic_code': 624})

    ds['epic_time2'].attrs.update({'units': 'msec since 0:00 GMT',
        'type': 'EVEN',
        'epic_code': 624})

    ds['u_1205'].attrs.update({'name': 'u',
        'long_name': 'Eastward Velocity',
        'generic_name': 'u',
        'epic_code': 1205})

    ds['v_1206'].attrs.update({'name': 'v',
        'long_name': 'Northward Velocity',
        'generic_name': 'v',
        'epic_code': 1206})

    ds['w_1204'].attrs.update({'name': 'w',
        'long_name': 'Vertical Velocity',
        'generic_name': 'w',
        'epic_code': 1204})

    ds['P_1'].attrs.update({'units': 'dbar',
        'name': 'P',
        'long_name': 'Pressure',
        'generic_name': 'depth',
        'epic_code': 1}) # TODO: is this generic name correct?

    if 'P_1ac' in ds:
        ds['P_1ac'].attrs.update({'units': 'dbar',
            'name': 'Pac',
            'long_name': 'Corrected pressure',
            'note': 'Corrected for variations in atmospheric pressure using nearby met station'})
        add_attributes(ds['P_1ac'], metadata, INFO)

    ds['depth'].attrs.update({'units': 'm',
        'long_name': 'mean water depth',
        'initial_instrument_height': metadata['initial_instrument_height'],
        'nominal_instrument_depth': metadata['nominal_instrument_depth'],
        '_FillValue': 1e35,
        'epic_code': 3})


    ds['bin_depth'].attrs.update({'units': 'm',
        'name': 'bin depth'})

    if 'P_1ac' in ds:
        ds['bin_depth'].attrs.update({'note': 'Actual depth time series of velocity bins. Calculated as corrected pressure(P_1ac) - bindist.'})
    else:
        ds['bin_depth'].attrs.update({'note': 'Actual depth time series of velocity bins. Calculated as pressure(P_1) - bindist.'})

    ds['Tx_1211'].attrs.update({'units': 'C',
        'name': 'Tx',
        'long_name': 'Instrument Transducer Temperature',
        'generic_name': 'temp',
        'epic_code': 1211})

    ds['AGC_1202'].attrs.update({'units': 'counts',
        'name': 'AGC',
        'long_name': 'Average Echo Intensity',
        'generic_name': 'AGC',
        'epic_code': 1202})

    ds['Hdg_1215'].attrs.update({'units': 'degrees',
        'name': 'Hdg',
        'long_name': 'Instrument Heading',
        'generic_name': 'hdg',
        'epic_code': 1215})

    if 'magnetic_variation_at_site' in metadata:
        ds['Hdg_1215'].attrs.update({'note': 'Heading is degrees true. Converted from magnetic with magnetic variation of ' + str(metadata['magnetic_variation_at_site'])})
    elif 'magnetic_variation' in metadata:
        ds['Hdg_1215'].attrs.update({'note': 'Heading is degrees true. Converted from magnetic with magnetic variation of ' + str(metadata['magnetic_variation'])})

    ds['Ptch_1216'].attrs.update({'units': 'degrees',
        'name': 'Ptch',
        'long_name': 'Instrument Pitch',
        'generic_name': 'ptch',
        'epic_code': 1216})
    ds['Roll_1217'].attrs.update({'units': 'degrees',
        'name': 'Roll',
        'long_name': 'Instrument Roll',
        'generic_name': 'roll',
        'epic_code': 1217})

    for v in ['P_1', 'Tx_1211', 'AGC_1202', 'Hdg_1215', 'Ptch_1216', 'Roll_1217', 'u_1205', 'v_1206', 'w_1204', 'bin_depth']:
        add_attributes(ds[v], metadata, INFO)

    for v in ['u_1205', 'v_1206', 'w_1204']:
        add_vel_attributes(ds[v], metadata, INFO)

    return ds
