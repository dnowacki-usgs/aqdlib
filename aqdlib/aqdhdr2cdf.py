from __future__ import division
import os
import numpy as np
import datetime as dt
import jdcal
import pytz
from netCDF4 import Dataset

def load(basefile, metadata):
    """
    Main load file
    """
    doublefill = 1e35
    shortfill = -32768

    # get instrument metadata from the HDR file
    instmeta = read_aqd_hdr(basefile)

    print "Loading ASCII files"

    # basePath = os.path.dirname(basefile)
    # baseName = os.path.basename(basefile)
    RAW = {}

    # Load sensor data
    RAW = load_sen(RAW, basefile, metadata)

    # Load amplitude and velocity data
    RAW = load_amp_vel(RAW, basefile)

    # Compute time stamps
    RAW = compute_time(RAW)

    # Deal with metadata peculiarities
    metadata = check_metadata(metadata, instmeta)

    RAW = check_orientation(RAW, metadata)

    # TODO: clock drift code
    # TODO: Move time to center of ensemble??
    # TODO: logmeta code

    RAW['instmeta'] = instmeta

    # Put in fill values
    print "about to insert fill values"
    RAW = insert_fill_values(RAW, doublefill)

    # configure file
    cdf_filename = '/Volumes/Backstaff/field/gb_proc/1076a/1076a1aqd/' + metadata['filename'] + '-raw.cdf' # TODO: fix the path
    print 'Opening %s' % cdf_filename

    define_aqd_cdf_file(cdf_filename, RAW, metadata)
    print 'Variables created'

    write_aqd_cdf_data(cdf_filename, RAW, metadata)
    print 'Variables written'

    add_min_max(cdf_filename)

    print 'Finished writing data to %s' % cdf_filename

    return RAW

def insert_fill_values(RAW, doublefill):
    print "Inserting fill values"
    for k in RAW:
        if k not in ['instmeta', 'time', 'time2', 'datetime'] and np.max(np.shape(RAW[k])) == np.max(np.shape(RAW['jd'])):
            nanind = np.where(np.isnan(RAW[k]))
            RAW[k][nanind] = doublefill

    return RAW

def check_orientation(RAW, metadata):
    print metadata['orientation']
    print 'Center_first_bin = %f\n' % metadata['center_first_bin']
    print 'bin_size = %f\n' % metadata['bin_size']
    print 'bin_count = %f\n' % metadata['bin_count']
    # TODO: these values are already in the HDR file...
    RAW['bindist'] = np.arange(metadata['center_first_bin'], (metadata['center_first_bin']+((metadata['bin_count']-1)*metadata['bin_size'])), metadata['bin_size'])

    if metadata['orientation'] == 'UP':
        print 'User instructed that instrument was pointing UP'
        # depth, or distance below surface, is a positive number below the
        # surface, negative above the surface, for CMG purposes and consistency with ADCP
        RAW['Depths'] = (metadata['WATER_DEPTH'] - metadata['transducer_offset_from_bottom']) - RAW['bindist']
        Depth_NOTE = 'user reports uplooking bin depths = water_depth - transducer offset from bottom - bindist' # TODO: this is never used
    elif metadata['orientation'] == 'DOWN':
        print 'User instructed that instrument was pointing DOWN'
        RAW['Depths'] = (metadata['WATER_DEPTH'] - metadata['transducer_offset_from_bottom']) + RAW['bindist']
        Depth_NOTE = 'user reports downlooking bin depths = water_depth - transducer_offset_from_bottom + bindist' # TODO: this is never used

    return RAW

def check_metadata(metadata, instmeta):

    # % now verify the global metadata for standard EPIC and cmg stuff
    # % everything in metadata and instmeta get written as global attributes
    # % these also get copied to the .nc file
    if 'initial_instrument_height' not in metadata or np.isnan(metadata['initial_instrument_height']):
        metadata['initial_instrument_height'] = 0

    for k in instmeta:
        print k
    metadata['serial_number'] = instmeta['AQDSerial_Number']

    # update metadata from Aquadopp header to CMG standard so that various
    # profilers have the same attribute wording.  Redundant, but necessary
    metadata['bin_count'] = instmeta['AQDNumberOfCells']
    metadata['bin_size'] = instmeta['AQDCellSize'] / 100 # from cm to m
    metadata['blanking_distance'] = instmeta['AQDBlankingDistance'] # already in m
    # Nortek lists the distance to the center of the first bin as the blanking
    # distance plus one cell size
    metadata['center_first_bin'] = metadata['blanking_distance'] + metadata['bin_size'] # in m


    metadata['salinity_set_by_user'] = instmeta['AQDSalinity']
    metadata['salinity_set_by_user_units'] = 'ppt'

    metadata['frequency'] = instmeta['AQDFrequency']
    metadata['beam_width'] = instmeta['AQDBeamWidth']
    metadata['beam_pattern'] = instmeta['AQDBeamPattern']
    metadata['beam_angle'] = instmeta['AQDBeamAngle']
    # instmeta['AQDHeadRotation'] = metadata.pop('head_rotation') # also deletes this key. Not sure why the Matlab file does this, maybe need TODO look into this

    # TODO: figure these out
    # metadata['insterr'] = instmeta['error']
    # metadata['inststat'] = instmeta['status']
    # metadata['instorient'] = instmeta['orient']

    metadata['INST_TYPE'] = 'Nortek Aquadopp Profiler';

    return metadata

def add_min_max(cdf_filename):
    """
    Add minimum and maximum values to variables in CDF file
    """

    rg = Dataset(cdf_filename, 'r+')
    exclude = rg.dimensions.keys()
    exclude.extend(('time2', 'TransMatrix'))
    for var in rg.variables:
        if var not in exclude:
            rg[var].minimum = rg[var][:].min()
            rg[var].maximum = rg[var][:].max()
    print "Assigned min and max values"
    rg.close()

def write_aqd_cdf_data(cdf_filename, RAW, metadata):
    """
    Write data to NetCDF file that has already been set up using
    define_aqd_cdf_file()
    """
    rg = Dataset(cdf_filename, 'r+')

    rg['lat'][:] = metadata['latitude']
    rg['lon'][:] = metadata['longitude']
    rg['time'][:] = RAW['time']
    rg['time2'][:] = RAW['time2']
    rg['depth'][:] = RAW['Depths']
    rg['bindist'][:] = RAW['bindist']

    rg['VEL1'][:] = RAW['V1'].T
    rg['VEL2'][:] = RAW['V2'].T
    rg['VEL3'][:] = RAW['V3'].T

    rg['AMP1'][:] = RAW['AMP1'].T
    rg['AMP2'][:] = RAW['AMP2'].T
    rg['AMP3'][:] = RAW['AMP3'].T

    rg['Temperature'][:] = RAW['temperature']
    rg['Pressure'][:] = RAW['pressure']
    rg['Battery'][:] = RAW['battery']

    rg['Pitch'][:] = RAW['pitch']
    rg['Roll'][:] = RAW['roll']
    rg['Heading'][:] = RAW['heading']

    rg['TransMatrix'][:] = RAW['instmeta']['AQDTransMatrix']

    if 'AnaInp1' in RAW:
        rg['AnalogInput1'][:] = RAW['AnaInp1']

    if 'AnaInp2' in RAW:
        rg['AnalogInput2'][:] = RAW['AnaInp2']

    rg.close()

def define_aqd_cdf_file(cdf_filename, RAW, metadata):
    """
    Define dimensions and variables in NetCDF file
    """
    double_fill_val = 1e35;

    # try:
    rg = Dataset(cdf_filename, 'w', format='NETCDF4', clobber=True)

    # write out EPIC metadata
    write_metadata(rg, metadata)
    write_metadata(rg, RAW['instmeta'])

    N, M = np.shape(RAW['V1'])

    # Time is the record dimension
    time = rg.createDimension('time', 0)
    depth = rg.createDimension('depth', M)
    lat = rg.createDimension('lat', 1)
    lon = rg.createDimension('lon', 1)
    Tmatrix = rg.createDimension('Tmatrix', 3)

    timeid = rg.createVariable('time', 'i', ('time',)) # 'i' == NC_INT
    timeid.FORTRAN_format = 'F10.2'
    timeid.units = 'True Julian Day'
    timeid.type = 'UNEVEN'
    timeid.epic_code = 624

    time2id = rg.createVariable('time2', 'i', ('time',)) # 'i' == NC_INT
    time2id.FORTRAN_format = 'F10.2'
    time2id.units = 'msec since 0:00 GMT'
    time2id.type ='UNEVEN'
    time2id.epic_code = 624

    latid = rg.createVariable('lat', 'f', ('lat',))
    latid.FORTRAN_format = 'F10.4'
    latid.units = 'degree_north'
    latid.type = 'EVEN'
    latid.epic_code = 500
    latid.minimum = double_fill_val
    latid.maximum = double_fill_val

    lonid = rg.createVariable('lon', 'f', ('lon',))
    lonid.FORTRAN_format = 'F10.4'
    lonid.units = 'degree_east'
    lonid.type = 'EVEN'
    lonid.epic_code = 502
    lonid.minimum = double_fill_val
    lonid.maximum = double_fill_val

    depthid = rg.createVariable('depth', 'f', ('depth',))
    depthid.units = 'm'
    depthid.long_name = 'mean water depth'
    depthid.bin_size = metadata['bin_size']
    depthid.center_first_bin = metadata['center_first_bin']
    depthid.bin_count = metadata['bin_count']
    depthid.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    bindistid = rg.createVariable('bindist', 'f', ('depth',))
    bindistid.units = 'm'
    bindistid.long_name = 'distance from transducer head'
    bindistid.bin_size = metadata['bin_size']
    bindistid.center_first_bin = metadata['center_first_bin']
    bindistid.bin_count = metadata['bin_count']
    bindistid.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    Tempid = rg.createVariable('Temperature', 'f', ('time',))
    Tempid.units = 'C'
    Tempid.long_name = 'TEMPERATURE (C)'
    Tempid.generic_name = 'temp'

    Pressid = rg.createVariable('Pressure', 'f', ('time',))
    Pressid.units = 'dbar'
    Pressid.long_name = 'Pressure (dbar)'
    Pressid.generic_name = 'press'

    VEL1id = rg.createVariable('VEL1', 'f', ('depth', 'time',))
    VEL1id.units = 'cm/s'
    VEL1id.Type = 'scalar'
    VEL1id.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    VEL2id = rg.createVariable('VEL2', 'f', ('depth', 'time',))
    VEL2id.units = 'cm/s'
    VEL2id.Type = 'scalar'
    VEL2id.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    VEL3id = rg.createVariable('VEL3', 'f', ('depth', 'time',))
    VEL3id.units = 'cm/s'
    VEL3id.Type = 'scalar'
    VEL3id.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    if RAW['instmeta']['AQDCoordinateSystem'] == 'ENU':
        VEL1id.long_name = 'Eastward current velocity'
        VEL2id.long_name = 'Northward current velocity'
        VEL3id.long_name = 'Vertical current velocity'
    elif RAW['instmeta']['AQDCoordinateSystem'] == 'XYZ':
        VEL1id.long_name = 'Current velocity in X Direction'
        VEL2id.long_name = 'Current velocity in Y Direction'
        VEL3id.long_name = 'Current velocity in Z Direction'
    elif RAW['instmeta']['AQDCoordinateSystem'] == 'BEAM':
        VEL1id.long_name = 'Beam 1 current velocity'
        VEL2id.long_name = 'Beam 2 current velocity'
        VEL3id.long_name = 'Beam 3 current velocity'

    AMP1id = rg.createVariable('AMP1', 'f', ('depth', 'time',))
    AMP1id.long_name = 'Beam 1 Echo Amplitude'
    AMP1id.units = 'counts'
    AMP1id.Type = 'scalar'
    AMP1id.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    AMP2id = rg.createVariable('AMP2', 'f', ('depth', 'time',))
    AMP2id.long_name = 'Beam 2 Echo Amplitude'
    AMP2id.units = 'counts'
    AMP2id.Type = 'scalar'
    AMP2id.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    AMP3id = rg.createVariable('AMP3', 'f', ('depth', 'time',))
    AMP3id.long_name = 'Beam 3 Echo Amplitude'
    AMP3id.units = 'counts'
    AMP3id.Type = 'scalar'
    AMP3id.transducer_offset_from_bottom = metadata['transducer_offset_from_bottom']

    Battid = rg.createVariable('Battery', 'f', ('time',))
    Battid.units = 'Volts'
    Battid.long_name = 'Battery Voltage'

    Pitchid = rg.createVariable('Pitch', 'f', ('time',))
    Pitchid.units = 'degrees'
    Pitchid.long_name = 'Instrument Pitch'

    Rollid  = rg.createVariable('Roll', 'f', ('time',))
    Rollid.units = 'degrees'
    Rollid.long_name = 'Instrument Roll'

    Headid  = rg.createVariable('Heading', 'f', ('time',))
    Headid.units = 'degrees'
    Headid.long_name = 'Instrument Heading'
    Headid.datum = 'magnetic north'

    Tmatid = rg.createVariable('TransMatrix', 'f', ('Tmatrix', 'Tmatrix',))
    Tmatid.long_name = 'Transformation Matrix for this Aquadopp'

    for n in ['1', '2']:
        if 'AnalogInput' + n in metadata:
            Anaid = rg.createVariable('AnalogInput' + n, 'f', ('time',))
            Anaid.units = 'Volts'
            Anaid.sensor_type = metadata['AnalogInput1']['sensor_type']
            Anaid.sensor_manufacturer = metadata['AnalogInput1']['sensor_manufacturer']
            Anaid.sensor_model = metadata['AnalogInput1']['sensor_model']
            Anaid.serial_number = metadata['AnalogInput1']['serial_number']

            if 'initial_sensor_height' in metadata['AnalogInput' + n]:
                # TODO
                # metadata.AnalogInput1.nominal_sensor_depth = metadata.WATER_DEPTH - metadata.AnalogInput1.initial_sensor_height;
                # netcdf.putAtt(ncid,Ana1id,'initial_sensor_height',metadata.AnalogInput1.initial_sensor_height);
                # netcdf.putAtt(ncid,Ana1id,'nominal_sensor_depth',metadata.AnalogInput1.nominal_sensor_depth);
                continue
            elif 'nominal_sensor_depth' in metadata['AnalogInput' + n]: # TODO: should be another if not elif??
                # netcdf.putAtt(ncid,Ana1id,'nominal_sensor_depth',metadata.AnalogInput1.nominal_sensor_depth);
                # metadata.AnalogInput1.initial_sensor_height = metadata.WATER_DEPTH - metadata.AnalogInput1.nominal_sensor_depth;
                # netcdf.putAtt(ncid,Ana1id,'initial_sensor_height',metadata.AnalogInput1.initial_sensor_height);
                continue

        # if isfield(metadata.AnalogInput1,'range'),
        #         netcdf.putAtt(ncid,Ana1id,'range',metadata.AnalogInput1.range);
        # end
        # if isfield(metadata.AnalogInput1.cals,'NTUcoef'),
        #         netcdf.putAtt(ncid,Ana1id,'NTUcoef',metadata.AnalogInput1.cals.NTUcoef);
        # end
        # if isfield(metadata.AnalogInput1.cals,'SEDcoef'),
        #         netcdf.putAtt(ncid,Ana1id,'SEDcoef',metadata.AnalogInput1.cals.SEDcoef);
        # end

    rg.close()

def write_metadata(rg, metadata):
    """
    Write out all metadata to CDF file
    """
    for k in metadata.keys():
        setattr(rg, k, metadata[k])

def hms2h(h,m,s):
    """
    Convert hour, minute, second to fractional hour
    """
    return h + m/60 + s/60/60

def julian(t):
    """
    Compute Julian date, relying heavily on jdcal package
    """
    y = t.year
    m = t.month
    d = t.day
    h = hms2h(t.hour, t.minute, t.second)
    return sum(jdcal.gcal2jd(y,m,d)) + h/24 + 0.5

def compute_time(RAW):
    """
    Compute Julian date and then time and time2 for use in NetCDF file
    """
    RAW['jd'] = np.array([julian(t) for t in RAW['datetime']])

    RAW['time'] = np.floor(RAW['jd'])
    # TODO: Hopefully this is correct... roundoff errors on big numbers...
    RAW['time2'] = np.round((RAW['jd'] - RAW['time'])*86400000)

    # TODO: start_time and stop_time into metadata
    # metadata['start_time'] =

    return RAW

def load_sen(RAW, basefile, metadata):
    """
    Load data from .sen file
    """
    senfile = basefile + '.sen'
    SEN = np.genfromtxt(senfile);

    RAW['heading'] = SEN[:,10]
    RAW['pitch'] = SEN[:,11]
    RAW['roll'] = SEN[:,12];
    RAW['pressure'] = SEN[:,13]
    RAW['temperature'] = SEN[:,14]
    RAW['battery'] = SEN[:,8]
    RAW['datetime'] = []
    for year, month, day, hour, minute, second in zip(SEN[:,2], SEN[:,0], SEN[:,1], SEN[:,3], SEN[:,4], SEN[:,5]):
        RAW['datetime'].append(dt.datetime(int(year), int(month), int(day), int(hour), int(minute), int(second), tzinfo=pytz.utc))
    RAW['datetime'] = np.array(RAW['datetime'])

    # Look for analog data
    if 'AnalogInput1' in metadata:
        RAW['AnaInp1'] = SEN[:,15] * 5 / 65535

    if 'AnalogInput2' in metadata:
        RAW['AnaInp2'] = SEN[:,16] * 5 / 65535

    return RAW

def load_amp_vel(RAW, basefile):
    """
    Load amplitude and velocity data from the .aN and .vN files
    """
    for n in [1, 2, 3]:
        afile = basefile + '.a' + str(n)
        RAW['AMP' + str(n)] = np.genfromtxt(afile)
        vfile = basefile + '.v' + str(n)
        RAW['Vvel' + str(n)] = np.genfromtxt(vfile)
        # convert to cm/s
        RAW['V' + str(n)] = RAW['Vvel' + str(n)] * 100

    return RAW

def read_aqd_hdr(basefile):
    """
    Get instrument metadata from .hdr file
    Was formerly readAQDprfHeader.m
    """
    #
    # % TODO read and save all of the instrument settings in the .hdr file
    # % replacing strncmp with strfind, strncmp need the exact length of string,
    # % many of those were wrong and lots of metadata was going missing.
    hdrFile = basefile + '.hdr'
    f = open(hdrFile, 'r')
    row = ''

    Instmeta = {}

    while 'Hardware configuration' not in row:
        row = f.readline().rstrip()
        if 'Profile interval' in row:
            idx = row.find(' sec')
            Instmeta['AQDProfileInterval'] = float(row[38:idx])
        elif 'Number of cells' in row:
            Instmeta['AQDNumberOfCells'] = float(row[38:])
        elif row.find('Cell size', 0, 9) != -1: # required here to differentiate from the wave cell size
            idx = row.find(' cm')
            Instmeta['AQDCellSize'] = float(row[38:idx])
        elif 'Average interval' in row:
            idx = row.find(' sec')
            Instmeta['AQDAverageInterval'] = float(row[38:idx])
        elif 'Measurement load' in row:
            idx = row.find(' %')
            Instmeta['AQDMeasurementLoad'] = float(row[38:idx])
        elif 'Transmit pulse length' in row:
            idx = row.find(' m')
            Instmeta['AQDTransmitPulseLength'] = float(row[38:idx])
        elif 'Blanking distance' in row:
            idx = row.find(' m')
            Instmeta['AQDBlankingDistance'] = float(row[38:idx])
        elif 'Compass update rate' in row:
            idx = row.find(' sec')
            Instmeta['AQDCompassUpdateRate'] = float(row[38:idx])
        elif 'Wave measurements' in row:
            Instmeta['WaveMeasurements'] = row[38:]
            # we are measuring waves
            # TODO: double check on this, the logic in the m-file seems wrong
        elif 'Analog input 1' in row:
            Instmeta['AQDAnalogInput1'] = row[38:]
        elif 'Analog input 2' in row:
            Instmeta['AQDAnalogInput2'] = row[38:]
        elif 'Power output' in row:
            Instmeta['AQDAnalogPowerOutput'] = row[38:]
        elif 'Powerlevel' in row: # TODO: WRONG, this is not analog powerlevel
            Instmeta['AQDAnalogPowerLevel'] = row[38:]
        elif 'Coordinate system' in row:
            Instmeta['AQDCoordinateSystem'] = row[38:]
        elif 'Sound speed' in row:
            Instmeta['AQDSoundSpeed'] = row[38:]
        elif 'Salinity' in row:
            Instmeta['AQDSalinity'] = row[38:]
        elif 'Number of beams' in row:
            Instmeta['AQDNumberOfBeams'] = float(row[38:])
        elif 'Number of pings per burst' in row:
            Instmeta['AQDNumberOfPingsPerBurst'] = float(row[38:])
        elif 'Software version' in row:
            Instmeta['AQDSoftwareVersion'] = row[38:] # can't be float, was wrong in m-file
        elif 'Deployment name' in row:
            Instmeta['AQDDeploymentName'] = row[38:]
        elif 'Deployment time' in row:
            Instmeta['AQDDeploymentTime'] = row[38:]
        elif 'Comments' in row:
            Instmeta['AQDComments'] = row[38:]

    while 'Head configuration' not in row:
        row = f.readline().rstrip()
        if 'Serial number' in row:
            Instmeta['AQDSerial_Number'] = row[38:]
        elif 'Hardware revision' in row:
            Instmeta['AQDHardwareRevision'] = row[38:]
        elif 'Revision number' in row:
            Instmeta['AQDRevisionNumber'] = row[38:]
        elif 'Recorder size' in row:
            Instmeta['AQDRecorderSize'] = row[38:]
        elif 'Firmware version' in row:
            Instmeta['AQDFirmwareVersion'] = row[38:]
        elif 'Velocity range' in row:
            Instmeta['AQDVelocityRange'] = row[38:]
        elif 'Power output' in row:
            Instmeta['AQDAnalogPowerOutput'] = row[38:]
        elif 'Analog input #1 calibration (a0, a1)' in row:
            Instmeta['AQDAnalogInputCal1'] = float(row[38:])
        elif 'Analog input #2 calibration (a0, a1)' in row:
            Instmeta['AQDAnalogInputCal2'] = float(row[38:])
        elif 'Sync signal data out delay' in row:
            Instmeta['AQDSyncOutDelay'] = row[38:]
        elif 'Sync signal power down delay' in row:
            Instmeta['AQDSyncPowerDelay'] = row[38:]

    while 'Current profile cell center distance from head (m)' not in row:
    # while 'Data file format' not in row:
        row = f.readline().rstrip()
        if 'Pressure sensor' in row:
            Instmeta['AQDPressureSensor'] = row[38:]
        elif 'Compass' in row:
            Instmeta['AQDCompass'] = row[38:]
        elif 'Tilt sensor' in row:
            Instmeta['AQDTilt'] = row[38:]
        elif 'Head frequency' in row:
            idx = row.find(' kHz')
            Instmeta['AQDFrequency'] = float(row[38:idx])
        elif 'Number of beams' in row:
            Instmeta['AQDNumBeams'] = float(row[38:])
        elif 'Serial number' in row:
            Instmeta['AQDHeadSerialNumber'] = row[38:]
        elif 'Transformation matrix' in row:
            Instmeta['AQDTransMatrix'] = np.zeros((3,3))
            for n in np.arange(3):
                Instmeta['AQDTransMatrix'][n,:] = [float(x) for x in row[38:].split()]
                row = f.readline().rstrip()
        elif 'Pressure sensor calibration' in row:
            Instmeta['AQDPressureCal'] = row[38:]

    #          elseif (strfind(str,'Wave measurements'))
    #              Instmeta.WaveMeasurements = (str(39:end));
    #             if strcmp(Instmeta.WaveMeasurements,'ENABLED')
    #              elseif (strfind(str,'Wave - Powerlevel'))
    #                  Instmeta.WavePower = (str(39:is-2));
    #              elseif (strfind(str,'Wave - Interval'))
    #                  is=findstr(str,'sec');
    #                  Instmeta.WaveInterval = str2num(str(39:is-2));
    #              elseif (strfind(str,'Wave - Number of samples'))
    #                  Instmeta.WaveNumberOfSamples = str2num(str(39:42));
    #              elseif (strfind(str,'Wave - Sampling rate'))
    #                  Instmeta.WaveSampleRate = str(39:42);
    #              elseif (strfind(str,'Wave - Cell size'))
    #                  Instmeta.WaveCellSize = str(39:42);
    #             end


    # % infer some things based on the Aquadopp brochure
    if Instmeta['AQDFrequency'] == 400:
        Instmeta['AQDBeamWidth'] = 3.7
    elif Instmeta['AQDFrequency'] == 600:
        Instmeta['AQDBeamWidth'] = 3.0
    elif Instmeta['AQDFrequency'] == 1000:
        Instmeta['AQDBeamWidth'] = 3.4
    elif Instmeta['AQDFrequency'] == 2000:
        Instmeta['AQDBeamWidth'] = 1.7
    else:
        Instmeta['AQDBeamWidth'] = np.nan

    Instmeta['AQDBeamPattern'] = 'convex'
    Instmeta['AQDBeamAngle'] = 25;
# Instmeta.AQDVelRange = 1000; % cm/s
# Instmeta.AQDTempRange = [-4 40];
# Instmeta.AQDPressRange = [0 100];
# % no tilt range given in AQD docs
#
# fclose(hdr);
#
# % %if waves data were not collected remove wave parameters from metadata
# % if strfind(Instmeta.AQDWaveStatus,'DISABLED',7)
# %     fields = {'AQDWavePower','AQDWaveInterval','AQDWaveSampleRate','AQDWaveNumberOfSamples'};
# %     Instmeta = rmfield(Instmeta,fields);
# % else
# % end

    return Instmeta
