import os
import numpy as np

def load(basefile, metadata):

    instmeta = readAQDprfHeader(basefile)

    doublefill = 1e35
    shortfill = -32768
    print "Loading ASCII files"

    # basePath = os.path.dirname(basefile)
    # baseName = os.path.basename(basefile)
    RAW = {}
    # Load sensor data
    RAW = load_sen(RAW, basefile)
    # Load amplitude and velocity data
    RAW = load_amp_vel(RAW, basefile)

    # update metadata from Aquadopp header to CMG standard so that various
    # profilers have the same attribute wording.  Redundant, but necessary
    metadata['bin_count'] = instmeta['AQDNumberOfCells']
    metadata['bin_size'] = instmeta['AQDCellSize']/100 # from cm to m
    metadata['blanking_distance'] = instmeta['AQDBlankingDistance'] # already in m
    # Nortek lists the distance to the center of the first bin as the blanking
    # distance plus one cell size
    metadata['center_first_bin'] = metadata['blanking_distance'] + metadata['bin_size'] # in m

    print metadata['orientation']
    print 'Center_first_bin = %f\n' % metadata['center_first_bin']
    print 'bin_size = %f\n' % metadata['bin_size']
    print 'bin_count = %f\n' % metadata['bin_count']
    RAW['bindist'] = np.arange(metadata['center_first_bin'], (metadata['center_first_bin']+((metadata['bin_count']-1)*metadata['bin_size'])), metadata['bin_size'])

    if metadata['orientation'] is 'UP':
        print 'User instructed that instrument was pointing UP'
        # depth, or distance below surface, is a positive number below the
        # surface, negative above the surface, for CMG purposes and consistency with ADCP
        RAW['Depths'] = (metadata['WATER_DEPTH'] - metadata['transducer_offset_from_bottom']) - RAW['bindist']
        Depth_NOTE = 'user reports uplooking bin depths = water_depth - transducer offset from bottom - bindist'
    elif metadata['orientation'] is 'DOWN':
        print 'User instructed that instrument was pointing DOWN'
        RAW['Depths'] = (metadata['WATER_DEPTH'] - metadata['transducer_offset_from_bottom']) + RAW['bindist']
        Depth_NOTE = 'user reports downlooking bin depths = water_depth - transducer_offset_from_bottom + bindist'

    RAW['instmeta'] = instmeta

    return RAW
    # jd = julian([SEN(:,3) SEN(:,1) SEN(:,2) SEN(:,4) SEN(:,5) SEN(:,6)]);

def load_sen(RAW, basefile):
    senfile = basefile + '.sen'
    SEN = np.genfromtxt(senfile);

    RAW['heading'] = SEN[:,10]
    RAW['pitch'] = SEN[:,11]
    RAW['roll'] = SEN[:,12];
    RAW['pressure'] = SEN[:,13]
    RAW['temperature'] = SEN[:,14]
    RAW['battery'] = SEN[:,8]

    return RAW

def load_amp_vel(RAW, basefile):
    for n in [1, 2, 3]:
        afile = basefile + '.a' + str(n)
        RAW['AMP' + str(n)] = np.genfromtxt(afile)
        vfile = basefile + '.v' + str(n)
        RAW['Vvel' + str(n)] = np.genfromtxt(vfile)
        # convert to cm/s
        RAW['V' + str(n)] = RAW['Vvel' + str(n)] * 100

    return RAW

def readAQDprfHeader(basefile):
    # %function to get instrument metadata about AWAC deployment from header file
    #
    # % TODO read and save all of the instrument settings in the .hdr file
    # % replacing strncmp with strfind, strncmp need the exact length of string,
    # % many of those were wrong and lots of metadata was going missing.
    hdrFile = basefile + '.hdr'
    with open(hdrFile) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    Instmeta = {}
    for row in content:
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

    return Instmeta

    # for k, v in Instmeta.iteritems():
    #     print k, v


    #          elseif (strfind(str,'Blanking distance'))
    #              is=findstr(str,'m');
    #              Instmeta.AQDBlankingDistance = str2num(str(39:is-2));
    #          elseif (strfind(str,'Compass update rate'))
    #              is=findstr(str,'sec');
    #              Instmeta.AQDCompassUpdateRate = str2num(str(39:is-2));
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
    #          elseif (strfind(str,'Analog input 1'))
    #              Instmeta.AQDAnalogInput1 = str(39:end);
    #          elseif (strfind(str,'Analog input 2'))
    #              Instmeta.AQDAnalogInput2 = str(39:end);
    #          elseif (strfind(str,'Power output'))
    #              Instmeta.AQDAnalogPowerOutput = str(39:end);
    #          elseif (strfind(str,'Powerlevel'))
    #              Instmeta.AQDAnalogPowerLevel = str(39:end);
    #          elseif (strfind(str,'Coordinate system'))
    #              Instmeta.AQDCoordinateSystem = str(39:end);
    #          elseif (strfind(str,'Sound speed'))
    #              Instmeta.AQDSoundSpeed = str(39:end);
    #          elseif (strfind(str,'Salinity'))
    #              Instmeta.AQDSalinity = str(39:end);
    #          elseif (strfind(str,'Number of beams'))
    #              Instmeta.AQDNumberOfBeams = str2num(str(39:end));
    #          elseif (strfind(str,'Number of pings per burst'))
    #              Instmeta.AQDNumberOfPingsPerBurst = str2num(str(39:end));
    #          elseif (strfind(str,'Software version'))
    #              Instmeta.AQDSoftwareVersion = str2num(str(39:end));
    #          elseif (strfind(str,'Deployment name'))
    #              Instmeta.AQDDeploymentName = str(39:end);
    #          elseif (strfind(str,'Deployment time'))
    #              Instmeta.AQDDeploymentTime = str(39:end);
    #          elseif (strfind(str,'Comments'))
    #              Instmeta.AQDComments = str(39:end);
    #          end
    # end
# % getting the hardware configuration
# while isempty(strfind(str,'Head configuration'));
#      str=fgetl(hdr);
#          if (strfind(str,'Compass'))
#              Instmeta.AQDCompass = str(39:end);
#          elseif (strfind(str,'Tilt sensor'))
#              Instmeta.AQDTilt = str(39:end);
#          elseif (strfind(str,'Pressure sensor calibration'))
#              Instmeta.AQDPressureCal = str2num(str(39:end));
#          elseif (strfind(str,'Number of beams'))
#              Instmeta.AQDNumBeams = str2num(str(39:end));
#          elseif (strfind(str,'Serial number'))
#              %Instmeta.SERIAL_NUMBER = str(39:end);
#              Instmeta.AQDSerial_Number = str(39:end);
#          elseif (strfind(str,'Internal code version'))
#              Instmeta.AQDInternalCodeVersion = str(39:end);
#          elseif (strfind(str,'Revision number'))
#              Instmeta.AQDRevisionNumber = str(39:end);
#          elseif (strfind(str,'Recorder size'))
#              Instmeta.AQDRecorderSize = str(39:end);
#          elseif (strfind(str,'Firmware version'))
#              Instmeta.AQDFirmwareVersion = str(39:end);
#          elseif (strfind(str,'Power output'))
#              Instmeta.AQDAnalogPowerOutput = str(39:end);
#          elseif (strfind(str,'Analog input #1 calibration (a0, a1)'))
#              Instmeta.AQDAnalogInputCal1 = str2num(str(39:end));
#          elseif (strfind(str,'Analog input #2 calibration (a0, a1)'))
#              Instmeta.AQDAnalogInputCal2 = str2num(str(39:end));
#          elseif (strfind(str,'Sync signal data out delay'))
#              Instmeta.AQDSyncOutDelay = str(39:end);
#          elseif (strfind(str,'Sync signal power down delay'))
#              Instmeta.AQDSyncPowerDelay = str(39:end);
#          end
# end
# % getting the head configuration
# while isempty(strfind(str,'Data file format'));
#      str=fgetl(hdr);
#          if (strfind(str,'Head frequency'))
#              is=findstr(str,'kHz');
#              Instmeta.AQDFrequency = str2num(str(39:is-2));
#          elseif (strfind(str,'Serial number'))
#              Instmeta.AQDHeadSerialNumber = str(39:end);
#          elseif (strfind(str,'Transformation matrix'))
#              Instmeta.AQDTransMatrix = zeros(3,3);
#              Instmeta.AQDTransMatrix(1,:) = strread(str(39:end));
#              str=fgetl(hdr);
#              Instmeta.AQDTransMatrix(2,:) = strread(str(39:end));
#              str=fgetl(hdr);
#              Instmeta.AQDTransMatrix(3,:) = strread(str(39:end));
#          elseif (strfind(str,'Pressure sensor calibration'))
#              Instmeta.AQDPressureCal = str2num(str(39:end));
#          end
# end
#
# % infer some things based on the Aquadopp brochure
# switch Instmeta.AQDFrequency
#     case 400
#         Instmeta.AQDBeamWidth = 3.7;
#     case 600
#         Instmeta.AQDBeamWidth = 3.0;
#     case 1000
#         Instmeta.AQDBeamWidth = 3.4;
#     case 2000
#         Instmeta.AQDBeamWidth = 1.7;
#     otherwise
#         Instmeta.AQDBeamWidth = NaN;
# end
# Instmeta.AQDBeamPattern = 'convex';
# Instmeta.AQDBeamAngle = 25;
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
