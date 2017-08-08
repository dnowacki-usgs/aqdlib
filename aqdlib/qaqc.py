from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from dateutil import parser
import pytz
import jdcal
import netCDF4 as nc
import aqdlib
import datetime as dt

def time_time2_to_datetime(time, time2):
    times = []

    for t, t2 in zip(time, time2):
        year, mon, day, frac = jdcal.jd2gcal(t - 0.5, t2/86400000)
        hour, minute, second = day2hms(frac)
        times.append(dt.datetime(year, mon, day, hour, minute, second, tzinfo=pytz.utc))

    return np.array(times)

def load_cdf(cdf_filename, varis):
    try:
        rg = nc.Dataset(cdf_filename, 'r')
        RAW = {}
        for var in varis:
            RAW[var] = rg[var][:]

        return RAW

    finally:
        rg.close()

def add_min_max(cdf_filename):
    """
    Add minimum and maximum values to variables in NC or CDF files
    """
    try:
        rg = nc.Dataset(cdf_filename, 'r+')
        exclude = rg.dimensions.keys()
        exclude.extend(('time2', 'TIM'))
        for var in rg.variables:
            if var not in exclude:
                rg[var].minimum = np.nanmin(rg[var][:])
                rg[var].maximum = np.nanmax(rg[var][:])

    finally:
        rg.close()

def add_fill_values(cdf_filename):
    """
    Assign _FillValue to each variable
    """
    try:
        rg = nc.Dataset(cdf_filename, 'r+')
        exclude = rg.dimensions.keys()
        exclude.extend(('time2', 'TIM'))
        for var in rg.variables:
            if var not in exclude:
                print(var, rg[var].dtype)
                if rg[var].dtype == 'float32':
                    rg[var]._FillValue = aqdlib.DOUBLE_FILL

    finally:
        rg.close()

def plot_inwater(RAW, inwater_time, outwater_time):
    plt.figure(figsize=(12,8))
    plt.plot(RAW['datetime'], RAW['pressure'])
    stind = np.where(RAW['datetime'] > pytz.utc.localize(parser.parse(inwater_time)))[0][0]
    stend = np.where(RAW['datetime'] < pytz.utc.localize(parser.parse(outwater_time)))[0][-1]
    for s in [stind, stend]:
        plt.plot(RAW['datetime'][s], RAW['pressure'][s], 'r.')
    plt.title(str(stind) + ' - ' + str(stend))
    plt.show()

def coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, cs):
    N, M = np.shape(vel1)

    u = np.zeros((N,M))
    v = np.zeros((N,M))
    w = np.zeros((N,M))

    if cs == 'ENU':
        print('Data already in Earth coordinates; doing nothing')

        u = vel1
        v = vel2
        w = vel3
        # u =  vel1 * math.cos(magvar) + vel2 * math.sin(magvar);
        # v = -vel1 * math.sin(magvar) + vel2 * math.cos(magvar);
        # w = vel3;
    elif cs == 'XYZ':
        # TODO: add XYZ
        print("Data are in XYZ coordinates; transforming to Earth coordinates")
    elif cs == 'BEAM':
        print('Data are in BEAM coordinates; transforming to Earth coordinates')

        for i in range(N):
            hh = np.pi * (heading[i] - 90) / 180;
            pp = np.pi * pitch[i] / 180;
            rr = np.pi * roll[i] / 180;

            H = np.array([[ np.cos(hh), np.sin(hh), 0],
                          [-np.sin(hh), np.cos(hh), 0],
                          [ 0,          0,          1]])

            # make tilt matrix
            P = np.array([[np.cos(pp), -np.sin(pp) * np.sin(rr), -np.cos(rr) * np.sin(pp)],
                          [0,           np.cos(rr),              -np.sin(rr)],
                          [np.sin(pp),  np.sin(rr) * np.cos(pp),  np.cos(pp) * np.cos(rr)]])

            # resulting transformation matrix
            R = np.dot(np.dot(H, P), T)

            for j in range(M):
                vel = np.dot(R, np.array([vel1[i,j], vel2[i,j], vel3[i,j]]).T)
                u[i,j] = vel[0]
                v[i,j] = vel[1]
                w[i,j] = vel[2]

    return u, v, w

def set_orientation(VEL, T, metadata, INFO):
    # TODO: this code seems too complicated. also should we really be modifying the trans matrix?
    # TODO: deal with atmos pressure

    # print(INFO)

    N, M = np.shape(VEL['U'])

    if 'press_ac' in VEL:
        Wdepth = np.nanmean(VEL['press_ac']) + INFO['transducer_offset_from_bottom']
    else:
        Wdepth = np.nanmean(VEL['pressure']) + INFO['transducer_offset_from_bottom']

    blank2 = INFO['AQDBlankingDistance'] + INFO['transducer_offset_from_bottom']
    binn = INFO['bin_size']
    blank3 = INFO['transducer_offset_from_bottom'] - INFO['AQDBlankingDistance']
    binc = INFO['bin_count']

    # print(Wdepth, binn, blank2)

    if INFO['orientation'] == 'UP':
        print('User instructed that instrument was pointing UP')
        # VEL['depths'] = np.flipud(np.arange(Wdepth - (binn * (M - 1) + blank2 + binn), Wdepth - (blank2 + binn), binn)) # need to use flipud because 1d array
        VEL['depths'] = np.flipud(np.linspace(Wdepth - (binn * (M - 1) + blank2 + binn), Wdepth - (blank2 + binn), num=binc))
        # print(VEL['depths'])
    elif INFO['orientation'] == 'DOWN':
        print('User instructed that instrument was pointing DOWN')
        T[1,:] = -T[1,:]
        T[2,:] = -T[2,:]
        # VEL['depths'] = np.arange(Wdepth - blank3 + binn, Wdepth - blank3 + binn * M, binn)
        VEL['depths'] = np.linspace(Wdepth - blank3 + binn, Wdepth - blank3 + binn * M, num=binc)

    return VEL, T

def make_bin_depth(VEL, metadata):

    N, M = np.shape(VEL['U'])

    if 'press_ac' in VEL:
        VEL['bin_depth'] = np.tile(VEL['press_ac'], (M, 1)) - np.tile(VEL['bindist'], (N, 1)).T;
    else:
        VEL['bin_depth'] = np.tile(VEL['pressure'], (M, 1)) - np.tile(VEL['bindist'], (N, 1)).T;

    return VEL

def magvar_correct(VEL, metadata):

    if 'magnetic_variation_at_site' in metadata:
        magvardeg = metadata['magnetic_variation_at_site']
    elif 'magnetic_variation' in metadata:
        magvardeg = metadata['magnetic_variation']
    else:
        print('No magnetic variation information provided; using zero for compass correction')
        magvardeg = 0

    print('Rotating heading by %f degrees' % magvardeg)

    VEL['heading'] = VEL['heading'] + magvardeg
    VEL['heading'][VEL['heading'] >= 360] = VEL['heading'][VEL['heading'] >= 360] - 360
    VEL['heading'][VEL['heading'] < 0] = VEL['heading'][VEL['heading'] < 0] + 360

    vel1 = VEL['U'].copy()
    vel2 = VEL['V'].copy()

    magvar = magvardeg * np.pi / 180

    print('Rotating horizontal velocities by %f degrees' % magvardeg)

    VEL['U'] =  vel1 * np.cos(magvar) + vel2 * np.sin(magvar)
    VEL['V'] = -vel1 * np.sin(magvar) + vel2 * np.cos(magvar)

    return VEL

def create_water_depth(VEL, metadata):
    if 'initial_instrument_height' in metadata:
        if 'press_ac' in VEL:
            metadata['nominal_instrument_depth'] = np.nanmean(VEL['press_ac'])
            VEL['Depth'] = metadata['nominal_instrument_depth']
            wdepth = metadata['nominal_instrument_depth'] + metadata['initial_instrument_height']
            metadata['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor, atmospherically corrected'
            metadata['WATER_DEPTH_datum'] = 'MSL'
        elif 'press' in VEL:
            metadata['nominal_instrument_depth'] = np.nanmean(VEL['press'])
            VEL['Depth'] = metadata['nominal_instrument_depth']
            wdepth = metadata['nominal_instrument_depth'] + metadata['initial_instrument_height']
            metadata['WATER_DEPTH_source'] = 'water depth = MSL from pressure sensor'
            metadata['WATER_DEPTH_datum'] = 'MSL'
        else:
            wdepth = metadata['WATER_DEPTH']
            metadata['nominal_instrument_depth'] = metadata['WATER_DEPTH'] - metadata['initial_instrument_height']
            VEL['Depth'] = metadata['nominal_instrument_depth']
        metadata['WATER_DEPTH'] = wdepth # TODO: why is this being redefined here? Seems redundant
    elif 'nominal_instrument_depth' in metadata:
        metadata['initial_instrument_height'] = metadata['WATER_DEPTH'] - metadata['nominal_instrument_depth']
        VEL['Depth'] = metadata['nominal_instrument_depth']

    if 'initial_instrument_height' not in metadata:
        metadata['initial_instrument_height'] = 0 # TODO: do we really want to set to zero?

    return VEL, metadata

def trim_vel(VEL, metadata, INFO):
    N, M = np.shape(VEL['U'])

    if 'press_ac' in VEL:
        WL = VEL['press_ac'] + INFO['transducer_offset_from_bottom']
    else:
        WL = VEL['pressure'] + INFO['transducer_offset_from_bottom']


    if 'trim_method' in metadata:
        blank = INFO['AQDBlankingDistance'] + metadata['transducer_offset_from_bottom'] # TODO: check this logic
        binn = INFO['bin_size']
        binc = INFO['bin_count']
        if metadata['trim_method'].lower() == 'water level' or metadata['trim_method'].lower() == 'water level sl':
            print('User instructed to trim data at the surface using pressure data')
            if metadata['trim_method'].lower() == 'water level':
                print('Trimming using water level')
                dist2 = np.linspace(blank + binn/2, (binn * (M-1)) + blank + binn/2, num=binc)
            elif metadata['trim_method'].lower() == 'water level sl':
                print('Trimming using water level and sidelobes')
                dist2 = np.linspace(blank + binn, (binn * (M-1)) + blank + binn, num=binc)
            # TODO: Why does matlab incorporate cosd of the beam angle??

            d2 = np.tile(dist2, (np.shape(WL)[0], 1))
            WL2 = np.tile(WL, (np.shape(d2)[1], 1)).T

            if metadata['trim_method'].lower() == 'water level':
                bads = d2 >= WL2
            elif metadata['trim_method'].lower() == 'water level sl':
                bads = d2 >= WL2 * np.cos(np.deg2rad(INFO['AQDBeamAngle']))

            VEL['U'][bads] = np.nan
            VEL['V'][bads] = np.nan
            VEL['W'][bads] = np.nan
            VEL['AGC'][bads] = np.nan

            # find first bin that is all bad values
            lastbin = np.argmin(np.all(bads, axis=0) == False)

            VEL['U'] = VEL['U'][:,0:lastbin]
            VEL['V'] = VEL['V'][:,0:lastbin]
            VEL['W'] = VEL['W'][:,0:lastbin]
            VEL['AGC'] = VEL['AGC'][:,0:lastbin]
            VEL['depths'] = VEL['depths'][0:lastbin]
            VEL['bindist'] = VEL['bindist'][0:lastbin]

            # TODO: need to add histcomment

    return VEL

def day2hms(d):
    frachours = d * 24
    h = np.int(np.floor(frachours))
    fracmins = (frachours - h) * 60
    m = np.int(np.floor(fracmins))
    fracsecs = (fracmins - m) * 60
    s = np.int(fracsecs)

    return h, m, s

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
