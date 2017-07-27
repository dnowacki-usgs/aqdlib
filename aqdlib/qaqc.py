from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
from dateutil import parser
import pytz
import calendar
import jdcal

def plot_inwater(RAW, inwater_time, outwater_time):
    plt.figure(figsize=(12,8))
    plt.plot(RAW['datetime'], RAW['pressure'])
    stind = np.where(RAW['datetime'] > pytz.utc.localize(parser.parse(inwater_time)))[0][0]
    stend = np.where(RAW['datetime'] < pytz.utc.localize(parser.parse(outwater_time)))[0][-1]
    for s in [stind, stend]:
        plt.plot(RAW['datetime'][s], RAW['pressure'][s], 'r.')
    plt.title(str(stind) + ' - ' + str(stend))
    plt.show()

def tots(d):
    return calendar.timegm(d.timetuple())

def atmcomp(aqddatetime, aqdpress, metdatetime, metpress, offset=0):
    """
    Atmospheric pressure compensation of pressure records.
    Inputs:
        aqddatetime: array of datetimes from Aquadopp [decibars]
        aqdpress: array of pressure (depth) values from Aquadopp

        metdatetime: array of datetimes from met station
        metpress: array of pressure (atmospheric) *** NOTE [decibars] ***

        offset: offset for when sensor is out of water

    Outputs:
        aqdpress_ac: atmospherically corrected pressure record
    """

    ptime = np.array([tots(x) for x in aqddatetime])
    mtime = np.array([tots(x) for x in metdatetime])

    metinterp = np.interp(ptime, mtime, metpress)

    aqdpress_ac = aqdpress - (metinterp - offset)

    return aqdpress_ac

def plot_atmcomp(aqddatetime, aqdpress, aqdpress_ac):
    plt.figure()
    plt.plot(aqddatetime, aqdpress)
    plt.plot(aqddatetime, aqdpress_ac)
    plt.show()

def coord_transform(vel1, vel2, vel3, heading, pitch, roll, T, cs):
    N, M = np.shape(vel1)

    u = np.zeros((N,M))
    v = np.zeros((N,M))
    w = np.zeros((N,M))

    if cs == 'ENU':
        print('Data already in earth coordinates; applying magnetic correction')
        u =  vel1 * math.cos(magvar) + vel2 * math.sin(magvar);
        v = -vel1 * math.sin(magvar) + vel2 * math.cos(magvar);
        w = vel3;
    elif cs == 'XYZ':
        # TODO: add XYZ
        print("xyz")
    elif cs == 'BEAM':
        print('Data are in BEAM coordinates; transforming to earth coordinates...')

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

    return (u, v, w)
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
