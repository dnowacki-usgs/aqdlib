from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from dateutil import parser
import pytz
import calendar

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
