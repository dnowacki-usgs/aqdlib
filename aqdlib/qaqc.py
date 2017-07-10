from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from dateutil import parser
import pytz


def plot_inwater(RAW, inwater_time, outwater_time):
    plt.figure(figsize=(12,8))
    plt.plot(RAW['datetime'], RAW['pressure'])
    stind = np.where(RAW['datetime'] > pytz.utc.localize(parser.parse(inwater_time)))[0][0]
    stend = np.where(RAW['datetime'] < pytz.utc.localize(parser.parse(outwater_time)))[0][-1]
    for s in [stind, stend]:
        plt.plot(RAW['datetime'][s], RAW['pressure'][s], 'r.')
    plt.title(str(stind) + ' - ' + str(stend))
    plt.show()
