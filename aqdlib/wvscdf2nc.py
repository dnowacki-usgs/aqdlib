from __future__ import division, print_function
import numpy as np
import datetime as dt
import pytz
from .aqdhdr2cdf import compute_time, read_aqd_hdr, check_metadata, check_orientation, define_aqd_cdf_file, write_aqd_cdf_data
import pandas as pd
import aqdcdf2nc

def wvs_cdf_to_nc(cdf_filename, metadata, p_1ac=False):

    nc_filename = metadata['filename'] + '-wvsb-cal.nc' # TODO: why is a "b" in there?

    VEL = {}

    VEL, INFO = aqdcdf2nc.load_cdf_amp_vel(cdf_filename, VEL, metadata, p_1ac=p_1ac)
    #
    # define_aqd_nc_file(nc_filename, VEL, metadata, INFO)
    #
    # write_aqd_nc_file(nc_filename, VEL, metadata)
    #
    # # TODO: Need to add all global attributes from CDF to NC file (or similar)
    # qaqc.add_min_max(nc_filename)
    # print('Added min/max values')
    #
    # qaqc.add_final_metadata(nc_filename)
    # print('Added final metadata')
    #
    # print('Done writing NetCDF file', nc_filename)

    return VEL
