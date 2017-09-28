# aqdlib - Process Nortek Aquadopp data in Python

This module contains code to process Nortek Aquadopp Profiler data as exported from AquaPro in text format, consistent with the procedures of the Sediment Transport Group at the USGS Woods Hole Coastal and Marine Science Center.

Processing consists of two main steps:

1. Convert from text to a raw netCDF file with `.cdf` extension
2. Convert the raw `.cdf` data into an EPIC-compliant netCDF file with `.nc` extension, optionally including atmospheric correction of the pressure data

## Text to raw netCDF (.cdf)

This step will generally be completed by using the import statement `from aqdlib import aqdhdr2cdf` and calling `aqdhdr2cdf.hdr_to_cdf()`, or by running `aqdhdr2cdf.py` from the command line.`

## Raw netCDF (.cdf) to EPIC-compliant and processed netCDF (.nc)

This step will generally be completed by using the import statement `from aqdlib import aqdcdf2nc` and calling `aqdcdf2nc.cdf_to_nc()`, or by running `aqdcdf2nc.py` from the command line. When calling `cdf_to_nc()`, the user may provide the path to a netCDF file consisting of atmospheric pressure, which will be used to atmospherically correct the pressure data. This path can also be passed as a command-line argument to `aqdcdf2nc.py`
