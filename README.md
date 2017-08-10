# aqdlib - Process Nortek Aquadopp data in Python

This module contains code to process Nortek Aquadopp Profiler data as exported from AquaPro in text format. It consists of two main steps:

1. Convert from text to a raw netCDF file with .cdf extension
2. Convert the raw .cdf data into an EPIC-compliant netCDF file with .nc extension

## Text to raw netCDF (.cdf)

This step will generally be completed via the runhdr2cdf.py script, which relies heavily on this module.

## Raw netCDF (.cdf) to EPIC-compliant and processed netCDF (.nc)

This step will generally be completed via the runcdf2nc.py script, which also relies heavily on this module.