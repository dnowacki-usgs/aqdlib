from .aqdhdr2cdf import prf_to_cdf
from .aqdcdf2nc import cdf_to_nc
from .globalatts import read_globalatts
from .qaqc import plot_inwater, coord_transform
from .atmcomp import atmcomp, plot_atmcomp

DOUBLE_FILL = 1e35
SHORT_FILL = -32768
