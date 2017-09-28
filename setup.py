from setuptools import setup

setup(name='aqdlib',
      version='0.0',
      description='Process Nortek Aquadopp data in Python',
      author='Dan Nowacki',
      author_email='dnowacki@usgs.gov',
      url='https://github.com/dnowacki-usgs/aqdlib',
      install_requires=['numpy', 'netCDF4', 'xarray'],
      packages=['aqdlib'],
     )
