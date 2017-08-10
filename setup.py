from setuptools import setup

setup(name='aqdlib',
      version='0.0',
      description='Aquadopp Processing Utilities',
      author='Dan Nowacki',
      author_email='dnowacki@usgs.gov',
      url='https://github.com/dnowacki-usgs/aqdlib',
      install_requires=['numpy', 'netCDF4'],
      scripts=['scripts/runaqd2cdf.py', 'scripts/runcdf2nc.py'],
      packages=['aqdlib'],
     )
