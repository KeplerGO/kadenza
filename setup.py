#!/usr/bin/env python
import os
import sys
from setuptools import setup

if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()

# Load the __version__ variable without importing the package already
exec(open('kadenza/version.py').read())

# Command-line tools
entry_points = {'console_scripts': [
    'kadenza-tpf = kadenza.kadenza:kadenza_tpf_main',
    'kadenza-ffi = kadenza.kadenza:kadenza_ffi_main'
]}

setup(name='kadenza',
      version=__version__,
      description="Converts raw cadence data from the Kepler spacecraft "
                  "into astronomer-friendly FITS files.",
      long_description=open('README.md').read(),
      author='Geert Barentsen',
      author_email='hello@geert.io',
      license='MIT',
      packages=['kadenza'],
      install_requires=['numpy',
                        'astropy>=1.1'],
      entry_points=entry_points,
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          ],
      )
