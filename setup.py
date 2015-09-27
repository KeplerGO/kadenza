#!/usr/bin/env python
from setuptools import setup

# Command-line tools
entry_points = {'console_scripts': [
    'kadenza-tpf = kadenza.kadenza:kadenza_tpf_main'
]}

setup(name='kadenza',
      version='1.0.0',
      description="Converts raw cadence data from the Kepler spacecraft "
                  "into astronomer-friendly FITS files.",
      long_description=open('README.md').read(),
      author='Geert Barentsen',
      author_email='hello@geert.io',
      license='MIT',
      packages=['kadenza'],
      install_requires=['numpy',
                        'astropy'],
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
