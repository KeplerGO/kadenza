# Kadenza: Kepler/K2 Cadence Data Reader 
***Converts raw cadence data from the Kepler spacecraft into astronomer-friendly FITS files.***

[![MIT license](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/barentsen/k2flix/blob/master/LICENSE) [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

Kadenza allows the raw cadence pixel data from the Kepler space telescope
to be inspected as soon as they arrive on Earth;
bypassing the (time-consuming) standard pipeline processing
to enable time-critical data analyses to be carried out quickly.

The primary motivation for creating this tool is to
enable the K2 Campaign 9 microlensing science team to use the raw
spacecraft data to identify microlensing events as soon as possible.
It also enables access to the dithering data which was obtained during
the original comissioning of the Kepler mission.
Most users of Kepler/K2 data will not require this tool however,
because it does ***not*** produce properly calibrated files.

Kadenza can be used both as a command-line tool or using its Python API.

### Example
Convert the raw cadence data into TPF files:
```
$ kadenza-tpf --target 200049143 cadence-data-list.txt pixel-mapping-file.fits
```

Produce a sparse FFI frame for a given channel:
```
$ kadenza-ffi cadence-data-file.fits pixel-mapping-file.fits
``` 

### Installation
If you have a working installation of Python on your system,
you can download and install the latest version as follows:
```
$ git clone https://github.com/KeplerGO/kadenza.git
$ cd kadenza
$ python setup.py install
```
Kadenza has only been tested under Linux with Python 3.4 at present.
Support for legacy Python will be added soon.

### Using kadenza
```
$ kadenza-tpf --help
usage: kadenza-tpf [-h] [-t [target_id]] cadencefile_list pixelmap_file

Turn raw Kepler Cadence Data into uncalibrated Target Pixel Files (TPF).

positional arguments:
  cadencefile_list      path to a text file that lists the '*_lcs-targ.fits'
                        cadence data files to use
  pixelmap_file         path to the '*_lcm.fits' pixel mapping reference file

optional arguments:
  -h, --help            show this help message and exit
  -t [target_id], --target [target_id]
                        only produce a TPF file for a specific EPIC/KIC
                        target_id
```

```
$ kadenza-ffi --help
usage: kadenza-ffi [-h] cadence_file pixelmap_file

Turns a raw Kepler Cadence Data file into an uncalibrated Full Frame Image
(FFI).

positional arguments:
  cadence_file   path to the '*_lcs-targ.fits' cadence data file
  pixelmap_file  path to the '*_lcm.fits' pixel mapping reference file

optional arguments:
  -h, --help     show this help message and exit
```

### Caveats

There are two main caveats to be aware of:
 * The present version does not set all header keywords exactly as they would
   appear in an official product.  In particular, the WCS keywords are
   untested.
 * This tool does *not* calibrate the data at all, it merely serves to
   transform the raw pixel counts into a FITS format that is similar to,
   but not at all identical, to the official pipeline products.

### Contributing
To report bugs and request features, please use the [issue tracker](https://github.com/KeplerGO/kadenza/issues) or open a pull request.

### License
Kadenza is made available under the MIT License.
For details see the LICENSE file.
