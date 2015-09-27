# Kadenza: Kepler/K2 Cadence Data Reader 
***Converts raw cadence data from the Kepler spacecraft into astronomer-friendly FITS files.***

[![MIT license](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/barentsen/k2flix/blob/master/LICENSE) [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)

Kadenza allows the raw cadence pixel data from the Kepler space telescope
to be inspected as soon as it arrives on Earth;
bypassing the (somewhat time-consuming) standard pipeline processing
to enable time-critical data analyses to be carried out.

The primary motivation for creating this tool is to
enable the K2 Campaign 9 microlensing science team to use the raw
spacecraft data to identify microlensing events as quickly as possible.

Kadenza can be used both as a command-line tool or using its Python API.

### Example
Convert the raw cadence data into TPF files:
```
$ kadenza-tpf --target 200049143 cadence-data-list.txt pixel-mapping-file.fits
```

Produce a sparse FFI frame for a given channel:
```
$ kadenza-ffi --channel 10 cadence-data-file.fits pixel-mapping-file.fits
``` 

### Installation
If you have a working installation of Python on your system,
you can download and install the latest version as follows:
```
$ git clone https://github.com/barentsen/k2flix.git
$ cd k2flix
$ python setup.py install
```
Kadenza has only been tested under Linux with Python 3.4 at present.
Support for legacy Python will be added soon.

### Using kadenza
TODO: insert a few usage examples here, including the output of `--help`.

### Contributing
To report bugs and request features, please use the [issue tracker](https://github.com/KeplerGO/kadenza/issues) or open a pull request.

### License
Kadenza is made available under the MIT License.
For details see the LICENSE file.
