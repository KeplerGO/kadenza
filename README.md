# Kadenza: Kepler/K2 Cadence Data Reader 
***Converts raw cadence data from the Kepler spacecraft into astronomer-friendly FITS files.***

[![Travis status](https://travis-ci.org/KeplerGO/kadenza.svg)](https://travis-ci.org/KeplerGO/kadenza) [![MIT license](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/barentsen/k2flix/blob/master/LICENSE) [![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.344973.svg)](https://doi.org/10.5281/zenodo.344973) [![ADS Bibcode](https://img.shields.io/badge/NASA%20ADS-2018ascl.soft03005B-blue.svg)](http://adsabs.harvard.edu/abs/2018ascl.soft03005B)

Kadenza enables time-critical data analyses to be carried out using NASA's Kepler Space Telescope
by enabling users to convert the "level 0" raw data files into user-friendly Target Pixel Files.

The primary motivation for this tool was to enable the microlensing, supernova, and exoplanet communities
to create quicklook lightcurves for transient events and exoplanets which require rapid follow-up.

Kadenza can be used both as a command-line tool or using its Python API.
Most users of Kepler/K2 data will not require this tool,
because it does not produce fully-calibrated data products.


### Installation
If you have a working installation of Python on your system,
you can download and install the latest version as follows:
```
$ git clone https://github.com/KeplerGO/kadenza.git
$ cd kadenza
$ python setup.py install
```
Kadenza has only been tested under Linux with Python 3 at present.

### Usage

#### Example 1: Creating a (sparse) Full Frame Image

For K2 Campaign 9, the [data archive at MAST](https://archive.stsci.edu/pub/k2/)
provides access to the raw *cadence data files*.
There is one such file per long cadence,
named `kplrYYYYDDDHHMMSS_lcs-targ.fits`,
which provides the pixel counts in that cadence
for all those pixels which were pre-selected to be downlinked
from the spacecraft (roughly 3% of the 95-megapixel camera).
The timestamp that is recorded in the filename
refers to the end of the long cadence,
which is a period of 0.4904h during which the spacecraft
summed 270 exposures of 6.02s each.

To reconstruct two-dimensional images from the cadence data files,
we need an additional file called the *pixel mapping reference file*.
This file specifies the (column, row) CCD coordinates for each value
in the one-dimensional cadence data arrays.
This information is kept in a separate file to avoid having to repeat
the pixel coordinates in each cadence file.
For K2 Campaign 9a the mapping file you need is called
`kplr2016068153039-085-085_lcm.fits` and can be obtained from the archive. 

One you have obtained the raw data, you can then use Kadenza
to convert them into two-dimensional images.
Once intalled (see below), Kadenza adds the `kadenza-ffi` tool
to your command-line which you simply call as follows:
```
$ kadenza-ffi cadence-data-file pixel-mapping-file
``` 

In our example, we'd execute:
```
$ kadenza-ffi kplrYYYYDDDHHMMSS_lcs-targ.fits kplr2016068153039-085-085_lcm.fits
```

This will create a new file called `kplrYYYYDDDHHMMSS_kadenza_ffi_raw.fits`
which is a FITS file with 84 image extensions, 
each corresponding to the different Kepler CCD channels.
The units are counts and unobserved pixels are set to "-1".
Timestamps refer to the end of each cadence can be obtained from the
filename or the DATE-END keyword in the header.


#### Example 2: Creating a Target Pixel File

You can also convert *all* the raw cadence data for a given campaign and target into a Target Pixel File (TPF) file using the `kadenza-tpf` tool.  For example, to create a Target Pixel File for the target with EPIC ID 247520207 observed in Campaign 13, you would first download all the cadence data, e.g. using `wget`:
```
$ wget -t1 -e robots=off --no-parent --continue --recursive https://archive.stsci.edu/pub/k2/raw_cadence_data/c13/
```
Next, create a text file which lists the paths of all the cadence files you downloaded sorted by filename, e.g. using `find`:
```
$ find /location/where/you/downloaded/the/raw/cadence/data -name '*lcs-targ.fits' | sort > raw-long-cadence-files.txt
```

Next, download the long cadence pixel mapping reference file from MAST:
```
$ wget https://archive.stsci.edu/pub/k2/pmrfs/c13/kplr2017032193929-091-091_lcm.fits
```

And finally, run `kadenza-tpf` to create an *uncalibrated* target pixel file:
```
$ kadenza-tpf --target 247520207 raw-long-cadence-files.txt kplr2017032193929-091-091_lcm.fits
```

This will create an unofficial Target Pixel File called `ktwo247520207-unofficial-tpf.fits`.

### Caveats

There are two main caveats to be aware of:
 * The present version does not set all header keywords exactly as they would
   appear in an official product.  In particular, the WCS keywords are
   untested.
 * This tool does *not* calibrate the data, it merely serves to
   transform the raw pixel counts into a FITS format that is similar to,
   but not at all identical, to the official pipeline products.
   
### Commands provided
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

### Contributing
To report bugs and request features, please use the [issue tracker](https://github.com/KeplerGO/kadenza/issues) or open a pull request.

### License
Kadenza is made available under the MIT License.
For details see the LICENSE file.

### Citation

`kadenza` was created by [Geert Barentsen](https://github.com/barentsen)
and [José Vinícius de Miranda Cardoso](https://github.com/mirca) for NASA's Kepler/K2 Guest Observer Office.
If this tool aided your research, please consider offering co-authorship to the authors,
or at the very least, cite this tool using both the [DOI identifier](https://doi.org/10.5281/zenodo.344973)
and the [ASCL entry](https://ascl.net/code/v/1891) using the following BibTeX entry:

```
@MISC{2018ascl.soft03005B,
  author        = {{Barentsen}, G. and {Cardoso}, J.~V.~d.~M.},
  title         = "{Kadenza: Kepler/K2 Raw Cadence Data Reader}",
  keywords      = {Software },
  howpublished  = {Astrophysics Source Code Library},
  year          = 2018,
  month         = mar,
  archivePrefix = "ascl",
  eprint        = {1803.005},
  adsurl        = {http://adsabs.harvard.edu/abs/2018ascl.soft03005B},
  adsnote       = {Provided by the SAO/NASA Astrophysics Data System},
  doi           = {10.5281/zenodo.344973},
  url           = {https://doi.org/10.5281/zenodo.344973}
}
```
