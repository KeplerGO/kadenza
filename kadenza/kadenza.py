"""
Converts raw Kepler Cadence Data into astronomer-friendly TPF/FFI FITS files.

Author: Geert Barentsen
"""
import os
import re
import datetime
import argparse

import numpy as np

from astropy import log
from astropy.io import fits
from astropy.utils.console import ProgressBar
from tqdm import tqdm

from . import __version__
from . import calibration


class PixelMappingFile(object):
    """Wraps a Kepler Pixel Mapping Reference file.

    A pixel mapping reference file describes the relationship between the
    pixel values recorded in a Cadence Data File and the pixel positions
    on the Kepler CCDs.  These files tend to have the suffix
        "*_lcm.fits"
    for long cadence, or
        "*_scm.fits"
    for short cadence data.

    Parameters
    ----------
    filename : str
        Path to the Pixel Mapping Reference file (*_lcm.fits or *_scm.fits).
    """
    def __init__(self, filename):
        self.filename = filename
        self.targets = self._targets()

    def _targets(self):
        """Returns a dict mapping target_id => extension/channel number."""
        fts = fits.open(self.filename)
        targets = {}
        for ext in range(1, len(fts)):
            for target_id in np.unique(fts[ext].data['target_id']):
                targets[target_id] = {
                                        "ext": ext,
                                        "channel": fts[ext].header['CHANNEL'],
                                        "module": fts[ext].header['MODULE'],
                                        "output": fts[ext].header['OUTPUT']
                }
        fts.close()
        return targets

    def get_mapping(self, target_id):
        """Returns the pixel mapping information for a specific target.

        Parameters
        ----------
        target_id : str
            The KIC or EPIC identifier of the target.

        Returns
        -------
        idx : `np.ndarray`
            Indexes of the pixel values associated with the target.

        row : `np.ndarray`
            Row numbers for the pixels indexed by pixel_idx.

        column : `np.ndarray`
            Column number for the pixels indexed by pixel_idx.
        """
        fts = fits.open(self.filename)
        ext = self.targets[target_id]['ext']
        idx = np.where(fts[ext].data['target_id'] == target_id)[0]
        row = fts[ext].data['row'][idx]
        column = fts[ext].data['column'][idx]
        fts.close()
        return idx, row, column


class TargetPixelFileFactory(object):
    """
    Parameters
    ----------
    cadence_pixel_files : list of str, or str
        List of paths to the cadence pixel files,
        or path to a text file containing the paths.

    pixel_mapping_file : str
        Path to the pixel mapping file.
    """
    def __init__(self, cadence_pixel_files, pixel_mapping_file, correct_smear=False):
        if type(cadence_pixel_files) is str:
            filenames = [fn.strip() for fn
                         in open(cadence_pixel_files, "r").readlines()]
            self.cadence_pixel_files = filenames
        else:
            self.cadence_pixel_files = cadence_pixel_files
        if correct_smear:
            self.collateral_files = ([fn.replace('-targ.fits',
                                                  '-col.fits')
                                      for fn in self.cadence_pixel_files])
            self.collateral_mapping_fn =  pixel_mapping_file.replace('\?\?\?-\?\?\?-lcm',
                                                                     '000-000-lcc.fits')
        else:
            self.collateral_files = None
        self.pixel_mapping = PixelMappingFile(pixel_mapping_file)
        self.no_cadences = len(self.cadence_pixel_files)

    def get_header_template(self, extension):
        """Returns a template `fits.Header` object for a given extension."""
        template_fn = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   "header-templates",
                                   "tpf-ext{}-header.txt".format(extension))
        return fits.Header.fromtextfile(template_fn)

    def make_tpf(self, target_id):
        tpf = fits.HDUList(self._make_primary_hdu(target_id))
        for ext in self._make_extensions(target_id):
            tpf.append(ext)
        return tpf

    def write_tpf(self, target_id, output_fn=None):
        """Creates and writes a TPF file to disk."""
        if output_fn is None:
            output_fn = 'ktwo{:09d}-unofficial-tpf.fits'.format(target_id)
        log.info("Writing {}".format(output_fn))
        try:
            self.make_tpf(target_id).writeto(output_fn,
                                             overwrite=True,
                                             checksum=True)
        except MemoryError:
            print('MemoryError: {}'.format(output_fn))

    def write_all_tpfs(self):
        """Produce TPF files for all targets in the cadence data."""
        target_ids = list(self.pixel_mapping.targets.keys())
        log.info("Writing {} Target Pixel Files.".format(len(target_ids)))
        ProgressBar.map(self.write_tpf, target_ids, multiprocess=True)

    def _make_primary_hdu(self, target_id):
        """Returns the primary extension of a Target Pixel File."""
        hdu = fits.PrimaryHDU()
        # Copy the default keywords from a template file from the MAST archive
        tmpl = self.get_header_template(0)
        for kw in tmpl:
            hdu.header[kw] = (tmpl[kw], tmpl.comments[kw])
        # Override the defaults where necessary
        hdu.header['ORIGIN'] = "Unofficial data product"
        hdu.header['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d")
        hdu.header['CREATOR'] = "kadenza"
        hdu.header['PROCVER'] = "{}".format(__version__)
        hdu.header['FILEVER'] = "0.0"
        hdu.header['OBJECT'] = target_name(target_id)
        hdu.header['KEPLERID'] = target_id
        hdu.header['CHANNEL'] = self.pixel_mapping.targets[target_id]['channel']
        hdu.header['MODULE'] = self.pixel_mapping.targets[target_id]['module']
        hdu.header['OUTPUT'] = self.pixel_mapping.targets[target_id]['output']
        # Empty a bunch of keywords rather than having incorrect info
        for kw in ["TIMVERSN", "CAMPAIGN", "DATA_REL", "TTABLEID",
                   "RA_OBJ", "DEC_OBJ"]:
            hdu.header[kw] = ""
        return hdu

    def _make_extensions(self, target_id):
        """Create the 'TARGETTABLES' extension."""
        # First: understand where to find and place the pixels
        pixel_idx, row_coords, column_coords = self.pixel_mapping.get_mapping(target_id)
        dcol = max(column_coords) - min(column_coords) + 1
        drow = max(row_coords) - min(row_coords) + 1
        crval1p = min(column_coords)
        crval2p = min(row_coords)

        # Initialize the data tables
        raw_cnts = np.zeros((self.no_cadences, drow, dcol), dtype='int')
        flux = np.zeros((self.no_cadences, drow, dcol), dtype='float32')
        flux_err = np.zeros((self.no_cadences, drow, dcol), dtype='float32')
        flux_bkg = np.zeros((self.no_cadences, drow, dcol), dtype='float32')
        flux_bkg_err = np.zeros((self.no_cadences, drow, dcol), dtype='float32')
        cosmic_rays = np.zeros((self.no_cadences, drow, dcol), dtype='float32')

        raw_cnts[:, :, :] = -1
        flux[:, :, :] = np.nan
        flux_err[:, :, :] = np.nan
        flux_bkg[:, :, :] = np.nan
        flux_bkg_err[:, :, :] = np.nan
        cosmic_rays[:, :, :] = np.nan

        mjd = np.zeros((self.no_cadences), dtype='float64')
        time = np.zeros((self.no_cadences), dtype='float64')
        timecorr = np.zeros((self.no_cadences), dtype='float32')
        cadenceno = np.zeros((self.no_cadences), dtype='int')
        quality = np.zeros((self.no_cadences), dtype='int')
        pos_corr1 = np.zeros((self.no_cadences), dtype='float32')
        pos_corr2 = np.zeros((self.no_cadences), dtype='float32')

        # Open the cadence data files and copy data across
        channel = self.pixel_mapping.targets[target_id]['channel']
        for cad_idx, fn in tqdm(enumerate(self.cadence_pixel_files), desc='Reading cadences', total=len(self.cadence_pixel_files)):
            log.debug("Opening {}".format(fn))
            try:
                cadfile = fits.open(fn)
            except IOError as e:
                log.warning('WARNING: Could not open {}: {}'.format(fn, e))
                continue  # try opening the next cadence file

            # If this is the first file; remember a few keywords we'll need later
            if cad_idx == 0:
                dateobs = cadfile[0].header['DATE-OBS']
                timeobs = cadfile[0].header['TIME-OBS']
                lcfxdoff = cadfile[0].header['LCFXDOFF']
                scfxdoff = cadfile[0].header['SCFXDOFF']
                int_time = cadfile[0].header['INT_TIME']
                crpix1 = cadfile[channel].header['CRPIX1']
                crpix2 = cadfile[channel].header['CRPIX2']
                crval1 = cadfile[channel].header['CRVAL1']
                crval2 = cadfile[channel].header['CRVAL2']
                gain = cadfile[channel].header['GAIN']
                readnois = cadfile[channel].header['READNOIS']
                timslice = cadfile[channel].header['TIMSLICE']
                meanblck = cadfile[channel].header['MEANBLCK']

            # Determine the raw pixel count offsets and number of readouts
            if cadfile[0].header['DATATYPE'].strip() == 'long cadence':
                fixed_offset = cadfile[0].header['LCFXDOFF']
                nreadout = 270
                # Long cadence number
                cadenceno[cad_idx] = cadfile[0].header['LC_INTER']
            else:
                fixed_offset = cadfile[0].header['SCFXDOFF']  # short cadence
                nreadout = 9
                # Short cadence number
                cadenceno[cad_idx] = cadfile[0].header['SC_INTER']

            # Populate cadence time and number
            mjd[cad_idx] = (cadfile[1].header['BSTRTIME'] + cadfile[1].header['BSTPTIME']) / 2.
            time[cad_idx] = mjd[cad_idx] + 2400000.5 - 2454833.0

            # Get smear values
            colldata = CollateralData(self.collateral_mapping_files[cad_idx],
                                      self.collateral_mapping_fn)
            smear_values = colldata.get_smear_at_columns(column_coords, channel)

            # Determine pixel values
            pixelvalues_raw = cadfile[channel].data['orig_value'][:]
            pixelvalues_adu = calibration.raw_counts_to_adu(pixelvalues_raw,
                                                            fixed_offset, meanblck, nreadout)
            pixelvalues_adu -= smear_values
            # Rough calibration: uses mean black instead of observed black!
            exposure_time = int_time * nreadout
            pixelvalues_flux = (pixelvalues_adu - nreadout*meanblck) * gain / exposure_time

            # Populate pixel arrays
            for idx in range(len(row_coords)):

                i, j = ccd2mask(1, 1, crval1p, crval2p,
                                1, 1, column_coords[idx], row_coords[idx])
                raw_cnts[cad_idx, j, i] = pixelvalues_raw[pixel_idx[idx]]
                flux[cad_idx, j, i] = pixelvalues_flux[pixel_idx[idx]]
                flux_err[cad_idx, j, i] = np.sqrt(pixelvalues_flux[pixel_idx[idx]])

            cadfile.close()
            del cadfile

        # Turn the data arrays into fits columns and initialize the HDU
        coldim = '({},{})'.format(dcol, drow)
        eformat = '{}E'.format(drow * dcol)
        jformat = '{}J'.format(drow * dcol)
        cols = []
        cols.append(fits.Column(name='TIME', format='D', unit='BJD - 2454833',
                                array=time))
        cols.append(fits.Column(name='TIMECORR', format='E', unit='D',
                                array=timecorr))
        cols.append(fits.Column(name='CADENCENO', format='J', array=cadenceno))
        cols.append(fits.Column(name='RAW_CNTS', format=jformat, unit='count',
                                dim=coldim, array=raw_cnts))
        cols.append(fits.Column(name='FLUX', format=eformat, unit='e-/s',
                                dim=coldim, array=flux))
        cols.append(fits.Column(name='FLUX_ERR', format=eformat, unit='e-/s',
                                dim=coldim, array=flux_err))
        cols.append(fits.Column(name='FLUX_BKG', format=eformat, unit='e-/s',
                                dim=coldim, array=flux_bkg))
        cols.append(fits.Column(name='FLUX_BKG_ERR', format=eformat, unit='e-/s',
                                dim=coldim, array=flux_bkg_err))
        cols.append(fits.Column(name='COSMIC_RAYS', format=eformat, unit='e-/s',
                                dim=coldim, array=cosmic_rays))
        cols.append(fits.Column(name='QUALITY', format='J', array=quality))
        cols.append(fits.Column(name='POS_CORR1', format='E', unit='pixels',
                                array=pos_corr1))
        cols.append(fits.Column(name='POS_CORR2', format='E', unit='pixels',
                                array=pos_corr2))
        coldefs = fits.ColDefs(cols)
        hdu = fits.BinTableHDU.from_columns(coldefs)

        # Set the header with defaults
        tmpl = self.get_header_template(1)
        for i, kw in enumerate(tmpl):
            if kw in ['XTENSION', 'KEPLERID', 'NAXIS1', 'NAXIS2']:
                continue
            hdu.header[kw] = (tmpl[kw], tmpl.comments[kw])
        # Override the defaults where necessary
        for n in [5, 6, 7, 8, 9]:
            hdu.header["TFORM{}".format(n)] = eformat
            hdu.header["TDIM{}".format(n)] = coldim
        hdu.header['TFORM4'] = jformat
        hdu.header['TDIM4'] = coldim

        # Set a WCS. TODO: check if this makes sense!
        hdu.header['1CRV4P'] = crval1p
        hdu.header['2CRV4P'] = crval2p
        hdu.header['1CRPX4'] = crpix1
        hdu.header['2CRPX4'] = crpix2
        hdu.header['1CRVL4'] = crval1
        hdu.header['2CRVL4'] = crval2

        hdu.header['1CRV5P'] = crval1p
        hdu.header['2CRV5P'] = crval2p
        hdu.header['1CRPX5'] = crpix1
        hdu.header['2CRPX5'] = crpix2
        hdu.header['1CRVL5'] = crval1
        hdu.header['2CRVL5'] = crval2

        hdu.header['1CRV6P'] = crval1p
        hdu.header['2CRV6P'] = crval2p

        hdu.header['1CRPX6'] = crpix1
        hdu.header['2CRPX6'] = crpix2
        hdu.header['1CRVL6'] = crval1
        hdu.header['2CRVL6'] = crval2

        hdu.header['1CRV7P'] = crval1p
        hdu.header['2CRV7P'] = crval2p

        hdu.header['1CRPX7'] = crpix1
        hdu.header['2CRPX7'] = crpix2
        hdu.header['1CRVL7'] = crval1
        hdu.header['2CRVL7'] = crval2

        hdu.header['1CRV8P'] = crval1p
        hdu.header['2CRV8P'] = crval2p

        hdu.header['1CRPX8'] = crpix1
        hdu.header['2CRPX8'] = crpix2
        hdu.header['1CRVL8'] = crval1
        hdu.header['2CRVL8'] = crval2

        hdu.header['1CRV9P'] = crval1p
        hdu.header['2CRV9P'] = crval2p

        hdu.header['1CRPX9'] = crpix1
        hdu.header['2CRPX9'] = crpix2
        hdu.header['1CRVL9'] = crval1
        hdu.header['2CRVL9'] = crval2

        hdu.header['OBJECT'] = target_name(target_id)
        hdu.header['KEPLERID'] = target_id

        hdu.header['EXPOSURE'] = (time[-1] - time[0]) * 0.92063492
        hdu.header['TELAPSE'] = time[-1] - time[0]
        hdu.header['LIVETIME'] = (time[-1] - time[0]) * 0.92063492
        hdu.header['TSTART'] = time[0]
        hdu.header['TSTOP'] = time[-1]
        hdu.header['LC_START'] = mjd[0]
        hdu.header['LC_END'] = mjd[-1]

        hdu.header['DATE-OBS'] = dateobs + ':' + timeobs + 'Z'
        hdu.header['DATE-END'] = ''  # TODO?

        hdu.header['GAIN'] = gain
        hdu.header['READNOIS'] = readnois
        hdu.header['TIMSLICE'] = timslice
        hdu.header['MEANBLCK'] = meanblck
        hdu.header['LCFXDOFF'] = lcfxdoff
        hdu.header['SCFXDOFF'] = scfxdoff

        hdu.header['DBTRESH'] = ''
        hdu.header['BLKALGO'] = ''

        # Now make the apperture mask extension
        mask = np.ones((drow, dcol), dtype='int32')
        aper_hdu = fits.ImageHDU(mask)

        # Set the header from the template TPF again
        tmpl = self.get_header_template(2)
        for i, kw in enumerate(tmpl):
            aper_hdu.header[kw] = (tmpl[kw], tmpl.comments[kw])

        aper_hdu.header['OBJECT'] = target_name(target_id)
        aper_hdu.header['KEPLERID'] = target_id

        aper_hdu.header['CRPIX1'] = crpix1
        aper_hdu.header['CRPIX2'] = crpix2
        aper_hdu.header['CRVAL1'] = crval1
        aper_hdu.header['CRVAL2'] = crval2

        aper_hdu.header['CRVAL1P'] = crval1p
        aper_hdu.header['CRVAL2P'] = crval2p

        return (hdu, aper_hdu)


class FullFrameImageFactory(object):
    """
    Parameters
    ----------
    cadence_pixel_file : str
        Path to the cadence data file.

    pixel_mapping_file : str
        Path to the pixel mapping file.
    """
    def __init__(self, cadence_pixel_file, pixel_mapping_file):
        self.cadence_pixel_file = cadence_pixel_file
        self.pixel_mapping_file = pixel_mapping_file

    def get_header_template(self, extension):
        """Returns a template `fits.Header` object for a given extension."""
        if extension > 1:
            extension = 1
        template_fn = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                   "header-templates",
                                   "ffi-ext{}-header.txt".format(extension))
        return fits.Header.fromtextfile(template_fn)

    def make_ffi(self):
        return self._make_hdulist()

    def write_ffi(self, output_fn=None):
        """Creates and writes a TPF file to disk."""
        basename = os.path.basename(self.cadence_pixel_file)
        lastcad = re.sub('_lcs-targ.fits', '', basename).strip('kplr')
        if output_fn is None:
            basename = os.path.basename(self.cadence_pixel_file)
            if "lcs-targ" in basename:
                output_fn = basename.replace("lcs-targ", "kadenza_ffi_raw")
            else:
                output_fn = "sparse_cadence_ffi_raw.fits"
        log.info("Writing {}".format(output_fn))
        self.make_ffi().writeto(output_fn, overwrite=True, checksum=True)

    def _make_hdulist(self):
        hdu0 = fits.PrimaryHDU()

        # Set the header with defaults
        tmpl = self.get_header_template(0)
        for i, kw in enumerate(tmpl):
            hdu0.header[kw] = (tmpl[kw], tmpl.comments[kw])
        # Override the primary extension defaults
        hdu0.header['ORIGIN'] = "Unofficial data product"
        hdu0.header['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d")
        hdu0.header['CREATOR'] = "kadenza"
        hdu0.header['PROCVER'] = "{}".format(__version__)

        hdulist = fits.HDUList(hdu0)

        # open the Cadence Pixel File and Pixel Mapping File
        incpf = fits.open(self.cadence_pixel_file, memmap=True)
        inpmrf = fits.open(self.pixel_mapping_file, memmap=True)

        # grab the dates and times from the header
        card = incpf[0].header
        startkey = card['STARTIME']
        endkey = card['END_TIME']
        dateobs = card['DATE-OBS']
        timeobs = card['TIME-OBS']
        tstart = startkey + 2400000.5
        tstop = endkey + 2400000.5
        deadc = 0.92063492
        telapse = tstop - tstart
        exposure = telapse * deadc

        for ext in ProgressBar(range(1, 85)):
            card = incpf[ext].header
            gain = card['GAIN']
            readnois = card['READNOIS']
            timslice = card['TIMSLICE']
            meanblck = card['MEANBLCK']

            card = inpmrf[ext].header
            chankey = card['CHANNEL']
            modkey = card['MODULE']
            outkey = card['OUTPUT']
            r = inpmrf[ext].data.field('row')[:]
            c = inpmrf[ext].data.field('column')[:]
            k = inpmrf[ext].data.field('target_id')[:]
            f = incpf[ext].data.field('orig_value')[:]

            # initialize FFI image
            ffimage = np.zeros((1070, 1132), dtype="int32")
            ffimage[:, :] = -1

            # populate FFI image
            ffimage[r[:len(k)], c[:len(k)]] = f[:len(k)]

            # create image extension headers
            hdu = fits.ImageHDU(ffimage)

            # Set the header with defaults
            tmpl = self.get_header_template(ext)
            for i, kw in enumerate(tmpl):
                hdu.header[kw] = (tmpl[kw], tmpl.comments[kw])
            # Overrides:
            hdu.header['EXTNAME'] = 'MOD.OUT %d.%d' % (modkey, outkey)
            hdu.header['CHANNEL'] = chankey
            hdu.header['MODULE'] = modkey
            hdu.header['OUTPUT'] = outkey
            hdu.header['MJDSTART'] = startkey
            hdu.header['MJDEND'] = endkey
            hdu.header['BJDREFI'] = 2454833
            hdu.header['BJDREFF'] = 0.0
            hdu.header['TSTART'] = tstart
            hdu.header['TSTOP'] = tstop
            hdu.header['TELAPSE'] = telapse
            hdu.header['EXPOSURE'] = exposure
            hdu.header['LIVETIME'] = exposure
            hdu.header['TIMEDEL'] = telapse
            hdu.header['DATE-OBS'] = ''
            hdu.header['DATE-END'] = '%sT%sZ' % (dateobs, timeobs)
            hdu.header['BUNIT'] = "count"
            hdu.header['BARYCORR'] = -1.0
            hdu.header['BACKAPP'] = False
            hdu.header['DEADAPP'] = True
            hdu.header['VIGNAPP'] = False
            hdu.header['GAIN'] = gain
            hdu.header['READNOIS'] = readnois
            hdu.header['TIMSLICE'] = timslice
            hdu.header['MEANBLCK'] = meanblck
            hdu.header['NPIXSAP'] = ""

            hdulist.append(hdu)

        incpf.close()
        inpmrf.close()
        return hdulist


""" Helper functions """


def target_name(target_id):
    """Returns the target name, 'KIC {target_id}' or 'EPIC {target_id}'."""
    if int(target_id) < 2e8:
        return "KIC {}".format(target_id)
    return "EPIC {}".format(target_id)


def ccd2mask(crpix1, crpix2, crval1, crval2,
             cdelt1, cdelt2, ccd_column, ccd_row):
    mask_column = (ccd_column - crval1) * cdelt1 + crpix1 - 1
    mask_row = (ccd_row - crval2) * cdelt2 + crpix2 - 1
    return mask_column, mask_row


""" Command-line interface """


def kadenza_tpf_main(args=None):
    parser = argparse.ArgumentParser(
                description="Turn raw Kepler Cadence Data into "
                            "uncalibrated Target Pixel Files (TPF).")
    parser.add_argument("-t", "--target", metavar='target_id',
                        nargs="?", type=int,
                        help="only produce a TPF file "
                             "for a specific EPIC/KIC target_id")
    parser.add_argument('cadencefile_list', nargs=1,
                        help="Path to a text file that lists the cadence data "
                             "files to use (one file per line). "
                             "These files are named '*_lcs-targ.fits' for "
                             "long cadence and '*_scs-targ.fits' for short "
                             "cadence.")
    parser.add_argument('pixelmap_file', nargs=1,
                        help="Path to the pixel mapping reference file. "
                             "This file is named '*_lcm.fits' for long "
                             "cadence and '*_scm.fits' for short cadence.")
    args = parser.parse_args(args)

    # Allow cadence file to be given rather than a list
    if args.cadencefile_list[0].endswith("lcs-targ.fits"):
        cflist = args.cadencefile_list
    else:
        cflist = args.cadencefile_list[0]
    factory = TargetPixelFileFactory(cflist,
                                     args.pixelmap_file[0])
    if args.target is None:
        factory.write_all_tpfs()
    else:
        factory.write_tpf(args.target)


def kadenza_ffi_main(args=None):
    parser = argparse.ArgumentParser(
                description="Converts a raw Kepler Cadence Data file into "
                            "an uncalibrated, sparse Full Frame Image (FFI). "
                            "The output file contains 84 image extensions "
                            "which correspond to the different Kepler CCDs. "
                            "The units are counts and all unobserved pixels "
                            "are set to -1.")
    parser.add_argument('cadence_file', nargs=1,
                        help="path to the '*_lcs-targ.fits' cadence data file")
    parser.add_argument('pixelmap_file', nargs=1,
                        help="path to the '*_lcm.fits' "
                             "pixel mapping reference file")
    parser.add_argument('correct_smear', nargs=1,
                        default=False)
    args = parser.parse_args(args)

    factory = FullFrameImageFactory(args.cadence_file[0],
                                    args.pixelmap_file[0])
    factory.write_ffi()


if __name__ == "__main__":
    # Example use
    cadence_files = "../sandbox/filenames.txt"
    pixel_mapping = "../sandbox/pixel-mappings/kplr2009115065205-013-013_lcm.fits"
    factory = TargetPixelFileFactory(cadence_files, pixel_mapping)
    # Get an arbitrary target_id
    target_id = list(factory.pixel_mapping.targets.keys())[5]
    factory.write_tpf(target_id)
    #factory.write_all_tpfs()
