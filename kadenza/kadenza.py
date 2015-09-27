"""
Converts raw Kepler Cadence Data into astronomer-friendly TPF/FFI FITS files.

Author: Geert Barentsen
"""
import os
import re
import datetime

import numpy as np

from astropy import log
from astropy.io import fits
from astropy.utils.console import ProgressBar


class PixelMappingFile():
    """Wraps a Kepler Pixel Mapping Reference file.

    A pixel mapping reference file describes the relationship between th
    pixel values recorded in a Cadence Data File and the pixel positions
    on the Kepler CCDs.  These files tend to have the suffix
        "*_lcm.fits"
    for long cadence and
        "*_scm.fits"
    for short cadence data.

    Parameters
    ----------
    filename : str
        Path to the Pixel Mapping Reference file (*_lcm.fits or *_scm.fits).
    """
    def __init__(self, filename):
        self.filename = filename
        self.hdulist = fits.open(filename)
        self.targets = self._targets()

    def _targets(self):
        """Returns a dict mapping target_id => extension/channel number."""
        targets = {}
        for ext in range(1, len(self.hdulist)):
            target_ids = np.unique(self.hdulist[ext].data['target_id'])
            targets.update(zip(target_ids, [ext]*len(target_ids)))
        return targets

    def get_extension(self, target_id):
        """Returns the extension/channel number for a target."""
        return self.targets[target_id]

    def get_channel(self, target_id):
        # Note that channel == extension number ... or so we believe!
        return self.hdulist[self.get_extension(target_id)].header['CHANNEL']

    def get_module(self, target_id):
        return self.hdulist[self.get_extension(target_id)].header['MODULE']

    def get_output(self, target_id):
        return self.hdulist[self.get_extension(target_id)].header['OUTPUT']

    def get_mapping(self, target_id):
        """Returns the pixel mapping information for a specific target.

        Parameters
        ----------
        target_id : str
            The KIC or EPIC identifier of the target.

        Returns
        -------
        idx : `numpy.ndarray`
            Indexes of the pixel values associated with the target.

        row : `numpy.ndarray`
            Row numbers for the pixels indexed by pixel_idx.

        column : `numpy.ndarray`
            Column number for  the pixels indexed by pixel_idx.
        """
        ext = self.get_extension(target_id)
        idx = np.where(self.hdulist[ext].data['target_id'] == target_id)[0]
        row = self.hdulist[ext].data['row'][idx]
        column = self.hdulist[ext].data['column'][idx]
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
    def __init__(self, cadence_pixel_files, pixel_mapping_file):
        if type(cadence_pixel_files) is str:
            filenames = [fn.strip() for fn
                         in open(cadence_pixel_files, "r").readlines()]
            self.cadence_pixel_files = filenames
        else:
            self.cadence_pixel_files = cadence_pixel_files
        self.pixel_mapping = PixelMappingFile(pixel_mapping_file)

        self.no_cadences = len(self.cadence_pixel_files)
        self.template = fits.open("tpf-template.fits.gz")

    def make_tpf(self, target_id):
        tpf = fits.HDUList(self._make_primary_hdu(target_id))
        for ext in self._make_extensions(target_id):
            tpf.append(ext)
        return tpf

    def make_all_tpfs(self):
        pass

    def save_tpf(self, target_id, output_fn=None):
        """Creates and writes a TPF file to disk."""
        basename = os.path.basename(self.cadence_pixel_files[-1])
        lastcad = re.sub('_lcs-targ.fits', '', basename).strip('kplr')
        if output_fn is None:
            output_fn = 'mod{}.{}_kplr{:09d}-{}_lpd-targ.fits'.format(
                                self.pixel_mapping.get_module(target_id),
                                self.pixel_mapping.get_module(target_id),
                                target_id,
                                lastcad)
        log.info("Writing {}".format(output_fn))
        self.make_tpf(target_id).writeto(output_fn,
                                         clobber=True,
                                         checksum=True)

    def _make_primary_hdu(self, target_id):
        """Returns the primary extension of a Target Pixel File."""
        hdu = fits.PrimaryHDU()
        # Copy the default keywords from a template file from the MAST archive
        tmpl = self.template[0].header
        for kw in tmpl:
            if kw in ['CHECKSUM']:  # TODO: can this be removed?
                continue
            hdu.header[kw] = (tmpl[kw], tmpl.comments[kw])
        # Override the defaults where necessary
        hdu.header['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d")
        hdu.header['CREATOR'] = "kadenza"
        hdu.header['PROCVER'] = "1.0"
        hdu.header['FILEVER'] = "0.0"
        hdu.header['TIMVERSN'] = ""
        hdu.header['OBJECT'] = target_name(target_id)
        hdu.header['KEPLERID'] = target_id
        hdu.header['CHANNEL'] = self.pixel_mapping.get_channel(target_id)
        hdu.header['MODULE'] = self.pixel_mapping.get_module(target_id)
        hdu.header['OUTPUT'] = self.pixel_mapping.get_output(target_id)
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
        channel = self.pixel_mapping.get_channel(target_id)
        for cad_idx, fn in enumerate(self.cadence_pixel_files):
            log.debug("Opening {}".format(fn))
            cadfile = fits.open(fn)
            mjd[cad_idx] = cadfile[0].header['MID_TIME']
            time[cad_idx] = mjd[cad_idx] + 2400000.5 - 2454833.0
            cadenceno[cad_idx] = cadfile[0].header['LC_COUNT']

            raw = cadfile[channel].data['orig_value'][:]

            for idx in range(len(row_coords)):
                i, j = ccd2mask(1, 1, crval1p, crval2p,
                                1, 1, column_coords[idx], row_coords[idx])
                raw_cnts[cad_idx, j, i] = raw[pixel_idx[idx]]
                flux[cad_idx, j, i] = raw[pixel_idx[idx]]
                flux_err[cad_idx, j, i] = np.sqrt(raw[pixel_idx[idx]])

            # Remember a few keywords we'll need for the header
            if cad_idx == 0:
                dateobs = cadfile[0].header['DATE-OBS']
                timeobs = cadfile[0].header['TIME-OBS']
                lcfxdoff = cadfile[0].header['LCFXDOFF']
                scfxdoff = cadfile[0].header['SCFXDOFF']
                crpix1 = cadfile[channel].header['CRPIX1']
                crpix2 = cadfile[channel].header['CRPIX2']
                crval1 = cadfile[channel].header['CRVAL1']
                crval2 = cadfile[channel].header['CRVAL2']
                gain = cadfile[channel].header['GAIN']
                readnois = cadfile[channel].header['READNOIS']
                timslice = cadfile[channel].header['TIMSLICE']
                meanblck = cadfile[channel].header['MEANBLCK']

            cadfile.close()

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
        cols.append(fits.Column(name='RAW_CNTS', format=jformat, unit='counts',
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
        tmpl = self.template[1].header
        for i, kw in enumerate(tmpl):
            if kw in ['XTENSION', 'CHECKSUM', 'KEPLERID', 'NAXIS1']:
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
        for i, kw in enumerate(tmpl):
            if kw in ['CHECKSUM']:
                continue
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


"""Helper functions"""


def target_name(target_id):
    """Returns the target name, 'KIC {target_id}' or 'EPIC {target_id}'."""
    if int(target_id) < 2e9:
        return "KIC {}".format(target_id)
    return "EPIC {}".format(target_id)


def ccd2mask(crpix1, crpix2, crval1, crval2,
             cdelt1, cdelt2, ccd_column, ccd_row):
    mask_column = (ccd_column - crval1) * cdelt1 + crpix1 - 1
    mask_row = (ccd_row - crval2) * cdelt2 + crpix2 - 1
    return mask_column, mask_row


if __name__ == "__main__":
    # Example use
    cadence_files = "filelist-all.txt"
    pixel_mapping = "../sandbox/pixel-mappings/kplr2009115065205-013-013_lcm.fits"
    factory = TargetPixelFileFactory(cadence_files, pixel_mapping)
    # Get an arbitrary target_id
    target_id = list(factory.pixel_mapping.targets.keys())[5]
    factory.save_tpf(target_id)
