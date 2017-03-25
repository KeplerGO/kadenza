"""Implements tools for reading collateral data and calibrating fluxes.

TODO:
Implement black correction
Correct smear for black
"""
import numpy as np
from astropy.io import fits


# Collateral pixel types
# Source: https://archive.stsci.edu/k2/manuals/KADN-26315.pdf
PIXELTYPE = {'BlackLevel': 1,
             'Masked': 2,
             'VirtualSmear': 3,
             'BlackMasked': 4,
             'BlackVirtual': 5}


def raw_counts_to_adu(raw_counts, fixed_offset, meanblck, nreadout):
    """Convert Kepler raw pixel counts photometer ADU.

    Converts raw cadence pixel counts to Analog to Digital Units read off
    the photometer following the Kepler Archive Manual (p23).

    This is necessary because the spacecraft software subtracts
    a mean black value from each channel to have the pixel values
    of all channels center around zero, so that a single compression
    table can be used to compress all pixel values.
    However to avoid negative values, a large positive offset value
    ('LCFXDOFF' or 'SCFXDOFF') is then added that lifts the pixel
    values to a level that matches the range in the compression
    table designed for the pixel counts.
    """
    return raw_counts + nreadout*meanblck - fixed_offset


class CollateralData(object):
    """Read collateral pixels from raw cadence files.

    Parameters
    ----------
    collateral_file : str
        Path to the collateral data file (*-col.fits).

    collateral_mapping_file : str
        Path to the collateral pixel mapping file (*-lcc.fits).
    """
    def __init__(self, collateral_fn, mapping_fn):
        self.collateral_fn = collateral_fn
        self.mapping_fn = mapping_fn
        self.collateral = fits.open(collateral_fn)
        self.mapping = fits.open(mapping_fn)
        if self.collateral[0].header['DATATYPE'].strip() == 'long cadence':
            self.fixed_offset = self.collateral[0].header['LCFXDOFF']
            self.nreadout = 270
        else:
            self.fixed_offset = self.collateral[0].header['SCFXDOFF']
            self.nreadout = 9

    def black_frame(self, channel):
        frame = np.zeros((1070, 1132), dtype="int32")
        mask = self.mapping[channel].data['col_pixel_type'] == PIXELTYPE['BlackLevel']
        black_values = raw_counts_to_adu(self.collateral[channel].data['orig_value'][mask],
                                         fixed_offset=self.fixed_offset,
                                         meanblck=self.collateral[channel].header['MEANBLCK'],
                                         nreadout=self.nreadout)
        # The smear values of multiple rows are binned on board
        black_values = black_values / self.collateral[0].header['NCOLBLCK']
        rows = self.mapping[channel].data['pixel_offset'][mask]
        for idx in range(len(black_values)):
            frame[rows[idx], :] = black_values[idx]
        return frame

    def smear_frame(self, channel):
        frame = np.zeros((1070, 1132), dtype="int32")
        mask = self.mapping[channel].data['col_pixel_type'] == PIXELTYPE['Masked']
        smear_values = raw_counts_to_adu(self.collateral[channel].data['orig_value'][mask],
                                         fixed_offset=self.fixed_offset,
                                         meanblck=self.collateral[channel].header['MEANBLCK'],
                                         nreadout=self.nreadout)
        # The smear values of multiple rows are binned on board
        smear_values = smear_values / self.collateral[0].header['NROWMASK']
        columns = self.mapping[channel].data['pixel_offset'][mask]
        for idx in range(len(smear_values)):
            frame[:, columns[idx]] = smear_values[idx]
        return frame

    def get_smear_at_columns(self, columns, channel):
        """
        columns: list-like of integers
        """
        columns = list(columns)
        mask = self.mapping[channel].data['col_pixel_type'] == PIXELTYPE['Masked']
        smear_at_columns = raw_counts_to_adu(self.collateral[channel].data['orig_value'][mask][np.array(columns) - 12],
                                            fixed_offset=self.fixed_offset,
                                            meanblck=self.collateral[channel].header['MEANBLCK'],
                                            nreadout=self.nreadout)
        smear_at_columns = smear_at_columns / self.collateral[0].header['NROWMASK']
        return smear_at_columns

    def save_smear(self, channel, output_fn='smear.fits'):
        smear = self.smear_frame(channel)
        print('Writing {}'.format(output_fn))
        fits.ImageHDU(smear).writeto(output_fn)
        from matplotlib.image import imsave
        print('Writing {}.png'.format(output_fn))
        imsave(output_fn + '.png', smear, cmap='gray')

    def save_black(self, channel, output_fn='black.fits'):
        black = self.black_frame(channel)
        print('Writing {}'.format(output_fn))
        fits.ImageHDU(black).writeto(output_fn)
        from matplotlib.image import imsave
        print('Writing {}.png'.format(output_fn))
        imsave(output_fn + '.png', black, cmap='gray')


if __name__ == '__main__':
    # Demo code:
    EXAMPLE_COLLATERAL = 'tests/data/kplr2016130002211_lcs-col.fits'
    EXAMPLE_COLLATERAL_MAPPING = 'tests/data/kplr2016068153039-000-000_lcc.fits'
    CHANNEL = 31
    col = CollateralData(EXAMPLE_COLLATERAL, EXAMPLE_COLLATERAL_MAPPING)
    col.save_smear(channel=CHANNEL, output_fn='smear-{}.fits'.format(CHANNEL))
    col.save_black(channel=CHANNEL, output_fn='black-{}.fits'.format(CHANNEL))
