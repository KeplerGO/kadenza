import os

from .. import PACKAGEDIR
from .. import calibration

EXAMPLE_COLLATERAL = os.path.join(PACKAGEDIR, 'tests', 'data', 'kplr2016130002211_lcs-col.fits')
EXAMPLE_COLLATERAL_MAPPING = os.path.join(PACKAGEDIR, 'tests', 'data', 'kplr2016068153039-000-000_lcc.fits')


def test_calibration():
    col = calibration.CollateralData(EXAMPLE_COLLATERAL, EXAMPLE_COLLATERAL_MAPPING)
    frame = col.smear_frame(channel=31)
