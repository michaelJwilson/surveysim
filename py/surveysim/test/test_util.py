"""Test surveysim.util.
"""
import unittest
import numpy as np
from astropy.table import Table, Column
from ..util import add_calibration_exposures

class TestUtil(unittest.TestCase):
    """Test surveysim.util.
    """

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_add_calibration_exposures(self):
        """Test adding calibration exposures to science exposures.
        """
        exposures = Table()
        exposures['TILEID'] = Column(np.array([0, 1], dtype=np.int32))
        exposures['PASS'] = Column(np.array([0, 0], dtype=np.int16))
        exposures['RA'] = Column(np.array([0.0, 1.0], dtype=np.float64))
        exposures['DEC'] = Column(np.array([0.0, 1.0], dtype=np.float64))
        exposures['EBMV'] = Column(np.array([0.0, 1.0], dtype=np.float64))
        exposures['NIGHT'] = Column(np.array(['20200101', '20200102'], dtype=(str, 8)))
        exposures['MJD'] = Column(np.array([58849.0, 58850.0], dtype=np.float64))
        exposures['EXPTIME'] = Column(np.array([0.0, 1.0], dtype=np.float64), unit='s')
        exposures['SEEING'] = Column(np.array([0.0, 1.0], dtype=np.float64), unit='arcsec')
        exposures['TRANSPARENCY'] = Column(np.array([0.0, 1.0], dtype=np.float64))
        exposures['AIRMASS'] = Column(np.array([0.0, 1.0], dtype=np.float64))
        exposures['MOONFRAC'] = Column(np.array([0.0, 1.0], dtype=np.float64))
        exposures['MOONALT'] = Column(np.array([0.0, 1.0], dtype=np.float64), unit='deg')
        exposures['MOONSEP'] = Column(np.array([0.0, 1.0], dtype=np.float64), unit='deg')
        exposures['PROGRAM'] = Column(np.array(['DARK', 'DARK'], dtype=(str, 6)))
        exposures['FLAVOR'] = Column(np.array(['science', 'science'], dtype=(str, 7)))
        output = add_calibration_exposures(exposures)
        self.assertEqual(len(output), 14)
        self.assertEqual('EXPID', output.colnames[0])
        self.assertTrue(np.all(output['EXPID'] == np.arange(14, dtype=np.int32)))
        self.assertEqual(output['RA'].unit, 'deg')
        print(output['MJD'])

def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
