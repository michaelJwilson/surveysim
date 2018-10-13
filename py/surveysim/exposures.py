"""Record simulated exposures and collect per-tile statistics.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import desisurvey.tiles


class ExposureList(object):

    def __init__(self, tiles_file=None, max_nexp=50000):
        self.tiles = desisurvey.tiles.get_tiles(tiles_file)
        self._exposures = np.empty(max_nexp, dtype=[
            ('mjd', np.float64),
            ('exptime', np.float32),
            ('tileid', np.int32),
            ('passnum', np.int32),
            ('snr2frac_start', np.float32),
            ('snr2frac_stop', np.float32),
            ('seeing', np.float32),
            ('transp', np.float32),
        ])
        self._tiledata = np.empty(self.tiles.ntiles, dtype=[
            ('exptime', np.float32),
            ('snr2frac', np.float32),
            ('nexp', np.int32)
        ])
        self.reset()

    def reset(self):
        self.nexp = 0
        self._tiledata[:] = 0

    def add(self, mjd, exptime, tileid, passnum, snr2frac_start, snr2frac_stop,
            seeing, transp):
        if self.nexp >= len(self._exposures):
            raise RuntimeError(
                'Need to increase max_nexp={}'.format(len(self._exposures)))
        self._exposures[self.nexp] = (
            mjd, exptime, tileid, passnum, snr2frac_start, snr2frac_stop,
            seeing, transp)
        self.nexp += 1
        tileinfo = self._tiledata[self.tiles.index(tileid)]
        tileinfo['exptime'] += exptime
        tileinfo['snr2frac'] = snr2frac_stop
        tileinfo['nexp'] += 1