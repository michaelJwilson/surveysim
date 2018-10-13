"""Record simulated exposures and collect per-tile statistics.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import astropy.io.fits

import desiutil.log

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
            ('sky', np.float32),
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
            seeing, transp, sky):
        if self.nexp >= len(self._exposures):
            raise RuntimeError(
                'Need to increase max_nexp={}'.format(len(self._exposures)))
        self._exposures[self.nexp] = (
            mjd, exptime, tileid, passnum, snr2frac_start, snr2frac_stop,
            seeing, transp, sky)
        self.nexp += 1
        tileinfo = self._tiledata[self.tiles.index(tileid)]
        tileinfo['exptime'] += exptime
        tileinfo['snr2frac'] = snr2frac_stop
        tileinfo['nexp'] += 1

    def save(self, name='exposures.fits', comment='', overwrite=True):
        hdus = astropy.io.fits.HDUList()
        header = astropy.io.fits.Header()
        header['TILES'] = self.tiles.tiles_file
        header['NEXP'] = self.nexp
        header['COMMENT'] = comment
        hdus.append(astropy.io.fits.PrimaryHDU(header=header))
        hdus.append(astropy.io.fits.BinTableHDU(self._exposures[:self.nexp], name='EXPOSURES'))
        hdus.append(astropy.io.fits.BinTableHDU(self._tiledata, name='TILEDATA'))
        config = desisurvey.config.Configuration()
        name = config.get_path(name)
        hdus.writeto(name, overwrite=overwrite)
        log = desiutil.log.get_logger()
        log.info('Saved {} exposures to {}'.format(self.nexp, name))
        if comment:
            log.info('Saved with comment "{}".'.format(header['COMMENT']))


def load(name, extra_nexp=0):
    config = desisurvey.config.Configuration()
    name = config.get_path(name)
    with astropy.io.fits.open(name) as hdus:
        header = hdus[0].header
        comment = header['COMMENT']
        nexp = header['NEXP']
        max_nexp = nexp + extra_nexp
        explist = ExposureList(tiles_file=header['TILES'], max_nexp=max_nexp)
        explist._exposures[:nexp] = hdus['EXPOSURES'].data
        explist._tiledata[:] = hdus['TILEDATA'].data
    log = desiutil.log.get_logger()
    log.info('Loaded {} exposures from {}'.format(nexp, name))
    if comment:
        log.info('Loaded with comment "{}".'.format(comment))
    return explist
