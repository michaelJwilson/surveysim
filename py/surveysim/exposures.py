"""Record simulated exposures and collect per-tile statistics.
"""
from __future__ import print_function, division, absolute_import

import numpy as np

import astropy.io.fits

import desiutil.log

import desisurvey.tiles


class ExposureList(object):

    def __init__(self, tiles_file=None, max_nexp=60000):
        self.tiles = desisurvey.tiles.get_tiles(tiles_file)
        self._exposures = np.empty(max_nexp, dtype=[
            ('MJD', np.float64),
            ('EXPTIME', np.float32),
            ('TILEID', np.int32),
            ('SNR2FRAC', np.float32),
            ('AIRMASS', np.float32),
            ('SEEING', np.float32),
            ('TRANSP', np.float32),
            ('SKY', np.float32),
        ])
        self.nexp = 0
        self._tiledata = np.empty(self.tiles.ntiles, dtype=[
            ('AVAIL', np.int32),
            ('PLANNED', np.int32),
            ('EXPTIME', np.float32),
            ('SNR2FRAC', np.float32),
            ('NEXP', np.int32)
        ])
        self._tiledata['AVAIL'] = -1
        self._tiledata['PLANNED'] = -1
        self._tiledata['EXPTIME'] = 0.
        self._tiledata['SNR2FRAC'] = 0.
        self._tiledata['NEXP'] = 0
        self.initial_night = None

    def update_tiles(self, night, available, planned):
        """Update tile availability and planning status.

        Parameters
        ----------
        night : datetime.date
            Night of initial observing.
        available : array
            Array of tile indices that are currently available.
        planned : array
            Array of tile indices that are currently planned (priority > 0).
        """
        if self.initial_night is None:
            self.initial_night = night
            night_index = 0
        else:
            night_index = (night - self.initial_night).days
            if night_index < 0:
                raise ValueError('night must be advancing.')
        self._tiledata['AVAIL'][available] = night_index
        self._tiledata['PLANNED'][planned] = night_index

    def add(self, mjd, exptime, tileID, snr2frac, airmass, seeing, transp, sky):
        if self.nexp >= len(self._exposures):
            raise RuntimeError(
                'Need to increase max_nexp={}'.format(len(self._exposures)))
        self._exposures[self.nexp] = (
            mjd, exptime, tileID, snr2frac, airmass, seeing, transp, sky)
        self.nexp += 1
        tileinfo = self._tiledata[self.tiles.index(tileID)]
        tileinfo['EXPTIME'] += exptime
        tileinfo['SNR2FRAC'] = snr2frac
        tileinfo['NEXP'] += 1

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
