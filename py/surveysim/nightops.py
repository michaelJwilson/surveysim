#! /usr/bin/python
from __future__ import print_function, division, absolute_import

from datetime import datetime, timedelta
import os
from shutil import copyfile

import numpy as np

import astropy.time
import astropy.io.fits as pyfits
from astropy.table import Table, vstack
import astropy.units as u

import desiutil.log

from desisurvey.exposurecalc import expTimeEstimator
from desisurvey.nextobservation import nextFieldSelector
import desisurvey.ephemerides
import desisurvey.config

from surveysim.observefield import observeField


LSTres = 1.0/144.0 # Should be the same as in afternoon planner and next field selector
MaxExpLen = 3600.0 # One hour
CRsplit = 1200.0   # 20 minutes
ReadOutTime = 120.0 # Should be the same as in next field selector


class obsCount:
    """
    Counter for observation number.  In real operations, each observation
    will have its own file with its number as part of the filename.
    """

    def __init__(self, start_val=0):
        """
        Initialise the counter to zero
        """
        self.obsNumber = start_val

    def update(self):
        """
        Adds 1 to the counter

        Returns:
            string containing the part of the filename with the observation number
        """
        self.obsNumber += 1
        return '{:08d}'.format(self.obsNumber)


def nightOps(night, date_string, obsplan, weather, ocnt, tilesObserved):
    """
    Carries out observations during one night and writes the output to disk

    Args:
        night: row of tabulated ephmerides data for this night
        date_string: string of the form YYYYMMDD
        obsplan: string, filename of today's afternoon plan
        weather: surveysim.weather.Weather object
        w: dictionnary containing the following keys
           'Seeing', 'Transparency', 'OpenDome', 'Clouds'
        ocnt: obsCount object
        tilesObserved: table with follwing columns: tileID, status

    Returns:
        Updated tilesObserved table
    """
    log = desiutil.log.get_logger()
    config = desisurvey.config.Configuration()

    nightOver = False
    # Start the night during bright twilight.
    mjd = night['brightdusk']
    time = astropy.time.Time(mjd, format='mjd')

    obsList = []

    # Test if the weather permits the dome to open tonight.
    if not weather.get(time)['open']:
        log.info('Bad weather forced the dome to remain shut for the night.')
        return tilesObserved

    slew = False
    ra_prev = 1.0e99
    dec_prev = 1.0e99
    while nightOver == False:
        # Get the current weather conditions.
        conditions = weather.get(time)
        seeing, transparency = conditions['seeing'], conditions['transparency']
        # Select the next target to observe.
        target, setup_time = nextFieldSelector(
            obsplan, mjd, conditions, tilesObserved, slew, ra_prev, dec_prev)
        if target != None:
            # Calculate the target's airmass.
            airmass = desisurvey.utils.get_airmass(
                time, target['RA'] * u.deg, target['DEC'] * u.deg)
            # Calculate the nominal exposure time required for this target.
            exposure = expTimeEstimator(
                seeing, transparency, airmass, target['Program'], target['Ebmv'],
                target['DESsn2'], night['moon_illum_frac'], target['MoonDist'],
                target['MoonAlt'])
            if exposure <= MaxExpLen:
                status, real_exposure, real_sn2 = observeField(target, exposure)
                real_exposure += ReadOutTime * np.floor(real_exposure/CRsplit)
                target['Status'] = status
                target['Exposure'] = real_exposure
                target['obsSN2'] = real_sn2
                mjd += (setup_time + real_exposure)/86400.0
                tilesObserved.add_row([target['tileID'], status])
                slew = True
                ra_prev = target['RA']
                dec_prev = target['DEC']
                # Prepare output table.
                t = astropy.time.Time(mjd, format = 'mjd')
                tbase = str(t.isot)
                obsList.append((target['tileID'],  target['RA'], target['DEC'], target['PASS'], target['Program'], target['Ebmv'],
                               target['maxLen'], target['moon_illum_frac'], target['MoonDist'], target['MoonAlt'], conditions['seeing'], conditions['transparency'],
                               airmass, target['DESsn2'], target['Status'],
                               target['Exposure'], target['obsSN2'], tbase, mjd))
            else:
                # Try another target?
                # Observe longer split into modulo(max_len)
                mjd += LSTres
                slew = False # Can slew to new target while waiting.
        else:
            mjd += LSTres
            slew = False
        # Check time
        if mjd > night['brightdawn']:
            nightOver = True

    if len(obsList) > 0:
        filename = config.get_path('obslist{0}.fits'.format(date_string))
        cols = np.rec.array(obsList,
                           names = ('TILEID  ',
                                    'RA      ',
                                    'DEC     ',
                                    'PASS    ',
                                    'PROGRAM ',
                                    'EBMV    ',
                                    'MAXLEN  ',
                                    'MOONFRAC',
                                    'MOONDIST',
                                    'MOONALT ',
                                    'SEEING  ',
                                    'LINTRANS',
                                    'AIRMASS ',
                                    'DESSN2  ',
                                    'STATUS  ',
                                    'EXPTIME ',
                                    'OBSSN2  ',
                                    'DATE-OBS',
                                    'MJD     '),
                            formats = ['i4', 'f8', 'f8', 'i4', 'a8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'i4', 'f8', 'f8', 'a24', 'f8'])
        tbhdu = pyfits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(filename, clobber=True)
        # This file is to facilitate plotting
        all_path = config.get_path('obslist_all.fits')
        if os.path.exists(all_path):
            obsListOld = Table.read(all_path, format='fits')
            obsListNew = Table.read(filename, format='fits')
            obsListAll = vstack([obsListOld, obsListNew])
            obsListAll.write(all_path, format='fits', overwrite=True)
        else:
            copyfile(filename, all_path)

    return tilesObserved
