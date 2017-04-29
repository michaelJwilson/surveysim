#! /usr/bin/python
from __future__ import print_function, division, absolute_import

import numpy as np
from datetime import datetime, timedelta
import os
from shutil import copyfile
from astropy.time import Time
import astropy.io.fits as pyfits
from astropy.table import Table, vstack
from desisurvey.exposurecalc import expTimeEstimator, airMassCalculator
from desisurvey.utils import mjd2lst
from desisurvey.nextobservation import nextFieldSelector
from surveysim.observefield import observeField
import desisurvey.ephemerides
import desisurvey.config
import desiutil.log


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

def nightOps(day_stats, date_string, obsplan, w, ocnt, tilesObserved,
             tableOutput=True):
    """
    Carries out observations during one night and writes the output to disk

    Args:
        day_stats: row of tabulated ephmerides data for today
        date_string: string of the form YYYYMMDD
        obsplan: string, filename of today's afternoon plan
        w: dictionnary containing the following keys
           'Seeing', 'Transparency', 'OpenDome', 'Clouds'
        ocnt: obsCount object
        tilesObserved: table with follwing columns: tileID, status
        tableOutput: bool, if True writes a table of all the night's observations
                     instead of one file per observation.

    Returns:
        Updated tilesObserved table
    """
    log = desiutil.log.get_logger()
    config = desisurvey.config.Configuration()

    nightOver = False
    # Start the night during bright twilight.
    mjd = day_stats['brightdusk']

    if tableOutput:
        obsList = []
    else:
        os.mkdir(date_string)

    conditions = w.getValues(mjd)
    f = open(config.get_path("nightstats.dat"), "a+")
    if conditions['OpenDome']:
        wcondsstr = "1 " + str(conditions['Seeing']) + " " + str(conditions['Transparency']) + " " + str(conditions['Clouds']) + "\n"
        f.write(wcondsstr)
    else:
        wcondsstr = "0 " + str(conditions['Seeing']) + " " + str(conditions['Transparency']) + " " + str(conditions['Clouds']) + "\n"
        f.write(wcondsstr)
    f.close()
    if conditions['OpenDome'] == False:
        log.info("Bad weather forced the dome to remain shut for the night.")
    else:
        log.info('Dome open conditions: seeing {0:.3f}", transparency {1:.3f}, '
                 .format(conditions['Seeing'], conditions['Transparency']) +
                 'cloud {0:.1f}%'.format(100 * conditions['Clouds']))

        slew = False
        ra_prev = 1.0e99
        dec_prev = 1.0e99
        while nightOver == False:
            conditions = w.updateValues(conditions, mjd)

            lst = mjd2lst(mjd)
            target, setup_time = nextFieldSelector(
                obsplan, mjd, conditions, tilesObserved, slew,
                ra_prev, dec_prev)
            if target != None:
                # Compute mean to apparent to observed ra and dec???
                airmass, tile_alt, tile_az = airMassCalculator(
                    target['RA'], target['DEC'], lst, return_altaz=True)
                exposure = expTimeEstimator(conditions, airmass, target['Program'], target['Ebmv'], target['DESsn2'], day_stats['moon_illum_frac'], target['MoonDist'], target['MoonAlt'])
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
                    if tableOutput:
                        t = Time(mjd, format = 'mjd')
                        tbase = str(t.isot)
                        obsList.append((target['tileID'],  target['RA'], target['DEC'], target['PASS'], target['Program'], target['Ebmv'],
                                       target['maxLen'], target['moon_illum_frac'], target['MoonDist'], target['MoonAlt'], conditions['Seeing'], conditions['Transparency'],
                                       airmass, target['DESsn2'], target['Status'],
                                       target['Exposure'], target['obsSN2'], tbase, mjd))
                    else:
                        # Output headers, but no data.
                        # In the future: GFAs (i, x, y + metadata for i=id, time, postagestampid) and fiber positions.
                        prihdr = pyfits.Header()
                        prihdr['TILEID  '] = target['tileID']
                        prihdr['RA      '] = target['RA']
                        prihdr['DEC     '] = target['DEC']
                        prihdr['PROGRAM '] = target['Program']
                        prihdr['EBMV    '] = target['Ebmv']
                        prihdr['MAXLEN  '] = target['maxLen']
                        prihdr['MOONFRAC'] = target['moon_illum_frac']
                        prihdr['MOONDIST'] = target['MoonDist']
                        prihdr['MOONALT '] = target['MoonAlt']
                        prihdr['SEEING  '] = conditions['Seeing']
                        prihdr['LINTRANS'] = conditions['Transparency']
                        prihdr['AIRMASS '] = airmass
                        prihdr['DESSN2  '] = target['DESsn2']
                        prihdr['STATUS  '] = target['Status']
                        prihdr['EXPTIME '] = target['Exposure']
                        prihdr['OBSSN2  '] = target['obsSN2']
                        t = Time(mjd, format = 'mjd')
                        tbase = str(t.isot)
                        nt = len(tbase)
                        prihdr['DATE-OBS'] = tbase
                        prihdr['MJD     '] = mjd
                        filename = config.get_path(
                            '{0}/desi-exp-{1}.fits'
                            .format(date_string, ocnt.update()))
                        prihdu = pyfits.PrimaryHDU(header=prihdr)
                        prihdu.writeto(filename, clobber=True)
                else:
                    # Try another target?
                    # Observe longer split into modulo(max_len)
                    mjd += LSTres
                    slew = False # Can slew to new target while waiting.
            else:
                mjd += LSTres
                slew = False
            # Check time
            if mjd > day_stats['brightdawn']:
                nightOver = True

    if tableOutput and len(obsList) > 0:
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
