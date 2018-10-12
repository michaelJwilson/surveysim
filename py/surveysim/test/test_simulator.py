from __future__ import print_function, division, absolute_import

import datetime
import unittest
import tempfile
import shutil
import os

import numpy as np

import astropy.table

import desimodel.io

import desisurvey.config
import desisurvey.progress
import desisurvey.rules
import desisurvey.plan
import desisurvey.ephemerides
import desisurvey.old.schedule

import surveysim.weather
import surveysim.simulator

class TestSimulator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create a temporary directory.
        cls.tmpdir = tempfile.mkdtemp()
        # Write output files to this temporary directory.
        config = desisurvey.config.Configuration()
        config.set_output_path(cls.tmpdir)

    @classmethod
    def tearDownClass(cls):
        # Remove the directory after the test.
        shutil.rmtree(cls.tmpdir)
        # Reset our configuration.
        desisurvey.config.Configuration.reset()
        desisurvey.utils._dome_closed_fractions = None

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_simulator(self):
        # Initialize reproducible random numbers.
        gen = np.random.RandomState(123)
        # Simulate the first week of the survey.
        config = desisurvey.config.Configuration()
        start = datetime.date(2019,12,1)
        stop = datetime.date(2019,12,8)
        config.first_day.set_value(start)
        config.last_day.set_value(stop)
        # Calculate ephemerides.
        ephem = desisurvey.ephemerides.Ephemerides(use_cache=False)
        # Precompute scheduler data.
        desisurvey.old.schedule.initialize(ephem)
        # Initialize an empty progress record.
        progress = desisurvey.progress.Progress()
        # Simulate weather.
        weather = surveysim.weather.Weather(gen=gen)
        # Initialize a table for efficiency stats tracking.
        stats = astropy.table.Table()
        num_nights = (config.last_day() - config.first_day()).days
        stats['available'] = np.zeros(num_nights)
        stats['overhead'] = np.zeros(num_nights)
        stats['delay'] = np.zeros(num_nights)
        stats['dawn'] = np.zeros(num_nights)
        stats['live'] = np.zeros(num_nights)
        # Initialize an observing plan.
        tiles = astropy.table.Table(desimodel.io.load_tiles(
            onlydesi=True, extra=False))
        design = tiles[['TILEID', 'RA', 'DEC', 'PASS']]
        design['HA'] = np.zeros(len(design))
        design['OBSTIME'] = np.full(len(design), 1000.)
        design.write(config.get_path('surveyinit.fits'))
        rules = desisurvey.rules.Rules()
        priorities = rules.apply(progress)
        plan = desisurvey.plan.create(design['HA'], priorities)
        plan.write(config.get_path('plan.fits'))
        # Run the simulation.
        sim = surveysim.simulator.Simulator(
            start, stop, progress, weather, stats, 'HA', 'plan.fits', gen=gen)
        while sim.next_day():
            pass
        #- Confirm that observations occur only on survey dates.
        obs_table = progress.get_summary('observed')
        for obs in obs_table:
            night = desisurvey.utils.get_date(obs['mjd_min'])
            assert start <= night < stop


def test_suite():
    """Allows testing of only this module with the command::

        python setup.py test -m <modulename>
    """
    return unittest.defaultTestLoader.loadTestsFromName(__name__)
