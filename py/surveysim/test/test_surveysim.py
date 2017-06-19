from __future__ import print_function, division, absolute_import

import datetime
import unittest
import tempfile
import shutil
import os

import numpy as np

import astropy.table

import desisurvey.config
import desisurvey.progress

import surveysim.weather
import surveysim.simulator

class TestSurveySim(unittest.TestCase):

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

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_simulator(self):
        gen = np.random.RandomState(123)
        start = datetime.date(2019,9,1)
        stop = datetime.date(2019,9,8)
        progress = desisurvey.progress.Progress()
        weather = surveysim.weather.Weather(start, stop, gen=gen)
        # Initialize a table for efficiency stats tracking.
        stats = astropy.table.Table()
        config = desisurvey.config.Configuration()
        num_nights = (config.last_day() - config.first_day()).days
        stats['available'] = np.zeros(num_nights)
        stats['overhead'] = np.zeros(num_nights)
        stats['delay'] = np.zeros(num_nights)
        stats['dawn'] = np.zeros(num_nights)
        stats['live'] = np.zeros(num_nights)
        sim = surveysim.simulator.Simulator(
            start, stop, progress, weather, stats, gen=gen)
        while sim.next_day():
            pass

        #- Confirm that observations occur only on survey dates.
        obs_table = progress.get_summary('observed')
        for obs in obs_table:
            night = desisurvey.utils.get_date(obs['mjd_min'])
            assert start <= night < stop


if __name__ == '__main__':
    unittest.main()
