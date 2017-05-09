from __future__ import print_function, division, absolute_import

import datetime
import unittest
import os
import shutil
import uuid

import numpy as np

from astropy.table import Table
from astropy import units
from astropy.time import Time

import desiutil.log

import desisurvey.config
import desisurvey.progress


class TestSurveySim(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.origdir = os.getcwd()
        cls.testdir = os.path.abspath('./test-{}'.format(uuid.uuid4()))
        os.mkdir(cls.testdir)
        os.chdir(cls.testdir)
        # Write all outputs to our test path.
        cls.config = desisurvey.config.Configuration()
        cls.config.set_output_path(cls.testdir)

    @classmethod
    def tearDownClass(cls):
        cls.config.reset()
        os.chdir(cls.origdir)
        if os.path.exists(cls.testdir):
            shutil.rmtree(cls.testdir)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_surveysim(self):
        from surveysim.simulator import Simulator
        start = datetime.date(2019,9,1)
        stop = datetime.date(2019,9,8)
        progress = desisurvey.progress.Progress()
        sim = Simulator(start, stop, progress, seed=123456)
        while sim.next_day():
            pass

        #- A plan should exist for every night
        for i in range(1,8):
            self.assertTrue(os.path.exists('obsplan201909{:02d}.fits'.format(i)))

        #- Confirm that observations occur only on survey dates.
        obs_table = progress.get_summary('observed')
        for obs in obs_table:
            night = desisurvey.utils.get_date(obs['mjd_min'])
            assert start <= night < stop


if __name__ == '__main__':
    unittest.main()
