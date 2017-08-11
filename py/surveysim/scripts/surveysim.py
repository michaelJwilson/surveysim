"""Script wrapper for running survey simulations.

This script simulates a sequence of observations until either the nominal
survey end date is reached, or else a predefined fiber-assignment trigger
condition occurs.

See doc/tutorial for an overview of running surveysim and surveyplan to
simulate the full DESI survey.

To run this script from the command line, use the ``surveysim``
entry point that is created when this package is installed, and
should be in your shell command search path.

To profile this script, use, for example::

    import cProfile
    import surveysim.scripts.surveysim as ssim
    cProfile.run("ssim.main(ssim.parse(['--stop', '2019-08-30']))", filename='profile.out')

    import pstats
    p = pstats.Stats('profile.out')
    p.sort_stats('time').print_stats(10)
"""
from __future__ import print_function, division, absolute_import

import argparse
import datetime
import os
import warnings
import sys

import numpy as np

import astropy.table

import desiutil.log

import desisurvey.config
import desisurvey.progress

import surveysim.weather
import surveysim.simulator


def parse(options=None):
    """Parse command-line options for running survey simulations.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
        help='display log messages with severity >= info')
    parser.add_argument('--debug', action='store_true',
        help='display log messages with severity >= debug (implies verbose)')
    parser.add_argument(
        '--start', type=str, default=None, metavar='DATE',
        help='survey starts on the evening of this day, formatted as YYYY-MM-DD')
    parser.add_argument(
        '--stop', type=str, default=None, metavar='DATE',
        help='survey stops on the morning of this day, formatted as YYYY-MM-DD')
    parser.add_argument(
        '--seed', type=int, default=123, metavar='N',
        help='random number seed for generating observing conditions')
    parser.add_argument(
        '--resume', action='store_true',
        help='resuming a previous simulation from its saved progress')
    parser.add_argument(
        '--strategy', default='HA',
        help='Next tile selector strategy to use')
    parser.add_argument(
        '--plan', default='plan.fits',
        help='Name of plan file to use')
    parser.add_argument(
        '--output-path', default=None, metavar='PATH',
        help='Output path where output files should be written')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    # Validate start/stop date args and covert to datetime objects.
    # Unspecified values are taken from our config.
    config = desisurvey.config.Configuration()
    if args.start is None:
        args.start = config.first_day()
    else:
        try:
            args.start = desisurvey.utils.get_date(args.start)
        except ValueError as e:
            raise ValueError('Invalid start: {0}'.format(e))
    if args.stop is None:
        args.stop = config.last_day()
    else:
        try:
            args.stop = desisurvey.utils.get_date(args.stop)
        except ValueError as e:
            raise ValueError('Invalid stop: {0}'.format(e))
    if args.start >= args.stop:
        raise ValueError('Expected start < stop.')

    return args


def main(args):
    """Command-line driver for running survey simulations.
    """
    # Set up the logger
    if args.debug:
        log = desiutil.log.get_logger(desiutil.log.DEBUG)
        args.verbose = True
    elif args.verbose:
        log = desiutil.log.get_logger(desiutil.log.INFO)
    else:
        log = desiutil.log.get_logger(desiutil.log.WARNING)

    # Freeze IERS table for consistent results.
    desisurvey.utils.freeze_iers()

    # Set the output path if requested.
    config = desisurvey.config.Configuration()
    if args.output_path is not None:
        config.set_output_path(args.output_path)

    # Initialize random numbers.
    gen = np.random.RandomState(args.seed)
    weather_name = 'weather_{0}.fits'.format(args.seed)

    if args.resume:
        # Load the previously generated random weather.
        weather = surveysim.weather.Weather(restore=weather_name)
        # Load the progress record to resume writing.
        progress = desisurvey.progress.Progress(restore='progress.fits')
        # Read the stats table that we will update.
        stats = astropy.table.Table.read(config.get_path('stats.fits'))
        # Resume from the last simulated date.
        with open(config.get_path('last_date.txt'), 'r') as f:
            args.start = desisurvey.utils.get_date(f.read().rstrip())
        if args.start >= args.stop:
            log.info('Reached stop date.')
            # Return a shell exit code so scripts can detect this condition.
            sys.exit(9)
    else:
        # Generate and save random weather.
        weather = surveysim.weather.Weather(gen=gen)
        weather.save(weather_name)
        # Create a new (empty) progress record.
        progress = desisurvey.progress.Progress()
        # Initialize a table for efficiency stats tracking.
        stats = astropy.table.Table()
        num_nights = (config.last_day() - config.first_day()).days
        stats['available'] = np.zeros(num_nights)
        stats['overhead'] = np.zeros(num_nights)
        stats['delay'] = np.zeros(num_nights)
        stats['dawn'] = np.zeros(num_nights)
        stats['live'] = np.zeros(num_nights)

    # Create the simulator.
    simulator = surveysim.simulator.Simulator(
        args.start, args.stop, progress, weather, stats,
        args.strategy, args.plan, gen)

    # Simulate one night of observing.
    simulator.next_day()

    # Save the current date.
    with open(config.get_path('last_date.txt'), 'w') as f:
        print(str(simulator.date), file=f)

    # Save the per-night efficiency stats.
    stats.write(config.get_path('stats.fits'), overwrite=True)

    # Save the survey progress after the simulation.
    progress.save('progress.fits')

    # Save the corresponding exposure sequence.
    exposures = progress.get_exposures()
    exposures.write(config.get_path('exposures.fits'), overwrite=True)

    return 0
