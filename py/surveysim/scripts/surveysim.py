"""Script wrapper for running survey simulations.

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

import desiutil.log

import desisurvey.config
import desisurvey.progress

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
        '--save', default='progress.fits', metavar='FILENAME',
        help='Name of FITS file where simulated observations are saved')
    parser.add_argument(
        '--resume', default=None, metavar='FILENAME',
        help='Name of saved observations for resuming a simulation')
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
            args.start = datetime.datetime.strptime(args.start, '%Y-%m-%d').date()
        except ValueError as e:
            raise ValueError('Invalid start: {0}'.format(e))
    if args.stop is None:
        args.stop = config.last_day()
    else:
        try:
            args.stop = datetime.datetime.strptime(args.stop, '%Y-%m-%d').date()
        except ValueError as e:
            raise ValueError('Invalid stop: {0}'.format(e))
    if args.start >= args.stop:
        raise ValueError('Expected start < stop.')

    if args.resume is not None and not os.path.exists(args.resume):
        raise ValueError('No resume file found: {0}'.format(args.resume))

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

    # Initialize the survey progress for this simulation.
    progress = desisurvey.progress.Progress(args.resume)

    # Create the simulator.
    simulator = surveysim.simulator.Simulator(
        args.start, args.stop, progress, args.seed)

    # Save simulated weather conditions.
    simulator.weather.save('weather.fits')

    # Simulate each night until the survey is complete or the last
    # day is reached.
    while simulator.next_day():
        pass

    # Save the survey progress after the simulation.
    progress.save(args.save)

    return 0
