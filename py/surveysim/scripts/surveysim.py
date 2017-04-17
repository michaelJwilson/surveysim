"""Script wrapper for running survey simulations.
"""
from __future__ import print_function, division, absolute_import

import argparse
import datetime
import os

import desiutil.log

import desisurvey.config

import surveysim.simulator


def parse(options=None):
    """Parse command-line options for running survey simulations.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true',
        help='provide verbose output on progress')
    parser.add_argument(
        '--start', type=str, default='2019-08-28', metavar='DATE',
        help='survey starts on the evening of this day, formatted as YYYY-MM-DD')
    parser.add_argument(
        '--stop', type=str, default='2020-07-14', metavar='DATE',
        help='survey stops on the morning of this day, formatted as YYYY-MM-DD')
    parser.add_argument(
        '--seed', type=int, default=123, metavar='N',
        help='random number seed for generating observing conditions')
    parser.add_argument(
        '--use-jpl', action='store_true',
        help='Use JPL Horizons to calculate ephemerides')
    parser.add_argument(
        '--save', default='tiles-observed.fits', metavar='FILENAME',
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
    try:
        args.start = datetime.datetime.strptime(args.start, '%Y-%m-%d').date()
    except ValueError as e:
        raise ValueError('Invalid start-date: {0}'.format(e))
    try:
        args.stop = datetime.datetime.strptime(args.stop, '%Y-%m-%d').date()
    except ValueError as e:
        raise ValueError('Invalid stop-date: {0}'.format(e))
    if args.start >= args.stop:
        raise ValueError('Expected start-date < stop-date.')

    if args.resume is not None and not os.path.exists(args.resume):
        raise ValueError('No resume file found: {0}'.format(args.resume))

    return args


def main(args):
    """Command-line driver for running survey simulations.
    """
    # Set up the logger
    if args.verbose:
        log = desiutil.log.get_logger(desiutil.log.DEBUG)
    else:
        log = desiutil.log.get_logger(desiutil.log.WARNING)

    # Remove any existing obslist_all.fits file since nightops
    # concatenates to it.
    if os.path.exists('obslist_all.fits'):
        os.rename('obslist_all.fits', 'obslist_all_save.fits')
        log.info('Renamed obslist_all.fits to obslist_all_save.fits')

    # Set the output path if requested.
    if args.output_path is not None:
        config = desisurvey.config.Configuration()
        config.set_output_path(args.output_path)

    # Create the simulator.
    simulator = surveysim.simulator.Simulator(
        args.start, args.stop, args.seed, tilesubset=None,
        use_jpl=args.use_jpl, tile_file=args.resume)

    # Simulate each night until the survey is complete or the last
    # day is reached.
    while simulator.next_day():
        pass

    if args.save is not None:
        log.info('Saving observed tiles to {0}'.format(args.save))
        simulator.tilesObserved.write(
            args.save, format='fits', overwrite=True)

    return 0
