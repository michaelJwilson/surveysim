"""Script wrapper for running survey simulations.
"""
from __future__ import print_function, division, absolute_import

import argparse
import datetime

import desiutil.log

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

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    # Validate start/stop date args and covert to datetime objects.
    try:
        args.start = datetime.datetime.strptime(args.start, '%Y-%m-%d')
    except ValueError as e:
        raise ValueError('Invalid start-date: {0}'.format(e))
    try:
        args.stop = datetime.datetime.strptime(args.stop, '%Y-%m-%d')
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

    simulator = surveysim.simulator.Simulator(
        args.start, args.stop, args.seed, tilesubset=None,
        use_jpl=args.use_jpl, tile_file=args.resume)

    return

    while simulator.next_day():
        pass

    if args.save is not None:
        log.info('Saving observed tiles to {0}'.format(args.save))
        simulator.tilesObserved.write(
            args.save, format='fits', overwrite=True)

    return 0
