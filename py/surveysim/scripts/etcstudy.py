"""Script to run simulations for exposure-time studies.
"""
from __future__ import print_function, division, absolute_import

import argparse
import sys

import numpy as np

import astropy.table
import astropy.units as u

import specsim.simulator

import desiutil.log


def parse(options=None):
    """Parse command-line options.

    Note that seeing is implicitly fixed by the fiberloss fractions we
    load and cannot be varied in this study.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
        help='display log messages with severity >= info')
    parser.add_argument('--debug', action='store_true',
        help='display log messages with severity >= debug (implies verbose)')
    parser.add_argument(
        '--seed', type=int, default=123, metavar='N',
        help='random number seed for reproducible results')
    parser.add_argument(
        '--load', default=None, metavar='FITS',
        help='name of file containing tabulated simulation parameters to load')
    parser.add_argument(
        '--nspec', type=int, default=50, metavar='NS',
        help='number of spectra to simulate for each observing condition')
    parser.add_argument(
        '--ncond', type=int, default=250, metavar='NC',
        help='number of observing conditions to simulate')
    parser.add_argument(
        '--interval', type=int, default=50, metavar='NP',
        help='log progress messages with this interval')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    return args


def main(args):
    """Command-line driver for running survey simulations.
    """
    if args.load is None:
        print('Missing required --load parameter.')
        sys.exit(-1)

    # Set up the logger
    if args.debug:
        log = desiutil.log.get_logger(desiutil.log.DEBUG)
        args.verbose = True
    elif args.verbose:
        log = desiutil.log.get_logger(desiutil.log.INFO)
    else:
        log = desiutil.log.get_logger(desiutil.log.WARNING)

    # Initialize random numbers.
    gen = np.random.RandomState(args.seed)

    # Load the tabulated simulation parameters to use.
    table = astropy.table.Table.read(args.load)
    z_array = table['z'].copy()
    source_fluxes=np.array(table['source_fluxes'])
    fiber_acceptance_fraction=np.array(table['fiber_acceptance_fraction'])
    log.info('Loaded {} targets with <z> = {:.2f} from {}'
             .format(len(z_array), np.median(z_array), args.load))
    if args.nspec > len(z_array):
        log.warn('Requested nspec > number of spectra available.')
        args.nspec = len(z_array)

    # Initialize the simulator.
    desi = specsim.simulator.Simulator('desi', num_fibers=args.nspec)
    log.info('Initialized specsim for {} fibers'.format(args.nspec))

    # Loop over simulated observing conditions.
    for i in range(args.ncond):

        # Pick a random subset of spectra to simulate.
        idx = gen.permutation(len(z_array))[:args.nspec]

        if i % args.interval == 0:
            log.info('Completed {} / {} simulations.'
                     .format(i + 1, args.ncond))
