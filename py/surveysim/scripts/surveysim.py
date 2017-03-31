"""Script wrapper for running survey simulations.
"""
from __future__ import print_function, division, absolute_import

import argparse

import desispec.log


def parse(options=None):
    """Parse command-line options for running survey simulations.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true',
        help='provide verbose output on progress')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    return args


def main(args):
    """Command-line driver for running survey simulations.
    """
    # Set up the logger
    if args.verbose:
        log = desispec.log.get_logger(desispec.log.DEBUG)
    else:
        log = desispec.log.get_logger()

    if args.verbose:
        log.info(args)

    return 0
