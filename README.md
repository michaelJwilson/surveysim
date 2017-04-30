# DESI Survey Simulations

This package simulates the nightly scheduling of observations during the DESI survey,
using randomly generated observing conditions. For example, to run the first year
of the nominal survey:

    % surveysim --start 2019-08-28 --end 2020-07-13 --verbose

The `--seed` option specifies the simulated weather conditions to use. For details
on all the available options use:

    % surveysim --help

The `surveysim` program reads the following inputs:

 - List of tiles in the DESI footprint from `$DESIMODEL/data/footprint/`) and
 - DESI survey configuration parameters from `data/config.yaml` contained within
   the current `desisurvey` package installation.
 - Optional: Survey planning data from `$DESISURVEY/planner.fits`

The `surveysim` program writes the following files:

 - `ephem_<start>_<stop>.fits`: solar system ephemerides for the simulated period.
 - `weather.fits`: randomly generated weather for the survey duration.
 - `progress.fits`: summary of simulated observations.
 - Optional: `obsplan<date>.fits` afternoon planning file generated for each
   night of simulated observing.

All output files are stored in `$DESISURVEY` by default, but this can be changed
using the `--output-path` option.
