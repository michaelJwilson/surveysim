=================================================
Tutorial Guide to Running DESI Survey Simulations
=================================================

Introduction
------------

The instructions below simulate the DESI survey as a sequence of tile exposures
and record the observing conditions (airmass, seeing, exposure time) for each
observation. Simulations are stochastic since the scheduling algorithms respond
to randomly generated weather.

These instructions do not perform fiber assignment or simulate any spectra,
but instead provide the necessary scheduling inputs for these tasks.

Please `create an issue <https://github.com/desihub/surveysim/issues/new>`_
with any corrections or suggestions for improvement to this tutorial.

Install Software
----------------

Requirements
~~~~~~~~~~~~

If this is your first exposure to DESI software,
`start here <https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted>`_.
We use `git for source control <https://desi.lbl.gov/trac/wiki/Computing/UsingGit>`_
and you will need to install the base DESI packages
`on your laptop <https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/Laptop>`_
or else `work at NERSC <https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC>`_.

The following DESI packages must be installed to run this tutorial:

- `specsim <https://github.com/desihub/specsim>`_
- `desiutil <https://github.com/desihub/desiutil>`_
- `desimodel <https://github.com/desihub/desimodel>`_
- `desisurvey <https://github.com/desihub/desisurvey>`_
- `surveysim <https://github.com/desihub/surveysim>`_

In addition, the following non-DESI packages must be installed via pip
since they are not included with the anaconda distribution:

- fitsio
- speclite
- ephem
- healpy

Note that these packages are already included in the custom DESI
anaconda distribution installed at NERSC, so only need to be installed
when running on your laptop or if you need to use non-default versions.

NESRC Installation
~~~~~~~~~~~~~~~~~~

The instructions below were tested on cori in Sep 2017::

    source /project/projectdirs/desi/software/desi_environment.sh
    mkdir -p $SCRATCH/desi/lib/python3.5/site-packages $SCRATCH/desi/bin $SCRATCH/desi/code
    export PYTHONPATH=$SCRATCH/desi/lib/python3.5/site-packages:$PYTHONPATH
    export PATH=$SCRATCH/desi/bin:$PATH
    cd $SCRATCH/desi/code
    for package in desimodel desisurvey surveysim; do
        git clone https://github.com/desihub/$package
        cd $package
        python setup.py develop --prefix $SCRATCH/desi/
        cd ..
    done

Only the first line should be necessary once the desimodel,
desisurvey and surveysim packages have been updated at NERSC.
Note that we use ``$SCRATCH`` for faster I/O but files are periodically
removed.
See `NERSC best practices <https://www.nersc.gov/users/data-analytics/data-analytics-2/python/best-practices/#toc-anchor-3>`_
for details.

Laptop Installation
~~~~~~~~~~~~~~~~~~~

The following instructions assume that you have installed the
`anaconda scientific python distribution <https://docs.continuum.io/anaconda/install>`_
and will create a new python environment for your DESI work.
Start from the directory you wish to install software into, then::

    conda create --name desi pip ipython jupyter numpy scipy astropy pyyaml requests h5py scikit-learn matplotlib basemap
    source activate desi
    pip install fitsio speclite ephem healpy
    for package in specsim desiutil desimodel desisurvey surveysim; do
        git clone https://github.com/desihub/$package
        cd $package
        pip install .
        cd ..
    done
    export DESIMODEL=$PWD/desimodel
    install_desimodel_data -d $DESIMODEL


Notes for experts:

- The instructions above assume that you are using the bash shell, and need to
  be modified slightly for (t)csh.
- The matplotlib and basemap packages are not required to follow the
  instructions below but are useful for plotting the outputs.

Setup Environment
-----------------

In General
~~~~~~~~~~

Create an output directory to hold all survey planning and simulation outputs
and create an environment variable pointing to it.

Ensure that your :envvar:`DESIMODEL` environment variable points to a valid data directory::

    ls $DESIMODEL/data

Also check that the relevant command-line scripts are in your path::

    surveyinit --help
    surveyplan --help
    surveysim --help

Note that all output from these commands goes into :envvar:`DESIMODEL` so they
can be run from any directory and will not write anything to ``$PWD``.

NERSC Environment
~~~~~~~~~~~~~~~~~

Save the output to the ``$SCRATCH`` volume::

    mkdir -p $SCRATCH/desi/output
    export DESISURVEY_OUTPUT=$SCRATCH/desi/output

Laptop Environment
~~~~~~~~~~~~~~~~~~

Enter the parent directory where you will save outputs, then::

    mkdir output
    export DESISURVEY_OUTPUT=$PWD/output

If you followed the installation recipe above then make sure you have activated your ``desi`` environment with::

    source activate desi

Initialize Survey Planning
--------------------------

Before starting the survey, we precompute some tabulated planning data using::

    surveyinit --verbose


This step takes ~35 minutes and writes the following files into output/:

- ephem_2019-12-01_2024-11-30.fits  (~1 min)
- scheduler.fits (~10 mins, ~1.3Gb)
- surveyinit.fits

The first file tabulates ephemerides of the sun, moon and planets.
The second file tabulates the observing efficiency over the footprint
and survey duration.  The last file contains optimized hour angle (HA)
assignments for each tile and an estimated exposure time.

These files take some time to generate, but are cached and not regenerated
after the first time you run this command. If you want to force these files
to be recalculated, add the ``--recalc`` option.

The dates appearing in the ephemerides filename are the nominal start and stop
dates of the five-year survey.  These parameters and many others are defined in
the `survey configuration <https://github.com/desihub/desisurvey/blob/master/py/desisurvey/data/config.yaml>`_,
which is well commented and provides a good overview of the assumptions used
when planning and scheduling observations.

Create Initial Observing Plan
-----------------------------

Next, we assign targets to each first-layer tile and determine the initial
observing priorities of each tile using::

    surveyplan --create --verbose

Note that this step does not currently run fiber assignment, but does keep
track of which tiles would have been assigned and are available for scheduling.

This step runs quickly and writes the following files into output/:

- progress.fits (empty initial progress record)
- plan.fits (~3 mins)
- plan_2019-12-01.fits (backup of plan.fits)

The plan is based on the initial hour-angle assignments computed by ``surveyinit``
and the observing priority rules specified in ``data/rules.yaml`` of the
``desisurvey`` package.  To experiment with different priority rules use, for example::

    surveyplan --create --verbose --rules $PWD/myrules.yaml

Simulate Initial Observing
--------------------------

::

    surveysim --seed 123 --verbose

This step takes ~2 minutes and writes the following files into output/:

- weather_123.fits
- stats.fits
- progress.fits
- last_date.txt

The generated progress.fits records the simulated exposures with all
parameters necessary to simulate spectra (exposure time, airmass, seeing,
moon brightness, etc). It is organized as a per-tile table, but can be
converted to a per-exposure table (which is more convenient for simulation) using::

    from desisurvey.progress import Progress
    Progress(restore='progress.fits').get_exposures().write('exposures.fits')

Refer to the ``get_exposures`` documentation to customize the per-exposure
data that is saved.

Iterate Planning and Observing
------------------------------

::

    #surveyplan --verbose
    surveysim --resume --verbose


Each pass of `surveyplan` takes ?? minutes and will write the following files
into output/ where YYYY-MM-DD is the next planned night of observing:

- plan.fits
- plan_YYYY-MM-DD.fits (backup of plan.fits)

Whenever the priorities change due to a change in the rules state machine,
the corresponding plan is also bookmarked with a symbolic link:

- plan_YYYY-MM-DD_bookmark.fits (symbolic link to plan_YYYY-MM-DD.fits)

Each pass of ``surveysim`` simulates one night's observing.  Jobs will
write the following files to output/, updating and overwriting the existing files:

- stats.fits
- progress.fits
- last_date.txt

Note that the simulated weather is only generated the first time ``surveysim``
is called and then read by subsequent passes, to ensure a consistent and
continuous weather model.

Automation
----------

You can wrap the commands above into a simple shell script, using the fact
that surveyplan exits with a non-zero error code when it detects that the
simulation has either run of out time or observed all tiles.  For example::

    surveyinit --verbose
    surveyplan --create ${PLAN_ARGS}
    surveysim ${SIM_ARGS}

    while :
    do
        (surveyplan ${PLAN_ARGS}) || break
        (surveysim --resume ${SIM_ARGS}) || break
    done

Look for complete examples of automation scripts in the ``surveysim/bin/`` directory.

Visualization
-------------

The `surveymovie` script reads simulation outputs and generates a movie with one
frame per exposure to visualize the scheduler algorithm and survey progress::

    surveymovie --verbose

`An example is available <https://www.youtube.com/watch?v=vO1QZD_aCIo>`_.
A key describing the information displayed in each frame is
`here <https://github.com/desihub/desisurvey/blob/master/doc/img/surveymovie-key.png>`_.
To generate a PNG of a single frame, use::

    surveymovie --expid 123 --save exposure123

to create ``exposure123.png``.

Directory Organization
----------------------

If you run simulations with different weather (random seed) or
scheduling strategies, it is a good idea to keep the outputs separate by
redefining the :envvar:`DESISURVEY_OUTPUT` environment variable. To save some time,
you can reuse the ephemerides and scheduler files generated when ``surveyinit``
is run for the first time in a new output directory, since these do depend on
the random seed or survey strategy. For example::

    mkdir output2
    cd output2
    ln $DESISURVEY_OUTPUT/ephem_2019-12-01_2024-11-30.fits .
    ln $DESISURVEY_OUTPUT/scheduler.fits .
    ln $DESISURVEY_OUTPUT/surveyinit.fits .
    cd ..
    export DESISURVEY_OUTPUT=$PWD/output2

Replace the soft links (``ln``) with copies (``cp``) above unless you will be
keeping the original :envvar:`DESISURVEY_OUTPUT` directory around.

To clean up an output directory before re-running a simulation use::

    rm -f $DESISURVEY_OUTPUT/plan*.fits $DESISURVEY_OUTPUT/scores*.fits $DESISURVEY_OUTPUT/progress*.fits $DESISURVEY_OUTPUT/stats.fits $DESISURVEY_OUTPUT/last_date.txt $DESISURVEY_OUTPUT/weather_*.fits
