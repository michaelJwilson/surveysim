# Tutorial Guide to Running DESI Survey Simulations

The instructions below simulate the DESI survey as a sequence of tile exposures and record the observing conditions (airmass, seeing, exposure time) for each observation. Simulations are stochastic since the scheduling algorithms respond to randomly generated weather.

These instructions do not perform fiber assignment or simulate any spectra, but instead provide the necessary scheduling inputs for these tasks.

Please [create an issue](https://github.com/desihub/surveysim/issues/new) with any corrections or suggestions for improvement to this tutorial.

## Install Software

If this is your first exposure to DESI software, [start here](https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted). We use [git for source control](https://desi.lbl.gov/trac/wiki/Computing/UsingGit) and you will need to install the base DESI packages [on your laptop](https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/Laptop)
or else [work at NERSC](https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC).  The instructions below have not yet been tested at NERSC.

The instructions below assume that you have already installed recent versions of the following DESI packages:
- specsim
- desiutil
- desimodel
- desisurvey
- surveysim

If you are not sure about how to do this or just want a specific recipe, the following instructions assume that you have installed the [anaconda scientific python distribution](https://docs.continuum.io/anaconda/install) and will create a new python environment for your DESI work. Start from the directory you wish to install software into, then:
```
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
```

Notes to experts:
- The instructions above assume that you are using the bash shell, and need to be modified slightly for (t)csh.
- The matplotlib and basemap packages are not required to follow the instructions below but are useful for plotting the outputs.

## Setup Environment

Create an output directory to hold all survey planning and simulation outputs:
```
mkdir output
export DESISURVEY=$PWD/output
```
If you followed the installation recipe above then make sure you have activated your `desi` environment with:
```
source activate desi
```
Ensure that your `$DESIMODEL` environment variable points to a valid data directory:
```
ls $DESIMODEL/data
```

## Initialize Survey Planning

Before starting the survey, we precompute some tabulated planning data using:
```
surveyinit --verbose
```

This step takes ~35 minutes and writes the following files into output/:
- ephem_2019-08-28_2024-07-13.fits  (~1 min)
- scheduler.fits (~10 mins, ~1.3Gb)
- surveyinit.fits

The first file tabulates ephemerides of the sun, moon and planets.  The second file tabulates the observing efficiency over the footprint and survey duration.  The last file contains optimized hour angle (HA) assignments for each tile and an estimated exposure time.

These files take some time to generate, but are cached and not regenerated after the first time you run this command. If you want to force these files to be recalculated, add the `--recalc` option.

The dates appearing in the ephemerides filename are the nominal start and stop dates of the five-year survey.  These parameters and many others are defined in the [survey configuration](https://github.com/desihub/desisurvey/blob/master/py/desisurvey/data/config.yaml),
which is well commented and provides a good overview of the assumptions used when planning and scheduling observations.

## Create Initial Observing Plan

Next, we assign targets to each first-layer tile and determine the initial observing priorities of each tile using:
```
surveyplan --create --verbose
```
Note that this step does not currently run fiber assignment, but does keep track of which tiles would have been assigned and are available for scheduling.

This step runs quickly and writes the following files into output/:
- plan.fits (~3 mins)
- plan_2019-08-28.fits (backup of plan.fits)
- progress_2019-08-28.fits (empty progress record)

The plan is based on the initial hour-angle assignments computed by `surveyinit` and the observing priority rules specified in `data/rules.yaml` of the `desisurvey` package.  To experiment with different priority rules use, for example:
```
surveyplan --create --verbose --rules $PWD/myrules.yaml
```

## Simulate Initial Observing

```
surveysim --seed 123 --verbose
```

This step takes ~2 minutes and writes the following files into output/:
- weather_123.fits
- stats.fits
- exposures.fits
- progress.fits
- last_date.txt

If you used the random seed above, the survey should stop after simulating
2019-11-30 when "Group 2 Priority 9" (DARK SCG) completes, but this is dependent
on the simulated weather so different seeds will generally give different results.

The generated exposures.fits contains a list of the simulated exposures with
all parameters necessary to simulate spectra (exposure time, airmass, seeing,
moon brightness, etc).

## Iterate Planning and Observing

```
#surveyplan --verbose
surveysim --resume --verbose
```

Each pass of `surveyplan` takes ?? minutes and will write the following files
into output/ where YYYY-MM-DD is the next planned night of observing:
- plan.fits
- plan_YYYY-MM-DD.fits (backup of plan.fits)
- progress_YYYY-MM-DD.fits (backup of progress.fits)

Each pass of `surveysim` simulates one night's observing.  Jobs will write the following files to output/, updating and overwriting the existing files:
- stats.fits
- exposures.fits
- progress.fits
- last_date.txt

Note that the simulated weather is only generated the first time `surveysim` is called and then read by subsequent passes, to ensure a consistent and continuous weather model.

## Automation

You can wrap the commands above into a simple shell script, using the fact
that surveyplan exits with a non-zero error code when it detects that the
simulation has either run of out time or observed all tiles.  For example:
```
surveyinit --verbose
surveyplan --create ${PLAN_ARGS}
surveysim ${SIM_ARGS}

while :
do
    (surveyplan ${PLAN_ARGS}) || break
    (surveysim --resume ${SIM_ARGS}) || break
done
```
Look for complete examples of automation scripts using different scheduling
strategies in the `surveysim/bin/` directory.

## Directory Organization

If you run simulations with different weather (random seed) or scheduling strategies, it is a good idea to keep the outputs separate by redefining the `$DESISURVEY` environment variable. To save some time, you can reuse the ephemerides and scheduler files generated when `surveyplan` is run for the first time in a new output directory, since these do depend on the random seed or survey strategy. For example:
```
mkdir output2
cd output2
ln ../output/ephem_2019-08-28_2024-07-13.fits .
ln ../output/scheduler.fits .
ln ../output/surveyinit.fits .
cd ..
export DESISURVEY=$PWD/output2
```
