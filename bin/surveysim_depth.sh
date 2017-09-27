#!/bin/bash
#############################################################################
# Simulate an alternate survey strategy that proceeds depth first.
# Note that this is one random realization of the observing conditions.
# Change the random seed for a different realization.
# This will take ~4 hours to run and writes ~5.1G to $DESISURVEY_OUTPUT.
#############################################################################

PLAN_ARGS='--verbose --fa-delay 0m --rules rules-depth.yaml'
SIM_ARGS='--verbose --scores --seed 123 --strategy HA+fallback'

surveyinit --verbose
surveyplan --create ${PLAN_ARGS}
surveysim ${SIM_ARGS}

while :
do
    (surveyplan ${PLAN_ARGS}) || break
    (surveysim --resume ${SIM_ARGS}) || break
done
