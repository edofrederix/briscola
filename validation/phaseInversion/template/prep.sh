#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

SOLVER=$1
MESHX=$2
NPROCX=$3
PSOLVER=$4
NORMALSCHEME=$5
CURVATURESCHEME=$6

VARS="\
    -DVARMESHX=$MESHX \
    -DVARNPROCX=$NPROCX \
    -DVARPSOLVER=$PSOLVER \
    -DVARNORMALSCHEME=$NORMALSCHEME \
    -DVARCURVATURESCHEME=$CURVATURESCHEME"

m4 $VARS system/briscolaMeshDict.m4 > system/briscolaMeshDict
m4 $VARS system/briscolaTwoPhaseDict.m4 > system/briscolaTwoPhaseDict
m4 $VARS system/briscolaSolverDict.m4 > system/briscolaSolverDict
m4 $VARS system/briscolaSchemeDict.m4 > system/briscolaSchemeDict

cp -r $BRISCOLA/cases/twoPhase/phaseInversion/code .
wclean code && wmake code > log.wmake 2>&1
