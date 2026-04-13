#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

SOLVER=$1
MESHX=$2
MESHY=$3
NPROCX=$4
NPROCY=$5
PSOLVER=$6
NORMALSCHEME=$7
CURVATURESCHEME=$8
RKSCHEME=$9
USOLVER=${10}

VARS="\
    -DVARMESHX=$MESHX \
    -DVARMESHY=$MESHY \
    -DVARNPROCX=$NPROCX \
    -DVARNPROCY=$NPROCY \
    -DVARPSOLVER=$PSOLVER \
    -DVARNORMALSCHEME=$NORMALSCHEME \
    -DVARCURVATURESCHEME=$CURVATURESCHEME \
    -DVARRKSCHEME=$RKSCHEME \
    -DVARUSOLVER=$USOLVER"

m4 $VARS system/briscolaMeshDict.m4 > system/briscolaMeshDict
m4 $VARS system/briscolaTwoPhaseDict.m4 > system/briscolaTwoPhaseDict
m4 $VARS system/briscolaSolverDict.m4 > system/briscolaSolverDict
m4 $VARS system/briscolaSchemeDict.m4 > system/briscolaSchemeDict

cp -r $BRISCOLA/cases/twoPhase/Hysing/code .
wclean code > log.wmake 2>&1
wmake code >> log.wmake 2>&1
