#!/bin/bash

SOLVER=$1
MESHX=$2
MESHY=$3
NPROCX=$4
NPROCY=$5
PSOLVER=$6
NORMALSCHEME=$7
CURVATURESCHEME=$8
REDUCEDPRESSURE=$9

VARS="\
    -DVARMESHX=$MESHX \
    -DVARMESHY=$MESHY \
    -DVARNPROCX=$NPROCX \
    -DVARNPROCY=$NPROCY \
    -DVARPSOLVER=$PSOLVER \
    -DVARNORMALSCHEME=$NORMALSCHEME \
    -DVARCURVATURESCHEME=$CURVATURESCHEME \
    -DVARREDUCEDPRESSURE=$REDUCEDPRESSURE"

m4 $VARS system/briscolaMeshDict.m4 > system/briscolaMeshDict
m4 $VARS system/briscolaTwoPhaseDict.m4 > system/briscolaTwoPhaseDict
m4 $VARS system/briscolaSolverDict.m4 > system/briscolaSolverDict
m4 $VARS system/briscolaSchemeDict.m4 > system/briscolaSchemeDict

if [ -d "code.$SOLVER" ]; then

    cp -r code.$SOLVER code

fi

wmake -silent code > log.wmake 2>&1
