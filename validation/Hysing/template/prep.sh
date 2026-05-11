#!/bin/bash

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
rm -fr code/build
(cmake -S code -B code/build && cmake --build code/build) > log.cmake 2>&1
