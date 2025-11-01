#!/bin/bash

MESH=$1
NPROCX=$2
NPROCY=$3
LSCHEME=$4
DSCHEME=$5
NU=$6
COARSEMODE=$7
RKSCHEME=$8

VARS="\
    -DVARNPROCX=$NPROCX \
    -DVARNPROCY=$NPROCY \
    -DVARLSCHEME=$LSCHEME \
    -DVARDSCHEME=$DSCHEME \
    -DVARNU=$NU \
    -DVARCOARSEMODE=$COARSEMODE \
    -DVARRKSCHEME=$RKSCHEME"

cp system/briscolaMeshDict.$MESH.m4 system/briscolaMeshDict.m4

m4 $VARS system/briscolaMeshDict.m4 > system/briscolaMeshDict
m4 $VARS system/briscolaSchemeDict.m4 > system/briscolaSchemeDict
m4 $VARS system/briscolaColocatedDict.m4 > system/briscolaColocatedDict
m4 $VARS system/briscolaStaggeredDict.m4 > system/briscolaStaggeredDict
m4 $VARS system/briscolaSolverDict.m4 > system/briscolaSolverDict
