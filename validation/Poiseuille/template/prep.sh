#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

MESH=$1
NPROCX=$2
NPROCY=$3
MODE=$4
PSOLVER=$5
RKSCHEME=$6

if [ "$MODE" == "normal" ]; then

    MESHX=$MESH
    MESHY=$MESH
    MESHFILE=briscolaMeshDict.normal
    SOURCE=0.03

elif [ "$MODE" == "mapped" ]; then

    MESHX=$MESH
    MESHY=$MESH
    MESHFILE=briscolaMeshDict.mapped
    SOURCE=0.0

else

    MESHX=$MESH
    MESHY=$(echo "print(int($MESH/2*3))" | python)
    MESHFILE=briscolaMeshDict.IBM
    SOURCE=0.03

fi

cp -r $BRISCOLA/cases/singlePhase/Poiseuille/code .
wmake -a -s code 2>&1 > log.wmake

IBMBC=${MODE}Dirichlet

VARS="\
    -DVARNPROCX=$NPROCX \
    -DVARNPROCY=$NPROCY \
    -DVARMESHX=$MESHX \
    -DVARMESHY=$MESHY \
    -DVARIBMBC=$IBMBC \
    -DVARPSOLVER=$PSOLVER \
    -DVARSOURCE=$SOURCE \
    -DVARRKSCHEME=$RKSCHEME"

m4 $VARS system/$MESHFILE.m4 > system/briscolaMeshDict
m4 $VARS system/briscolaSolverDict.m4 > system/briscolaSolverDict
m4 $VARS system/briscolaSchemeDict.m4 > system/briscolaSchemeDict
m4 $VARS system/briscolaSinglePhaseDict.m4 > system/briscolaSinglePhaseDict
m4 $VARS 0/U.m4 > 0/U
