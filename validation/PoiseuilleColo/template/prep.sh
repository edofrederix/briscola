#!/bin/bash

MESH=$1
NPROCX=$2
NPROCY=$3
IBM=$4
PSOLVER=$5

if [ "$IBM" == "none" ]; then

    MESHX=$MESH
    MESHY=$MESH
    MESHFILE=briscolaMeshDict.noIBM

else

    MESHX=$MESH
    MESHY=$(echo "print(int($MESH/2*3))" | python)
    MESHFILE=briscolaMeshDict.IBM

fi

IBMBC=${IBM}Dirichlet

VARS="\
    -DVARNPROCX=$NPROCX \
    -DVARNPROCY=$NPROCY \
    -DVARMESHX=$MESHX \
    -DVARMESHY=$MESHY \
    -DVARIBMBC=$IBMBC \
    -DVARPSOLVER=$PSOLVER"

m4 $VARS system/$MESHFILE.m4 > system/briscolaMeshDict
m4 $VARS system/briscolaSolverDict.m4 > system/briscolaSolverDict
m4 $VARS 0/U.m4 > 0/U

wmake -silent code > log.wmake 2>&1
