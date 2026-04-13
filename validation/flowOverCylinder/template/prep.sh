#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

MESH=$1
DSCHEME=$2
GRADSCHEME=$3
NU=$4
RKSCHEME=$5

D=1.0
F=4.0
LI=10.0
LO=40.0
H=10.0
W=0.5
GI=4.0
GO=6.0
GS=4.0

##

cp -r $BRISCOLA/cases/singlePhase/flowOverCylinder/code .
wclean code 2>&1 > log.wmake
wmake code 2>&1 >> log.wmake

##

R=$(echo "print($D/2.0)" | python)
D2=$(echo "print($D*$F)" | python)
R2=$(echo "print($D2/2.0)" | python)
RSQRT=$(echo "print($R/2.0**0.5)" | python)
R2SQRT=$(echo "print($R2/2.0**0.5)" | python)

RHO=$(echo "print($MESH/($R2-$R))" | python)

R2MR=$(echo "print($R2-$R)" | python)
R2SQRT2=$(echo "print(2*$R2SQRT)" | python)
HMR2SQRT=$(echo "print($H-$R2SQRT)" | python)

GIi=$(echo "print(1.0/$GI)" | python)
GOi=$(echo "print(1.0/$GO)" | python)
GSi=$(echo "print(1.0/$GS)" | python)

NR=$(python numberOfCells.py $R2MR $RHO 1)
NQ=$(python numberOfCells.py $R2SQRT2 $RHO 1)
NI=$(python numberOfCells.py $LI $RHO $GI)
NO=$(python numberOfCells.py $LO $RHO $GO)
NS=$(python numberOfCells.py $HMR2SQRT $RHO $GS)

m4  -DVARD=$D \
    -DVARR=$R \
    -DVARD2=$D2 \
    -DVARR2=$R2 \
    -DVARRSQRT=$RSQRT \
    -DVARR2SQRT=$R2SQRT \
    -DVARLI=$LI \
    -DVARLO=$LO \
    -DVARH=$H \
    -DVARW=$W \
    -DVARNR=$NR \
    -DVARNQ=$NQ \
    -DVARNI=$NI \
    -DVARNO=$NO \
    -DVARNS=$NS \
    -DVARGI=$GI \
    -DVARGO=$GO \
    -DVARGS=$GS \
    -DVARGIi=$GIi \
    -DVARGOi=$GOi \
    -DVARGSi=$GSi \
    system/briscolaMeshDict.m4 > system/briscolaMeshDict

m4  -DVARNU=$NU system/briscolaSinglePhaseDict.m4 > system/briscolaSinglePhaseDict

m4  -DVARGRADSCHEME=$GRADSCHEME \
    -DVARDSCHEME=$DSCHEME \
    -DVARRKSCHEME=$RKSCHEME \
    system/briscolaSchemeDict.m4 > system/briscolaSchemeDict

