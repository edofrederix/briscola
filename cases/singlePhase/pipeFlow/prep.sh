#!/bin/bash

cmake -S code -B code/build && cmake --build code/build

H=0.3
R=1.0
G=0.7071068
L=15.707963

NX=16
NY=16
NR=32
NZ=256
GRADING=4

GRADINGI=$(echo "print(1.0/$GRADING)" | python)

m4  -DVARH=$H \
    -DVARR=$R \
    -DVARG=$G \
    -DVARL=$L \
    -DVARNX=$NX \
    -DVARNY=$NY \
    -DVARNR=$NR \
    -DVARNZ=$NZ \
    -DVARGRAD=$GRADING \
    -DVARGRADI=$GRADINGI \
    system/briscolaMeshDict.m4 > system/briscolaMeshDict
