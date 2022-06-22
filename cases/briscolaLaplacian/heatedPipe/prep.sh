#!/bin/bash

H=0.1
R=0.5
G=0.3535534
L=1.0

NX=32
NY=32
NR=32
NZ=64

m4  -DVARH=$H \
    -DVARR=$R \
    -DVARG=$G \
    -DVARL=$L \
    -DVARNX=$NX \
    -DVARNY=$NY \
    -DVARNR=$NR \
    -DVARNZ=$NZ \
    system/briscolaMeshDict.m4 > system/briscolaMeshDict
