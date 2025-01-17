#!/bin/bash

N=${1:-16}
NP=${2:-2}

m4 -DVARN=$N -DVARNP=$NP system/briscolaMeshDict.m4 > system/briscolaMeshDict
