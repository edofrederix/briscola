#!/bin/bash

CURDIR=$(pwd)

for DIR in $(find . -mindepth 1 -maxdepth 1 -type d); do

    CASE=$(basename $DIR)

    cd $DIR

    rm -fr slurm* runs results.csv

    cd $CURDIR

done
