#!/bin/bash

CURDIR=$(pwd)

for DIR in $(find . -mindepth 1 -maxdepth 1 -type d); do

    CASE=$(basename $DIR)

    cd $DIR

    if [ -f "$CASE.job" ]; then

        sbatch $CASE.job

    fi

    cd $CURDIR

done
