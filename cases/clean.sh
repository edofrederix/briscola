#!/bin/bash

CURDIR=$(pwd)

for FILE in $(find . -mindepth 3 -maxdepth 3 -name clean.sh); do

    DIR=$(dirname $FILE)
    CASE=$(basename $DIR)

    cd $DIR

    ./clean.sh

    cd $CURDIR

done
