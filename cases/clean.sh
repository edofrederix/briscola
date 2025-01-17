#!/bin/bash

for FILE in $(find . -mindepth 2 -maxdepth 5 -name clean.sh); do
(
    DIR=$(dirname $FILE)

    echo $DIR

    cd $DIR

    ./clean.sh > /dev/null 2>&1
) &
done

wait
