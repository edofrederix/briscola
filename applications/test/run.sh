#!/bin/bash

# By default run 16 tests simultaneously. May be overwritten by providing an
# argument.

N=${1:-16}

for DIR in $(find . -maxdepth 1 -mindepth 1 -type d | sort); do

    # Wait as long as N jobs are running in the background

    while [ $(jobs | wc -l) -ge $N ]; do

        sleep 1

    done

    if [ -f "$DIR/build.sh" ]; then
    (
        echo Test $DIR

        cd $DIR

        ./build.sh
        ./run.sh

        cd ..
    ) &
    fi

done

wait
