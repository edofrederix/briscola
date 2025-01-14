#!/bin/bash

# By default run 4 jobs simultaneously

N=${WM_NCOMPPROCS:-4}

for DIR in ./*/; do

    TEST=$(echo $DIR | cut -c3- | rev | cut -c2- | rev)

    # Wait as long as N jobs are running in the background

    while [ $(jobs | wc -l) -ge $N ]; do

        sleep 1

    done

    if [ -f "$TEST/build.sh" ]; then
    (
        echo Test $TEST

        cd $DIR

        ./build.sh
        ./run.sh

        cd ..
    ) &
    fi

done

wait
