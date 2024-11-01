#!/bin/bash

for TEST in ./*/; do
(
    if [ -f "$TEST/clean.sh" ]; then

        echo $TEST

        cd $TEST
        ./clean.sh

    fi
) &
done

wait
