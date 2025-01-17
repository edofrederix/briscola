#!/bin/bash

TEST=FFTPoissonSolver2D

STRETCH=("(1 1 1)" "(2 1 1)" "(1 2 1)")

if [ -f build/Test-$TEST ]; then

    for s in "${STRETCH[@]}"
    do
        cat system/briscolaMeshDict.orig > system/briscolaMeshDict

        sed -i "s/VARSTRETCH/$s/" "system/briscolaMeshDict"

        OUTPUT=$(mpirun --oversubscribe -np 4 ./build/Test-$TEST -parallel > /dev/null 2>&1)

        RET=$?

        rm -f system/briscolaMeshDict

        if [ "$RET" != "0" ]; then

            echo Test $TEST failed

        else

            echo Test $TEST succeeded

        fi
    done

else

    echo Test $TEST is not compiled

fi
