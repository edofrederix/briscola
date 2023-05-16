#!/bin/bash

TEST=FFTPoissonSolver

STRETCH=("(1 1 1)" "(2 1 1)" "(1 2 1)" "(1 1 2)")

DECOMP=("(2 2 2)" "(1 2 4)" "(2 1 4)" "(4 2 1)" "(1 1 8)" "(1 8 1)" "(8 1 1)")

if [ -f build/Test-$TEST ]; then

    for s in "${STRETCH[@]}"
    do
        for d in "${DECOMP[@]}"
        do
            cat system/briscolaMeshDict.orig > system/briscolaMeshDict

            sed -i "s/VARSTRETCH/$s/" "system/briscolaMeshDict"

            sed -i "s/VARDECOMP/${DECOMP[0]}/" "system/briscolaMeshDict"

            OUTPUT=$(mpirun -np 8 --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

            rm -f system/briscolaMeshDict

            RET=$?

            if [ "$RET" != "0" ]; then

                echo Test $TEST failed

            else

                echo Test $TEST succeeded

            fi
        done
    done

else

    echo Test $TEST is not compiled

fi
