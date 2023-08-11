#!/bin/bash

TEST=FFTPoissonSolver

STRETCH=("(1 1 1)" "(2 1 1)" "(1 2 1)" "(1 1 2)")

DECOMP=("(2 2 2)" "(1 2 4)" "(2 1 4)" "(4 2 1)" "(1 1 8)" "(1 8 1)" "(8 1 1)")

if [ -f build/Test-$TEST ]; then

    for s in "${!STRETCH[@]}"
    do
        for d in "${!DECOMP[@]}"
        do
            cat system/briscolaMeshDict.orig > system/briscolaMeshDict

            sed -i "s/VARSTRETCH/${STRETCH[$s]}/" "system/briscolaMeshDict"

            sed -i "s/VARDECOMP/${DECOMP[$d]}/" "system/briscolaMeshDict"

            OUTPUT=$(mpirun -np 8 --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

            RET=$?

            rm -f system/briscolaMeshDict

            NTESTS=$((${#STRETCH[@]}*${#DECOMP[@]}))
            NTEST=$(($s*${#DECOMP[@]}+$d+1))

            if [ "$RET" != "0" ]; then

                echo Test $TEST $NTEST / $NTESTS failed

            else

                echo Test $TEST $NTEST / $NTESTS succeeded

            fi
        done
    done

else

    echo Test $TEST is not compiled

fi
