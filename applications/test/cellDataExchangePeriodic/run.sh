#!/bin/bash

TEST=cellDataExchangePeriodic
NPROCS=(1 2 3)

if [ -f build/Test-$TEST ]; then

    for NPROC in ${NPROCS[@]}; do

        m4 -DVARNPROC=$NPROC system/briscolaMeshDict.m4 > \
            system/briscolaMeshDict

        NPROCTOTAL=$(( $NPROC*$NPROC*$NPROC ))

        if [ "$NPROCTOTAL" == "1" ]; then

            OUTPUT=$(./build/Test-$TEST > /dev/null 2>&1)

        else

            OUTPUT=$(mpirun -np $NPROCTOTAL --oversubscribe \
                ./build/Test-$TEST -parallel > /dev/null 2>&1)

        fi

        RET=$?

        if [ "$RET" != "0" ]; then

            echo Test $TEST failed
            exit

        fi

    done

    echo Test $TEST succeeded

else

    echo Test $TEST is not compiled

fi
