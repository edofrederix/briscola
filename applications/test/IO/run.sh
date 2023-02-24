#!/bin/bash

TEST=IO

CONTROLDICTS=(system/controlDict.ascii system/controlDict.binary)
MESHDICTS=(system/briscolaMeshDict.box system/briscolaMeshDict.pipe)

if [ -f build/Test-$TEST ]; then

    for CONTROLDICT in ${CONTROLDICTS[@]}; do

        for MESHDICT in ${MESHDICTS[@]}; do

            cp $CONTROLDICT system/controlDict
            cp $MESHDICT system/briscolaMeshDict

            OUTPUT=$(mpirun -np 5 --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

            RET=$?

            if [ "$RET" != "0" ]; then

                echo Test $TEST with $CONTROLDICT $MESHDICT failed

            else

                echo Test $TEST with $CONTROLDICT $MESHDICT succeeded

            fi

            rm -fr 1 system/controlDict system/briscolaMeshDict briscola*.vtk.series

        done

    done

else

    echo Test $TEST is not compiled

fi
