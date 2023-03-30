#!/bin/bash

TEST=IO

FORMATS=(ascii binary)
MESHES=(box pipe pipeGrad bigbox 2x 2y 2z)
NPROCS=(4 5 5 16 4 4 4)

if [ -f build/Test-$TEST ]; then

    for FORMAT in ${FORMATS[@]}; do

        for I in ${!MESHES[@]}; do

            MESH=${MESHES[I]}
            NPROC=${NPROCS[I]}

            cp system/controlDict.$FORMAT system/controlDict
            cp ../meshDicts/briscolaMeshDict.$MESH system/briscolaMeshDict

            OUTPUT=$(mpirun -np $NPROC --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

            RET=$?

            if [ "$RET" != "0" ]; then

                echo "Test $TEST with format = $FORMAT, mesh = $MESH failed"

            else

                echo "Test $TEST with format = $FORMAT, mesh = $MESH succeeded"

            fi

            rm -fr [1-8] system/controlDict system/briscolaMeshDict briscola*.vtk.series

        done

    done

else

    echo Test $TEST is not compiled

fi
