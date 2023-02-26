#!/bin/bash

TEST=IO

FORMATS=(ascii binary)
MESHES=(box pipe bigbox)
NPROCS=(10 5 16)

if [ -f build/Test-$TEST ]; then

    for FORMAT in ${FORMATS[@]}; do

        for I in ${!MESHES[@]}; do

            MESH=${MESHES[I]}
            NPROC=${NPROCS[I]}

            cp system/controlDict.$FORMAT system/controlDict
            cp system/briscolaMeshDict.$MESH system/briscolaMeshDict

            OUTPUT=$(mpirun -np $NPROC --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

            RET=$?

            if [ "$RET" != "0" ]; then

                echo "Test $TEST with format = $FORMAT, mesh = $MESH failed"

            else

                echo "Test $TEST with format = $FORMAT, mesh = $MESH succeeded"

            fi

            rm -fr 1 system/controlDict system/briscolaMeshDict briscola*.vtk.series

        done

    done

else

    echo Test $TEST is not compiled

fi
