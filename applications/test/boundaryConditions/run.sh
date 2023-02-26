#!/bin/bash

TEST=boundaryConditions
MESHES=(2d 3d)
NPROCS=(8 16)

rm -f system/briscolaMeshDict

if [ -f build/Test-$TEST ]; then

    for I in ${!MESHES[@]}; do

        MESH=${MESHES[I]}
        NPROC=${NPROCS[I]}

        cp system/briscolaMeshDict.$MESH system/briscolaMeshDict

        OUTPUT=$(mpirun -np $NPROC --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

        RET=$?

        if [ "$RET" != "0" ]; then

            echo Test $TEST failed

        else

            echo Test $TEST succeeded

        fi

        rm system/briscolaMeshDict

    done

else

    echo Test $TEST is not compiled

fi
