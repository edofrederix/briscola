#!/bin/bash

TEST=meshField

MESHES=(box bigbox pipe)
NPROCS=(4 16 5)

if [ -f build/Test-$TEST ]; then

    for I in ${!MESHES[@]}; do

        MESH=${MESHES[I]}
        NPROC=${NPROCS[I]}

        cp system/briscolaMeshDict.$MESH system/briscolaMeshDict

        OUTPUT=$(mpirun -np $NPROC --oversubscribe ./build/Test-$TEST -parallel > /dev/null 2>&1)

        RET=$?

        if [ "$RET" != "0" ]; then

            echo Test $TEST mesh = $MESH failed

        else

            echo Test $TEST mesh = $MESH succeeded

        fi

        rm -f system/briscolaMeshDict

    done

else

    echo Test $TEST is not compiled

fi
