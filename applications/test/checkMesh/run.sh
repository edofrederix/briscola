#!/bin/bash

TEST=checkMesh

MESHES=(box bigbox pipe pipeGrad 2x 2y 2z 8)
NPROCS=(4 16 5 5 4 4 4 16)

for I in ${!MESHES[@]}; do

    MESH=${MESHES[I]}
    NPROC=${NPROCS[I]}

    cp ../meshDicts/briscolaMeshDict.$MESH system/briscolaMeshDict

    OUTPUT=$(mpirun -np $NPROC --oversubscribe briscolaCheckMesh -parallel > /dev/null 2>&1)

    RET=$?

    if [ "$RET" != "0" ]; then

        echo Test $TEST mesh = $MESH failed

    else

        echo Test $TEST mesh = $MESH succeeded

    fi

    rm -f system/briscolaMeshDict

done
