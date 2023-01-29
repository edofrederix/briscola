#!/bin/bash

TEST=allToAll

if [ -f build/Test-$TEST ]; then

    OUTPUT=$(mpirun -np 64 --oversubscribe ./build/Test-$TEST -parallel '(512 512 512)' '(4 4 4)' '(1 8 8)' 2 > /dev/null 2>&1)

    RET=$?

    if [ "$RET" != "0" ]; then

        echo Test $TEST failed

    else

        echo Test $TEST succeeded

    fi


else

    echo Test $TEST is not compiled

fi
