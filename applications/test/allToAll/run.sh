#!/bin/bash

TEST=allToAll

if [ -f build/Test-$TEST ]; then

    OUTPUT=$(mpirun -np 18 --oversubscribe ./build/Test-$TEST -parallel '(32 48 32)' '(9 2 1)' '(3 3 2)' 2 > /dev/null 2>&1)

    RET=$?

    if [ "$RET" != "0" ]; then

        echo Test $TEST failed

    else

        echo Test $TEST succeeded

    fi


else

    echo Test $TEST is not compiled

fi
