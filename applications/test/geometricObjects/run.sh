#!/bin/bash

TEST=geometricObjects

if [ -f build/Test-$TEST ]; then

    OUTPUT=$(./build/Test-$TEST > /dev/null 2>&1)

    RET=$?

    if [ "$RET" != "0" ]; then

        echo Test $TEST failed

    else

        echo Test $TEST succeeded

    fi


else

    echo Test $TEST is not compiled

fi
