#!/bin/bash

APP=directSolverTiming

if [ -f $APP ]; then

    ./prep.sh

    mpirun -np 8 --oversubscribe $APP -parallel > log

else

    echo Test $APP is not compiled

fi
