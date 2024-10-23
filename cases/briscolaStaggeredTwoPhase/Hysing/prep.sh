#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

CASE=${1:-1}

case "$CASE" in

    1)

        m4 -DVARRHO2=100 -DVARMU2=1 -DVARSIGMA=24.5 \
            system/briscolaTwoPhaseDict.m4 > system/briscolaTwoPhaseDict

        m4 -DVARPSOLVER=split \
            system/briscolaSolverDict.m4 > system/briscolaSolverDict

        ;;

    2)

        m4 -DVARRHO2=1 -DVARMU2=0.1 -DVARSIGMA=1.96 \
            system/briscolaTwoPhaseDict.m4 > system/briscolaTwoPhaseDict

        m4 -DVARPSOLVER=MG \
            system/briscolaSolverDict.m4 > system/briscolaSolverDict

        ;;

    *)

        echo "Invalid case"
        exit

esac

echo $CASE > case.txt

wmake -silent code
