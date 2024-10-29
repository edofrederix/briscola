#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

IBM=${1:-false}

case "$IBM" in

    true)

        cp system/briscolaMeshDict.IBM system/briscolaMeshDict

        ;;

    false)

        cp system/briscolaMeshDict.noIBM system/briscolaMeshDict

        ;;

    *)

        echo "Invalid IBM option (true or false)"
        exit

esac

wmake -silent code

