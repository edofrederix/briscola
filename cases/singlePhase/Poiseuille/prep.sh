#!/bin/bash

if [ -z "$BRISCOLA" ]; then

    echo "BRISCOLA environment variable not set"
    exit

fi

MODE=${1:-normal}

case "$MODE" in

    normal)

        cp system/briscolaMeshDict.normal system/briscolaMeshDict

        ;;

    IBM)

        cp system/briscolaMeshDict.IBM system/briscolaMeshDict

        ;;

    mapped)

        cp system/briscolaMeshDict.mapped system/briscolaMeshDict

        ;;

    *)

        echo "Invalid mode (normal, IBM or mapped)"
        exit

esac

wmake -silent code
