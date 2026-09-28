#!/bin/bash

# IBM mode (true or false)
IBM=${1:-false}

case "$IBM" in

    false)

        cp system/briscolaMeshDict.noIBM system/briscolaMeshDict

        ;;

    true)

        cp system/briscolaMeshDict.IBM system/briscolaMeshDict

        ;;

    *)

        echo "Invalid IBM mode (true or false)"
        exit

esac

cmake -S code -B code/build && cmake --build code/build
