#!/bin/bash

# Simulation mode (colocated or staggered)
MODE=${1:-colocated}

# Number of cells in the y-direction per brick (i.e., a quarter of the channel
# height)
NY=${2:-16}

##

if [[ ! $NY =~ ^[0-9]+$ ]]; then
    echo "Invalid number of cells"
    exit
fi

DY=$(echo "print(0.025/$NY/2.0)" | python)
DY2=$(echo "print($DY*2.0)" | python)

NX=$(echo "print($NY*4)" | python)
NZ=$(echo "print($NY*2)" | python)

case "$MODE" in

    colocated)

        # Damp the cell centers on both sides of the interface, so we need a
        # thickness of two cells

        m4 -DVARTHICKNESS=$DY2 -DVARNX=$NX -DVARNY=$NY -DVARNZ=$NZ \
            system/briscolaMeshDict.m4 > system/briscolaMeshDict

        ;;

    staggered)

        # Only the y-faces at the interface, so we need a thickness of just one
        # cell

        m4 -DVARTHICKNESS=$DY -DVARNX=$NX -DVARNY=$NY -DVARNZ=$NZ \
            system/briscolaMeshDict.m4 > system/briscolaMeshDict

        ;;

    *)

        echo "Invalid mode (colocated or staggered)"
        exit

esac

echo "Mode = $MODE"

cmake -S code -B code/build && cmake --build code/build
