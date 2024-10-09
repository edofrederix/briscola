#!/bin/bash

MESH=${1:-uniform}

if [ -f "system/briscolaMeshDict.$MESH" ]; then

    cp system/briscolaMeshDict.$MESH system/briscolaMeshDict

else

    echo "Invalid mesh"

fi
