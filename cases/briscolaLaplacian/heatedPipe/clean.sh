#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -f briscola*.pvd timeData

find . -name *.m4 | while read IN; do
    OUT=$(echo $IN | rev | cut -c 4- | rev)
    if [ -f "$OUT" ]; then
        rm $OUT
    fi
done
