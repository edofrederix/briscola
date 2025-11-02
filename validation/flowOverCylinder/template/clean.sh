#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

wclean code

cleanCase

rm -f briscola*.vtk.series system/briscolaMeshDict *.txt *.pdf \
    system/briscolaSinglePhaseDict \
    system/briscolaSchemeDict \
    code/libcode.so
