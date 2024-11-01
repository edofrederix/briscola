#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

wclean code

rm -f briscola*.vtk.series *.pdf *.txt code/libcode.so \
    system/briscolaMeshDict \
    system/briscolaSolverDict \
    system/briscolaSchemeDict \
    0/U
