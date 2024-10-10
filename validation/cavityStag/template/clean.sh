#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -f briscola*.vtk.series *.pdf \
    system/briscolaMeshDict.m4 \
    system/briscolaMeshDict \
    system/briscolaSchemeDict \
    system/briscolaStaggeredDict \
    system/briscolaSolverDict
