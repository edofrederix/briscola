#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -fr briscola*.vtk.series *.pdf *.txt \
    system/briscolaTwoPhaseDict \
    system/briscolaSolverDict \
    code/build
