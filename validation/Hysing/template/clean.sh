#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -f briscola*.vtk.series *.pdf data.txt result.txt

rm -fr code
wclean code.briscolaColocatedTwoPhase
wclean code.briscolaStaggeredTwoPhase

rm -f \
    system/briscolaMeshDict \
    system/briscolaTwoPhaseDict \
    system/briscolaSolverDict
