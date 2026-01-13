#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -fr code

rm -f briscola*.vtk.series *.pdf *.txt
