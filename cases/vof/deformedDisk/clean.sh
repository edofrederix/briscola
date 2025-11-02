#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

wclean code

rm -f briscola*.vtk.series *.pdf code/libcode.so
