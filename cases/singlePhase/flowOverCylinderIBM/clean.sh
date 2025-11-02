#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -f briscola*.pvd timeData code/libcode.so
rm -f *.vtk.series slurm*
