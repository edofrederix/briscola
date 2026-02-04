#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

wclean code

rm -f briscola*.pvd timeData
rm -f *.vtk.series
rm -f slurm*
