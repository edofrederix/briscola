#!/bin/bash

source $FOAM_SRC/../bin/tools/CleanFunctions

cleanCase

rm -f briscola*.pvd timeData
rm -f *.vtk.series
rm -f slurm*
