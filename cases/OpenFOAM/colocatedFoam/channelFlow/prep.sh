#!/bin/bash

cd ${0%/*} || exit 1

. $WM_PROJECT_DIR/bin/tools/RunFunctions

cp -r 0.org 0

wmake -silent code

runApplication blockMesh
runApplication code/initialCondition
runApplication decomposePar
