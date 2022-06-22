#!/bin/bash

cd ${0%/*} || exit 1

. $WM_PROJECT_DIR/bin/tools/CleanFunctions

wclean code
rm -f code/initialCondition

cleanCase
rm -rf 0
