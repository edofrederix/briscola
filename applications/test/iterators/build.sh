#!/bin/bash

if [ ! -d build ]; then
    mkdir build
fi

wmake -silent
