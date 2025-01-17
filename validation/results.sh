#!/bin/bash

find . -mindepth 2 -maxdepth 2 -name "results.csv" -exec sh -c "echo {} && cat {}" \;
