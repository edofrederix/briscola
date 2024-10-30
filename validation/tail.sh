#!/bin/bash

tail -f $(find . -mindepth 2 -maxdepth 2 -name "results.csv")
