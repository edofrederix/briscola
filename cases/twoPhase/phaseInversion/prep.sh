#!/bin/bash
 
if [ -z "$BRISCOLA" ]; then
 
    echo "BRISCOLA environment variable not set"
    exit
 
fi
 
wmake -silent code