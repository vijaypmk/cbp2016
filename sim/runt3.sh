#!/bin/bash

echo "Running make clean."
make clean

echo "Running make."
make

echo "Running Trace 3."
./predictor ../traces/SHORT_MOBILE-3.bt9.trace.gz
