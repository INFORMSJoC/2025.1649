#!/bin/bash
cd build
cmake ..
make
./CircleCovering Square 100 1 &
sleep 1
wait