#!/usr/bin/env bash
git pull
cd build
cmake .. && make
cd ..
./build/particle_filter
