#!/bin/bash

# Exit on error
set -e

# Create build directory and compile
echo "Building project..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_SYMLINKS=ON
cmake --build build -j