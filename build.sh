#!/bin/bash

# Exit on error
set -e

# Create build directory and compile
echo "Building project..."
cmake -S . -B build -D CMAKE_BUILD_TYPE=Debug
cmake --build build -j