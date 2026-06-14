#!/bin/bash

set -e

VERSION="1.0.2"
ARCHS="${CMAKE_OSX_ARCHITECTURES:-arm64;x86_64}"
PACKAGE_DIR="nlm"
PACKAGE_ZIP="nlm_v${VERSION}.zip"

echo "Building project for ${ARCHS}..."
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DUSE_SYMLINKS=OFF \
  -DCMAKE_OSX_ARCHITECTURES="${ARCHS}"
cmake --build build -j

echo "Refreshing package contents..."
cmake -E rm -rf "${PACKAGE_DIR}"
cmake -E make_directory "${PACKAGE_DIR}"
cmake -E copy_directory externals "${PACKAGE_DIR}/externals"
cmake -E copy_directory extras "${PACKAGE_DIR}/extras"
cmake -E copy_directory docs "${PACKAGE_DIR}/docs"
cmake -E copy_directory misc "${PACKAGE_DIR}/misc"
cmake -E copy_directory help "${PACKAGE_DIR}/help"
cmake -E copy build/package-info.json README.md LICENSE "${PACKAGE_DIR}"

echo "Creating ${PACKAGE_ZIP}..."
rm -f "${PACKAGE_ZIP}"
COPYFILE_DISABLE=1 zip -r "${PACKAGE_ZIP}" "${PACKAGE_DIR}" \
  -x "*/.DS_Store" \
  -x "__MACOSX/*"
