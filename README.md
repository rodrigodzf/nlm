<div align="center">

# NLM

</div>

A collection of non-linear modal resonators for Max/MSP.

## Overview

This package provides a set of Max/MSP externals that enable efficient real-time non-linear modal synthesis for strings, membranes, and plates. The externals, implemented in C++, offer interactive control of physical parameters, allow the loading of custom modal data, and provide multichannel output.

Physical modelling synthesis creates sound by simulating the physical properties of real-world objects. The `nlm` externals are based on modal methods, a family of numerical approaches that decompose physical systems into modes that can evolve independently in the linear case and interact in the presence of nonlinearities.

### Objects

- **nlm.plate~** - Non-linear modal plate synthesiser with Berger and von Kármán models
- **nlm.string~** - Non-linear modal string synthesiser (Kirchhoff carrier model)
- **mcs.nlm.plate~** - Multi-channel version of the plate synthesiser
- **mcs.nlm.string~** - Multi-channel version of the string synthesiser

## Installation

> **Note:** These externals have so far only been tested on macOS.

1. Download the latest release from the releases page
2. Unzip the package
3. Place the unzipped folder in your Max Packages folder:
   - macOS: `~/Documents/Max 8/Packages`
   - Windows: `C:\Users\<username>\Documents\Max 8\Packages`
4. Restart Max

## Building from Source

### Requirements

- CMake 3.19 or higher
- C++ compiler with C++17 support
- Max/MSP SDK (included as submodule)

### Install dependencies

```bash
brew install libmatio
brew install eigen
```

### Build Instructions

1. Clone the repository with submodules:

   ```bash
   git clone --recursive https://github.com/rodrigodzf/nlm.git
   ```

2. Create a build directory:

   ```bash
   cd nlm
   mkdir build
   cd build
   ```

3. Configure and build:

   ```bash
   cmake ..
   cmake --build .
   ```

4. Install:

   ```bash
   cmake --install .
   ```
