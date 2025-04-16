# modalsim

A collection of modal resonators for Max/MSP that simulate resonant objects using modal synthesis techniques.

## Overview

This package provides a set of Max/MSP externals for creating and manipulating modal resonators. Modal synthesis is a technique for simulating the vibrations of objects by modeling them as a collection of resonant modes.

### Objects

- **modal_resonator~** - A basic modal resonator with frequency, decay, and gain controls

## Installation

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

```
brew install libmatio
brew install eigen
```

### Build Instructions

1. Clone the repository with submodules:
   ```
   git clone --recursive https://github.com/yourusername/modalsim.git
   ```

2. Create a build directory:
   ```
   cd modalsim
   mkdir build
   cd build
   ```

3. Configure and build:
   ```
   cmake ..
   cmake --build .
   ```

4. Install:
   ```
   cmake --install .
   ```

## License

This project is licensed under the MIT License - see the LICENSE file for details. 