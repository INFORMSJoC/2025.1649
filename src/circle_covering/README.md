# Circle Covering

## Project Overview

Circle Covering is a C++ project for solving region coverage problems with identical circles.
This project is designed for solving both constrained and unconstrained region coverage problems with minimum circles.

Source code is located under `src/`, public headers under `include/`, region data under `data/`, and output result files under `results/`.

## Recommended Usage
A prebuilt executable may already be available in the `build/` directory.

In the environment described below, we recommend building and running this project with CMake. This ensures:

1. Consistency with the compiler and include paths.
2. Proper handling of third-party components such as `clipper` and `cg_descent`.
3. Correct relative-path resolution for `data/` and `results/`.

Direct manual compilation of individual `.cpp` files is not recommended, as it may lead to:

1. Missing include directories
2. Incorrect compiler flags
3. Linking issues with third-party components
4. Inconsistent runtime behavior

If rebuilding is required, please use CMake as described below.

## Prerequisites

- CMake (recommended >= 3.15)
- A C++ compiler supporting C++11 or later (MinGW-w64 `g++`, MSVC, or clang)

This project is organized as a CMake-based project. Use CMake to generate build files and compile the executable.

## Build 

```bash
mkdir build
cd build
cmake ..
cmake --build .
```
You can also build from an IDE ( e.g., Visual Studio).

Note: the CMake target/executable is named `CircleCovering` as defined in `CMakeLists.txt` (see `add_executable`).


## Running the Program

After building, run the generated executable from the `build/` directory so relative paths resolve correctly.

The program uses interactive input: run the executable from the `build/` directory and you will be prompted for the following values:

- Region (e.g. `C14H2`)
- n (number of circles, positive integer)
- CONSTRAINT (0 = unconstrained, 1 = constrained)

Example:

```bash
cd build
./CircleCovering C14H2 50 1
```

When running, the program reads the region data from `../data/<Region>.txt`, executes the search, and output the results under `../results/`.


## Project Layout (short)

- `src/` - source code (entry point at `src/main.cpp`)
- `include/` - public headers and third-party includes (e.g. `clipper/`, `cg_descent/`)
- `data/` - benmark regions
- `results/` - output result files and logs
- `CMakeLists.txt` - CMake build configuration

## Results explanation

- **Results file location**: Result files are stored under `results/<Region>/`, typically with names like `<Region>_<n>_<run>.txt`.
- **File format**: Each result `.txt` file begins with a header line containing two values: `n r`, where `n` is the number of circles and `r` is the computed radius. This is followed by `n` lines of circle center coordinates:

```
n r
x y
x y
...
```

- `results/<Region>/ResultSummary.txt` records summary statistics for the region.
- `results/<Region>/<n>/log.txt` contains per-run logs including run index, radius, and elapsed time.
