# cryptoLp

Implementation for the paper:  
**Secure and Efficient $L^p$-Norm Computation for Two-Party Learning Applications**

## Repository Structure

This implementation extends the [SCI](https://github.com/mpc-msri/EzPC/tree/master/SCI) library. 
The core protocols are in: `SCi/tests`.

- [`aux-protocols.cpp`](SCI/src/BuildingBlocks/aux-protocols.cpp)  
- [`math-functions.cpp`](SCI/src/Math/math-functions.cpp)

Benchmarks can be found in: `SCi/tests`.

## Dependencies

Install the following before building (not handled by `build.sh`):

```bash
sudo apt update
sudo apt install -y build-essential cmake git pkg-config libssl-dev libgmp-dev
```

## Build

```bash
bash build.sh
```

Compiled binaries are placed in `./build/bin`.

## Running Benchmarks

Each benchmark requires two terminals, one per party. For example, to run the Manhattan distance benchmark:

```bash
# Terminal 1
./build/bin/bench_md-OT r=1

# Terminal 2
./build/bin/bench_md-OT r=2
```

## Known Build Issues

| # | Error | Cause | Fix |
|---|-------|-------|-----|
| 1 | `Could NOT find OpenSSL` | Missing Linux path hints | `build.sh` now auto-detects and passes OpenSSL paths to CMake |
| 2 | `Could NOT find GMP` | Missing GMP hints for CMake | `build.sh` now detects and passes GMP/GMPXX paths to CMake |
| 3 | `std::unique_lock` not found in `locks.h` | Missing `<mutex>` include in SEAL | `build.sh` patches SEAL source; also raises `SEAL_POLY_MOD_DEGREE_MAX` to `65536` |
| 4 | `uint32_t was not declared` in cleartext float library | Missing `<cstdint>` in header | `build.sh` adds `#include <cstdint>` to `cleartext_library_float.h` |
| 5 | `Could NOT find Eigen3/SEAL` after cloning to a different path | Absolute machine-specific paths in `SCI/cmake/SCIConfig.cmake` | Dependency hints were changed to clone-location-independent relative paths in `SCI/cmake/SCIConfig.cmake.in` (and regenerated `SCIConfig.cmake`) |
