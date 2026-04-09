# cryptoLp

This is the proof-of-concept implementation for our paper:  
**Secure and Efficient $L^p$-Norm Computation for Two-Party Learning Applications**

## Repository Structure

This version modifies the [SCI](https://github.com/mpc-msri/EzPC/tree/master/SCI) library directly. Most of our implementation can be found in:

- [`aux-protocols.cpp`](SCI/src/BuildingBlocks/aux-protocols.cpp)  
- [`math-functions.cpp`](SCI/src/Math/math-functions.cpp)

We plan to further refactor and organize the repository.

## Dependencies

`build.sh` does not install OpenSSL or GMP automatically. Install dependencies before building.

### Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y \
	build-essential cmake git pkg-config \
	libssl-dev libgmp-dev
```

### macOS

```bash
brew install cmake git openssl@3 gmp
```

## Build Instructions

To compile the benchmark binaries:

```bash
bash build.sh
```

The compiled binaries will be located in `./build/bin`.

## Running Benchmarks

To run a benchmark (e.g., for Manhattan distance with OT), open **two terminals** and execute the following:

**Terminal 1:**

```bash
./build/bin/bench_md-OT r=1
```

**Terminal 2:**

```bash
./build/bin/bench_md-OT r=2
```

Here, `r=1` and `r=2` indicate the party number (party 1 and party 2).

## Common Build Issues

This section records build failures encountered while running `bash build.sh`, along with the corresponding fixes.

1. **CMake could not find OpenSSL on Linux**
	- **Error**: `Could NOT find OpenSSL (missing: OPENSSL_CRYPTO_LIBRARY OPENSSL_INCLUDE_DIR)`
	- **Cause**: The build script only configured explicit OpenSSL paths on macOS, and Linux path hints were missing.
	- **Fix**: Updated `build.sh` to auto-detect OpenSSL on Linux and pass `OPENSSL_ROOT_DIR`, `OPENSSL_INCLUDE_DIR`, `OPENSSL_CRYPTO_LIBRARY`, and `OPENSSL_SSL_LIBRARY` to CMake.

2. **CMake could not find GMP**
	- **Error**: `Could NOT find GMP (missing: GMP_INCLUDE_DIR GMP_LIBRARIES)`
	- **Cause**: CMake was not given explicit include/library hints for GMP on Linux.
	- **Fix**: Updated `build.sh` to detect GMP (`gmp.h`, `libgmp.so`, `libgmpxx.so`) from common locations and pass `GMP_INCLUDE_DIR`, `GMP_LIBRARIES`, and `GMPXX_LIBRARIES` to CMake.

3. **SEAL failed to compile with GCC (`std::unique_lock` not found)**
	- **Error**: In `locks.h`, `std::unique_lock<std::shared_mutex>` failed due to missing `<mutex>` include.
	- **Cause**: SCI's existing SEAL patch was not reliably applied in this build flow because the dependency was cloned directly by `build.sh`.
	- **Fix**: Added idempotent SEAL source patching in `build.sh`:
	  - Adds `#include <mutex>` to `SCI/extern/SEAL/native/src/seal/util/locks.h` when needed.
	  - Updates `SEAL_POLY_MOD_DEGREE_MAX` from `32768` to `65536` in `SCI/extern/SEAL/native/src/seal/util/defines.h` when needed.

4. **`uint32_t` not declared while compiling cleartext float library**
	- **Error**: Multiple errors in `SCI/src/cleartext_library_float.cpp` beginning with `uint32_t was not declared in this scope`.
	- **Cause**: `SCI/src/cleartext_library_float.h` uses fixed-width integer types without including `<cstdint>`.
	- **Fix**: Added `#include <cstdint>` to `SCI/src/cleartext_library_float.h`.

After applying the fixes above, `bash build.sh` completes successfully and installs artifacts into `./build/install`.