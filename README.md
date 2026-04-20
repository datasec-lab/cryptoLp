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
