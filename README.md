# cryptoLp

This is the proof-of-concept implementation for our paper:  
**Secure and Efficient $L^p$-Norm Computation for Two-Party Learning Applications**

## Repository Structure

This version modifies the [SCI](https://github.com/mpc-msri/EzPC/tree/master/SCI) library directly. Most of our implementation can be found in:

- [`aux-protocols.cpp`](SCI/src/BuildingBlocks/aux-protocols.cpp)  
- [`math-functions.cpp`](SCI/src/Math/math-functions.cpp)

We plan to further refactor and organize the repository.

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
