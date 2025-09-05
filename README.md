# cryptoLp

This is a proof-of-concept implementation for our paper, 'Secure and Efficient $L^p$-Norm Computation for Two-Party Learning Applications.'

In this version, we have made direct changes to the SCI library. Most of our implementations can be found in the files: [aux-protocols.cpp](SCI/src/BuildingBlocks/aux-protocols.cpp) and [math-functions.cpp](SCI/src/Math/math-functions.cpp). We plan to further organize the repository.

To build the binaries for our benchmarks, run [build.sh](build.sh):
```bash
bash build.sh
```
The binary files are located in `/build/bin`. To run the binaries and test our benchmarks, open two terminals. For instance, to test the Manhattan distance binary `bench_md-OT`, execute `bench_md-OT r=1` in one terminal and `bench_md-OT r=2` in the other, where `r` specifies the party number.
