#!/bin/bash

# Add and initialize SEAL submodule if missing
if [ ! -d "SCI/extern/SEAL" ]; then
    echo "Adding SEAL submodule..."
    git submodule add https://github.com/microsoft/SEAL.git SCI/extern/SEAL
    pushd SCI/extern/SEAL && git checkout 1d5c816 && popd
else
    echo "SEAL submodule already exists."
fi

# Add and initialize Eigen submodule if missing
if [ ! -d "SCI/extern/eigen" ]; then
    echo "Adding Eigen submodule..."
    git submodule add https://gitlab.com/libeigen/eigen.git SCI/extern/eigen
    pushd SCI/extern/eigen && git checkout 603e213 && popd
else
    echo "Eigen submodule already exists."
fi

# Update all submodules recursively
git submodule update --init --recursive


mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=./install ..
cmake --build . --target install --parallel
