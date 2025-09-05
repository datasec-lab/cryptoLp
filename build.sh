#!/bin/bash

# Check if current directory is a git repository
if git rev-parse --git-dir > /dev/null 2>&1; then
    echo "This directory is a Git repository."
else
    echo "This directory is NOT a Git repository."
    echo "Initializing a new Git repository and setting up submodules..."

    # Initialize a new Git repository
    git init

    # Add submodules
    git submodule add https://github.com/microsoft/SEAL.git SCI/extern/SEAL
    git submodule add https://gitlab.com/libeigen/eigen.git SCI/extern/eigen

    # Initialize and update submodules recursively
    git submodule update --init --recursive

    # Checkout specific commits
    pushd SCI/extern/SEAL && git checkout 1d5c816 && popd
    pushd SCI/extern/eigen && git checkout 603e213 && popd

    echo "Submodules are set up and checked out to specified commits."
fi

mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=./install ..
cmake --build . --target install --parallel
