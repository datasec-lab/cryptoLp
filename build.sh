#!/bin/bash

set -e  # Stop on first error

# --------- CONFIG ---------
SEAL_COMMIT=1d5c816
EIGEN_COMMIT=603e213
SEAL_DIR=SCI/extern/SEAL
EIGEN_DIR=SCI/extern/eigen
# --------------------------

# Ensure extern directory exists
mkdir -p SCI/extern

clone_if_missing() {
    local dir=$1
    local repo=$2
    local commit=$3

    if [ ! -d "$dir/.git" ]; then
        echo "Cloning $repo into $dir..."
        git clone "$repo" "$dir"
        pushd "$dir"
        git checkout "$commit"
        popd
    else
        echo "Repo already cloned: $dir"
    fi
}

clone_if_missing "$SEAL_DIR" https://github.com/microsoft/SEAL.git "$SEAL_COMMIT"
clone_if_missing "$EIGEN_DIR" https://gitlab.com/libeigen/eigen.git "$EIGEN_COMMIT"

# --------- OPENSSL CONFIG ---------
if [[ "$OSTYPE" == "darwin"* ]]; then
    if ! command -v brew &>/dev/null; then
        echo "Homebrew not found. Please install OpenSSL manually or install Homebrew."
        exit 1
    fi
    OPENSSL_ROOT=$(brew --prefix openssl@3)
    OPENSSL_FLAGS="-DOPENSSL_ROOT_DIR=$OPENSSL_ROOT -DOPENSSL_INCLUDE_DIR=$OPENSSL_ROOT/include -DOPENSSL_CRYPTO_LIBRARY=$OPENSSL_ROOT/lib"
else
    OPENSSL_FLAGS=""
fi
# ----------------------------------

# Build
mkdir -p build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=./install $OPENSSL_FLAGS
cmake --build . --target install --parallel
