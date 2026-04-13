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

# Apply SCI compatibility fixes to SEAL source (idempotent).
SEAL_DEFINES="$SEAL_DIR/native/src/seal/util/defines.h"
SEAL_LOCKS="$SEAL_DIR/native/src/seal/util/locks.h"

if grep -q "SEAL_POLY_MOD_DEGREE_MAX 32768" "$SEAL_DEFINES"; then
    echo "Patching SEAL poly modulus bound"
    sed -i 's/SEAL_POLY_MOD_DEGREE_MAX 32768/SEAL_POLY_MOD_DEGREE_MAX 65536/' "$SEAL_DEFINES"
fi

if grep -q "#include <shared_mutex>" "$SEAL_LOCKS" && ! grep -q "#include <mutex>" "$SEAL_LOCKS"; then
    echo "Patching SEAL locks header for GCC compatibility"
    sed -i '/#include <shared_mutex>/a #include <mutex>' "$SEAL_LOCKS"
fi

# --------- OPENSSL CONFIG ---------
if [[ "$OSTYPE" == "darwin"* ]]; then
    if ! command -v brew &>/dev/null; then
        echo "Homebrew not found. Please install OpenSSL manually or install Homebrew."
        exit 1
    fi
    OPENSSL_ROOT=$(brew --prefix openssl@3)
    OPENSSL_FLAGS="-DOPENSSL_ROOT_DIR=$OPENSSL_ROOT -DOPENSSL_INCLUDE_DIR=$OPENSSL_ROOT/include -DOPENSSL_CRYPTO_LIBRARY=$OPENSSL_ROOT/lib/libcrypto.dylib -DOPENSSL_SSL_LIBRARY=$OPENSSL_ROOT/lib/libssl.dylib"
else
    OPENSSL_FLAGS=""
    GMP_FLAGS=""

    # Prefer active conda environment when available.
    if [[ -n "$CONDA_PREFIX" ]] && [[ -f "$CONDA_PREFIX/include/openssl/ssl.h" ]]; then
        if [[ -f "$CONDA_PREFIX/lib/libcrypto.so" ]] && [[ -f "$CONDA_PREFIX/lib/libssl.so" ]]; then
            OPENSSL_FLAGS="-DOPENSSL_ROOT_DIR=$CONDA_PREFIX -DOPENSSL_INCLUDE_DIR=$CONDA_PREFIX/include -DOPENSSL_CRYPTO_LIBRARY=$CONDA_PREFIX/lib/libcrypto.so -DOPENSSL_SSL_LIBRARY=$CONDA_PREFIX/lib/libssl.so"
        fi
    fi

    if [[ -n "$CONDA_PREFIX" ]] && [[ -f "$CONDA_PREFIX/include/gmp.h" ]] && [[ -f "$CONDA_PREFIX/lib/libgmp.so" ]]; then
        GMP_FLAGS="-DGMP_INCLUDE_DIR=$CONDA_PREFIX/include -DGMP_LIBRARIES=$CONDA_PREFIX/lib/libgmp.so -DGMPXX_LIBRARIES=$CONDA_PREFIX/lib/libgmpxx.so"
    fi

    # Fallback to common system OpenSSL location on Debian/Ubuntu.
    if [[ -z "$OPENSSL_FLAGS" ]] && [[ -f "/usr/include/openssl/ssl.h" ]] && [[ -f "/usr/lib/x86_64-linux-gnu/libcrypto.so" ]] && [[ -f "/usr/lib/x86_64-linux-gnu/libssl.so" ]]; then
        OPENSSL_FLAGS="-DOPENSSL_ROOT_DIR=/usr -DOPENSSL_INCLUDE_DIR=/usr/include -DOPENSSL_CRYPTO_LIBRARY=/usr/lib/x86_64-linux-gnu/libcrypto.so -DOPENSSL_SSL_LIBRARY=/usr/lib/x86_64-linux-gnu/libssl.so"
    fi

    if [[ -z "$GMP_FLAGS" ]] && [[ -f "/usr/include/gmp.h" ]] && [[ -f "/usr/lib/x86_64-linux-gnu/libgmp.so" ]]; then
        GMP_FLAGS="-DGMP_INCLUDE_DIR=/usr/include -DGMP_LIBRARIES=/usr/lib/x86_64-linux-gnu/libgmp.so"
    fi
fi
# ----------------------------------

# Build
mkdir -p build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=./install $OPENSSL_FLAGS $GMP_FLAGS
cmake --build . --target install --parallel
