cmake_minimum_required(VERSION 3.5)

project(cryptoLp)
set(NAME "cryptoLp")

set(NO_REVEAL_OUTPUT ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(SCI)

add_subdirectory(tests)
