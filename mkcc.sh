#! /bin/sh

CC=clang
CCX=clang++

cmake . -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CCX -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
rm Makefile
rm -fr CMakeFiles
rm CMakeCache.txt
rm cmake_install.cmake
