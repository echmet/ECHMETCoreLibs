#! /bin/sh

CXX=g++
ECL_PATH="/home/madcat/Devel/ECHMET/ECHMETCoreLibs-bin"

${CXX} cpuCheck.cpp -o cpuCheck \
       -std=c++14 \
       -DECHMET_COMPILER_GCC_LIKE \
       -I"${ECL_PATH}/include/ECHMET/CoreLibs" \
       -L"${ECL_PATH}/lib" \
       -lECHMETShared
