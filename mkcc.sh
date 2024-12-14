#! /bin/sh

function print_usage() {
	echo "Usage:"
	echo "mkcc.sh path_to_eigen"
}

function clean_cmake_stuff() {
	rm Makefile
	rm -fr CMakeFiles
	rm cmake_install.cmake
	rm CTestTestfile.cmake
}

CC=gcc
CCX=g++;

IN_EIGEN_DIR="$1"
if [ -z $IN_EIGEN_DIR ]; then
	print_usage
	exit 1
fi

cmake . \
	-DCMAKE_CXX_COMPILER=$CCX \
	-DCMAKE_BUILD_TYPE=Debug \
	-DEIGEN_INCLUDE_DIR="$IN_EIGEN_DIR" \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=ON

rm CMakeCache.txt
clean_cmake_stuff
