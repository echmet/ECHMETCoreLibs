#! /bin/sh
source ./ref_tool_glob.sh

LD_LIBRARY_PATH=${ECL_BIN} ./ref_tool ${1} ${2}
