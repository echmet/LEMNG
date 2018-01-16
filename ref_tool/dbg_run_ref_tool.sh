#! /bin/sh
source ./ref_tool_glob.sh

LD_LIBRARY_PATH=${LEMNG_BIN}:${ECL_BIN} lldb -- ./ref_tool ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8}
