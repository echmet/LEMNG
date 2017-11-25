#! /bin/sh
ECL_BIN="/home/madcat/Devel/ECHMET/ECHMETCoreLibs_varReal/bin/lib"

LD_LIBRARY_PATH=../build:${ECL_BIN} lldb -- ./ref_tool ${1} ${2} ${3} ${4} ${5} ${6}
