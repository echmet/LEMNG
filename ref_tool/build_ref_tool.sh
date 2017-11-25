#! /bin/sh

LIBJANSSON_INCLUDE="/home/madcat/Devel/ECHMET/jansson-bin/include"
LIBJANSSON_BIN="/home/madcat/Devel/ECHMET/jansson-bin/lib/libjansson.a"
ECL_INCLUDE="/home/madcat/Devel/ECHMET/ECHMETCoreLibs_varReal/bin/include/ECHMET/CoreLibs"
ECL_BIN="/home/madcat/Devel/ECHMET/ECHMETCoreLibs_varReal/bin/lib/"

clang -c jsonloader/constituents_json_ldr.c \
	-I${LIBJANSSON_INCLUDE}
clang++ -std=c++98 -Wall -Wextra -pedantic -g -O0 \
	ref_tool.cpp json_input_processor.cpp jsonloader/inputreader.cpp \
	constituents_json_ldr.o ${LIBJANSSON_BIN} \
	-o ref_tool \
	-I../include \
	-I${ECL_INCLUDE} \
	-DECHMET_COMPILER_GCC_LIKE \
	-L${ECL_BIN} \
	-L../build \
	-lECHMETShared -lSysComp \
	-lLEMNG
	
