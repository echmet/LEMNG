#! /bin/sh
source ./ref_tool_glob.sh

clang -c jsonloader/constituents_json_ldr.c \
	-I${LIBJANSSON_INCLUDE}
clang++ -std=c++98 -Wall -Wextra -pedantic -g -O0 \
	ref_tool.cpp json_input_processor.cpp jsonloader/inputreader.cpp \
	constituents_json_ldr.o ${LIBJANSSON_BIN} \
	-o ref_tool \
	-I${LEMNG_INCLUDE} \
	-I${ECL_INCLUDE} \
	-DECHMET_COMPILER_GCC_LIKE \
	-Wl,-rpath,${ECL_BIN} \
	-L${ECL_BIN} \
	-L${LEMNG_BIN} \
	-lECHMETShared -lSysComp \
	-lLEMNG
	
