#! /bin/sh
source ../ref_tool/ref_tool_glob.sh

clang -c ../ref_tool/jsonloader/constituents_json_ldr.c \
	-I${LIBJANSSON_INCLUDE}
clang++ -std=c++98 -Wall -Wextra -pedantic -g -O0 \
	pansargryning_resgen.cpp ../ref_tool/json_input_processor.cpp ../ref_tool/jsonloader/inputreader.cpp \
	constituents_json_ldr.o ${LIBJANSSON_BIN} \
	-o pgn_resgen \
	-I${LEMNG_INCLUDE} \
	-I${ECL_INCLUDE} \
	-I../ref_tool \
	-DECHMET_COMPILER_GCC_LIKE \
	-Wl,-rpath,${ECL_BIN} \
	-L${ECL_BIN} \
	-L${LEMNG_BIN} \
	-lECHMETShared -lSysComp \
	-lLEMNG
	
