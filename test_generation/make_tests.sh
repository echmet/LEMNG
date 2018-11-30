#! /bin/bash

ECL_PATH="/home/madcat/Devel/ECHMET/ECHMETCoreLibs-bin/lib"
LEMNG_PATH="/home/madcat/Devel/ECHMET/LEMNG-bin/lib"
GEN_PATH="./ref_res_gen"

CMAKE_DIRECTIVES="cmake_dirs.txt"

TEST_RECIPES="$1"
TEST_TARGET_DIR="$2"

if [ -z "${TEST_RECIPES}" ]; then
	echo "No directory with test recipes"
	exit
fi

if [ -z "${TEST_TARGET_DIR}" ]; then
	echo "No target directory"
	exit
fi

if [ -f "$CMAKE_DIRECTIVES" ]; then
	rm "$CMAKE_DIRECTIVES"
fi

for f in "$TEST_RECIPES"/*.json; do
	[ -e "$f" ] || continue;

	BASENAME=${f##*/}
	BASENAME=${BASENAME%.json}

	./test_generator.py --ECL_path "${ECL_PATH}" --LEMNG_path "${LEMNG_PATH}" --generator_path "${GEN_PATH}" \
		--input "${f}" --output	"${TEST_TARGET_DIR}/${BASENAME}_nois.cpp" --silent >> "$CMAKE_DIRECTIVES"
	./test_generator.py --ECL_path "${ECL_PATH}" --LEMNG_path "${LEMNG_PATH}" --generator_path "${GEN_PATH}" \
		--input "${f}" --output	"${TEST_TARGET_DIR}/${BASENAME}_is.cpp" --debhue --onsfuo --silent >> "$CMAKE_DIRECTIVES"
done
