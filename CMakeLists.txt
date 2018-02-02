cmake_minimum_required(VERSION 2.8.12)

project(LEMNG LANGUAGES CXX)

option(SENSITIVE_NUMDERS "Use more precise input concentrations to calculate numerical derivatives at the expense of speed" OFF)
option(PARALLEL_NUM_OPS "Enable parallel numeric operations" ON)
option(DISABLE_TRACING "Disable tracing" OFF)

add_definitions("-D_USE_MATH_DEFINES -DECHMET_DLL_BUILD")

set(LEMNG_THREADING_LIBS)

if (WIN32)
    if (MINGW)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
        add_definitions("-DECHMET_COMPILER_MINGW")
    elseif (MSYS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
        add_definitions("-DECHMET_COMPILER_MSYS")
    elseif (MSVC)
        add_definitions("-DECHMET_COMPILER_MSVC")
    else ()
        message(FATAL "Unsupported or mistedected compiler")
    endif ()
elseif (UNIX)
    set(LEMNG_THREADING_LIBS pthread)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -Wno-c++98-compat -Wno-unused-function -pedantic -fvisibility=hidden")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
    else ()
        message(FATAL "Unsupported or misdetected compiler")
    endif ()
endif ()

if (SENSITIVE_NUMDERS)
    add_definitions("-DECHMET_LEMNG_SENSITIVE_NUMDERS")
endif ()
if (PARALLEL_NUM_OPS)
	add_definitions("-DECHMET_LEMNG_PARALLEL_NUM_OPS")
endIf ()
if (DISABLE_TRACING)
    add_definitions("-DECHMET_TRACER_DISABLE_TRACING")
endif ()

set(libIonProps_SRCS
    src/base_types.cpp
    src/lemng.cpp
    src/calculator_common.cpp
    src/calculator_linear.cpp
    src/calculator_nonlinear.cpp
    src/calculator_matrices.cpp
    src/calculator_types.cpp
    src/efg_plotter.cpp
    src/helpers.cpp
    src/results_maker.cpp)

include_directories(${INCLUDE_DIRECTORIES}
                    SYSTEM "${EIGEN_INCLUDE_DIR}"
                    SYSTEM "${ECHMET_CORE_LIBS_DIR}/include/ECHMET/CoreLibs"
                    "${CMAKE_CURRENT_SOURCE_DIR}/include")

link_directories("${ECHMET_CORE_LIBS_DIR}/lib")

add_library(LEMNG SHARED ${libIonProps_SRCS})

set_target_properties(LEMNG
                      PROPERTIES VERSION 0.1
                                 SOVERSION 0.1
                      PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/lemng.h")
target_link_libraries(LEMNG
                      PRIVATE ECHMETShared
                      PRIVATE SysComp
                      PRIVATE CAES
                      PRIVATE IonProps
                      PRIVATE ${LEMNG_THREADING_LIBS})

if (WIN32)
    install(TARGETS LEMNG
            RUNTIME DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/LEMNG)
else ()
    install(TARGETS LEMNG
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/LEMNG)
endif ()
if (WIN32)
    if (MSVC)
	    set(BUILD_CFG $<CONFIG>)
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${BUILD_CFG}/LEMNG.lib" DESTINATION lib)
	else ()
	    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/LEMNG.lib" DESTINATION lib)
	endif ()
endif ()
