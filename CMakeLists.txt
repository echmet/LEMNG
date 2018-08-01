cmake_minimum_required(VERSION 2.8.12)

project(LEMNG LANGUAGES CXX)

option(SENSITIVE_NUMDERS "Use more precise input concentrations to calculate numerical derivatives at the expense of speed" OFF)
option(PARALLEL_NUM_OPS "Enable parallel numeric operations" ON)
option(DISABLE_TRACING "Disable tracing" OFF)

if (WIN32)
    option(ECHMET_FORCE_WINXP "Define target Windows version to Windows XP")
endif ()

add_definitions("-D_USE_MATH_DEFINES -DECHMET_DLL_BUILD")

set(LEMNG_THREADING_LIBS)

if (WIN32)
    if (MINGW)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
        add_definitions("-DECHMET_COMPILER_MINGW")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (MSYS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
        add_definitions("-DECHMET_COMPILER_MSYS")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (MSVC)
        add_definitions("-DECHMET_COMPILER_MSVC")
        set(DEFAULT_SYMVER_LINK " ")

        if (ECHMET_FORCE_WINXP)
            add_definitions("-D_WIN32_WINNT=0x0501")
        endif ()
    else ()
        message(FATAL "Unsupported or mistedected compiler")
    endif ()
elseif (UNIX)
    set(LEMNG_THREADING_LIBS pthread)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.2)
            message(FATAL_ERROR "GCC 4.9.2 or newer is required to build ECHMETCoreLibs")
        endif ()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -g")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
        message(FATAL_ERROR "Clang 3.4 or newer is required to build ECHMETCoreLibs")
	endif ()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -std=c++14 -Wall -Wextra -Wno-c++98-compat -Wno-unused-function -pedantic -fvisibility=hidden -g")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
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
                      PROPERTIES VERSION 0.5
                                 SOVERSION 0.5
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK}
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
