cmake_minimum_required(VERSION 2.8.12)

project(LEMNG LANGUAGES CXX)

option(SENSITIVE_NUMDERS "Use more precise input concentrations to calculate numerical derivatives at the expense of speed" OFF)
option(PARALLEL_NUM_OPS "Enable parallel numeric operations" ON)
option(DISABLE_TRACING "Disable tracing" OFF)
option(BUILD_TESTS "Build unit tests" ON)

if (WIN32)
    option(ECHMET_FORCE_WINXP "Define target Windows version to Windows XP" OFF)
endif ()

if (UNIX)
    option(USE_ADDRESS_SANITIZER "Use address sanitizer. This requires sufficiently current version of GCC or Clang" OFF)
    option(USE_THREAD_SANITIZER "Use thread sanitizer. This requires sufficiently current version of GCC or Clang" OFF)
endif ()

if ((USE_THREAD_SANITIZER) AND (USE_ADDRESS_SANITIZER))
    message(FATAL_ERROR "Address and thread sanitizer cannot be used together")
endif ()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_definitions("-D_USE_MATH_DEFINES -DECHMET_DLL_BUILD")

if (UNIX)
    include(GNUInstallDirs)
elseif (WIN32)
    set(CMAKE_INSTALL_LIBDIR "lib")
    set(CMAKE_INSTALL_INCLUDEDIR "include")
    set(CMAKE_INSTALL_BINDIR "bin")
endif ()

set(LEMNG_THREADING_LIBS)

set(LEMNG_VERSION_MAJOR 0)
set(LEMNG_VERSION_MINOR 6)
set(LEMNG_VERSION_PATCH 2)

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

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -g")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")

        if (USE_ADDRESS_SANITIZER)
            add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize=leak -fstack-protector-strong)
            add_link_options(-fsanitize=address -fsanitize=undefined -fsanitize=leak)
        endif ()
        if (USE_THREAD_SANITIZER)
            add_compile_options(-fno-omit-frame-pointer -fsanitize=thread -fstack-protector-strong)
            add_link_options(-fsanitize=thread)
        endif ()
    elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
            message(FATAL_ERROR "Clang 3.4 or newer is required to build ECHMETCoreLibs")
        endif ()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -Wall -Wextra -Wno-c++98-compat -Wno-unused-function -pedantic -fvisibility=hidden -g")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")

        if (USE_ADDRESS_SANITIZER)
            add_compile_options(-fno-omit-frame-pointer -fsanitize=address -fsanitize=undefined -fsanitize=leak -fstack-protector-strong)
            add_link_options(-fsanitize=address -fsanitize=undefined -fsanitize=leak)
        endif ()
        if (USE_THREAD_SANITIZER)
            add_compile_options(-fno-omit-frame-pointer -fsanitize=thread -fstack-protector-strong)
            add_link_options(-fsanitize=thread)
        endif ()
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
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_BINARY_DIR}")

if (NOT "${EIGEN_INCLUDE_DIR}" STREQUAL "")
    include_directories(${INCLUDE_DIRECTORIES}
                        SYSTEM "${EIGEN_INCLUDE_DIR}")
endif ()

if (NOT "${ECHMET_CORE_LIBS_DIR}" STREQUAL "")
    include_directories(${INCLUDE_DIRECTORIES}
                        SYSTEM "${ECHMET_CORE_LIBS_DIR}/include/ECHMET/CoreLibs")
    link_directories(${LINK_DIRECTORIES}
                     "${ECHMET_CORE_LIBS_DIR}/${CMAKE_INSTALL_LIBDIR}")
endif ()

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/src/lemng_config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/lemng_config.h")

add_library(LEMNG SHARED ${libIonProps_SRCS})

set_target_properties(LEMNG
                      PROPERTIES VERSION ${LEMNG_VERSION_MAJOR}.${LEMNG_VERSION_MINOR}
                                 SOVERSION ${LEMNG_VERSION_MAJOR}.${LEMNG_VERSION_MINOR}
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK}
                      PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/lemng.h")
target_link_libraries(LEMNG
                      PRIVATE ECHMETShared
                      PRIVATE SysComp
                      PRIVATE CAES
                      PRIVATE IonProps
                      PRIVATE ${LEMNG_THREADING_LIBS})

if (BUILD_TESTS)
    enable_testing()

    add_executable(benzoiclibcd_empty_nois_exe src/tests/benzoiclibcd_empty_nois.cpp)
    target_link_libraries(benzoiclibcd_empty_nois_exe PRIVATE LEMNG
                                                      PRIVATE ECHMETShared
                                                      PRIVATE SysComp)
    add_test(benzoiclibcd_empty_nois benzoiclibcd_empty_nois_exe)

    add_executable(benzoiclibcd_empty_is_exe src/tests/benzoiclibcd_empty_is.cpp)
    target_link_libraries(benzoiclibcd_empty_is_exe PRIVATE LEMNG
                                                    PRIVATE ECHMETShared
                                                    PRIVATE SysComp)
    add_test(benzoiclibcd_empty_is benzoiclibcd_empty_is_exe)

    add_executable(formlibcddmbcd_rflu_nois_exe src/tests/formlibcddmbcd_rflu_nois.cpp)
    target_link_libraries(formlibcddmbcd_rflu_nois_exe PRIVATE LEMNG
                                                       PRIVATE ECHMETShared
                                                       PRIVATE SysComp)
    add_test(formlibcddmbcd_rflu_nois formlibcddmbcd_rflu_nois_exe)

    add_executable(formlibcddmbcd_rflu_is_exe src/tests/formlibcddmbcd_rflu_is.cpp)
    target_link_libraries(formlibcddmbcd_rflu_is_exe PRIVATE LEMNG
                                                     PRIVATE ECHMETShared
                                                     PRIVATE SysComp)
    add_test(formlibcddmbcd_rflu_is formlibcddmbcd_rflu_is_exe)

    add_executable(formlibcd_rflu_nois_exe src/tests/formlibcd_rflu_nois.cpp)
    target_link_libraries(formlibcd_rflu_nois_exe PRIVATE LEMNG
                                                  PRIVATE ECHMETShared
                                                  PRIVATE SysComp)
    add_test(formlibcd_rflu_nois formlibcd_rflu_nois_exe)

    add_executable(formlibcd_rflu_is_exe src/tests/formlibcd_rflu_is.cpp)
    target_link_libraries(formlibcd_rflu_is_exe PRIVATE LEMNG
                                                PRIVATE ECHMETShared
                                                PRIVATE SysComp)
    add_test(formlibcd_rflu_is formlibcd_rflu_is_exe)

    add_executable(formna_k_nois_exe src/tests/formna_k_nois.cpp)
    target_link_libraries(formna_k_nois_exe PRIVATE LEMNG
                                            PRIVATE ECHMETShared
                                            PRIVATE SysComp)
    add_test(formna_k_nois formna_k_nois_exe)

    add_executable(formna_k_is_exe src/tests/formna_k_is.cpp)
    target_link_libraries(formna_k_is_exe PRIVATE LEMNG
                                          PRIVATE ECHMETShared
                                          PRIVATE SysComp)
    add_test(formna_k_is formna_k_is_exe)

    add_executable(nacl_empty_nois_exe src/tests/nacl_empty_nois.cpp)
    target_link_libraries(nacl_empty_nois_exe PRIVATE LEMNG
                                              PRIVATE ECHMETShared
                                              PRIVATE SysComp)
    add_test(nacl_empty_nois nacl_empty_nois_exe)

    add_executable(nacl_empty_is_exe src/tests/nacl_empty_is.cpp)
    target_link_libraries(nacl_empty_is_exe PRIVATE LEMNG
                                            PRIVATE ECHMETShared
                                            PRIVATE SysComp)
    add_test(nacl_empty_is nacl_empty_is_exe)

    add_executable(oscillating_nois_exe src/tests/oscillating_nois.cpp)
    target_link_libraries(oscillating_nois_exe PRIVATE LEMNG
                                               PRIVATE ECHMETShared
                                               PRIVATE SysComp)
    add_test(oscillating_nois oscillating_nois_exe)

    add_executable(oscillating_is_exe src/tests/oscillating_is.cpp)
    target_link_libraries(oscillating_is_exe PRIVATE LEMNG
                                             PRIVATE ECHMETShared
                                             PRIVATE SysComp)
    add_test(oscillating_is oscillating_is_exe)

    add_executable(trtrbcd_snap_sibu_nois_exe src/tests/trtrbcd_snap_sibu_nois.cpp)
    target_link_libraries(trtrbcd_snap_sibu_nois_exe PRIVATE LEMNG
                                                     PRIVATE ECHMETShared
                                                     PRIVATE SysComp)
    add_test(trtrbcd_snap_sibu_nois trtrbcd_snap_sibu_nois_exe)

    add_executable(trtrbcd_snap_sibu_is_exe src/tests/trtrbcd_snap_sibu_is.cpp)
    target_link_libraries(trtrbcd_snap_sibu_is_exe PRIVATE LEMNG
                                                   PRIVATE ECHMETShared
                                                   PRIVATE SysComp)
    add_test(trtrbcd_snap_sibu_is trtrbcd_snap_sibu_is_exe)

    add_executable(phosphorical_empty_nois_exe src/tests/phosphorical_empty_nois.cpp)
    target_link_libraries(phosphorical_empty_nois_exe PRIVATE LEMNG
                                                      PRIVATE ECHMETShared
                                                      PRIVATE SysComp)
    add_test(phosphorical_empty_nois phosphorical_empty_nois_exe)

    add_executable(phosphorical_empty_is_exe src/tests/phosphorical_empty_is.cpp)
    target_link_libraries(phosphorical_empty_is_exe PRIVATE LEMNG
                                                    PRIVATE ECHMETShared
                                                    PRIVATE SysComp)
    add_test(phosphorical_empty_is phosphorical_empty_is_exe)

    add_executable(formlixs_analyte_nois_exe src/tests/formlixs_analyte_nois.cpp)
    target_link_libraries(formlixs_analyte_nois_exe PRIVATE LEMNG
                                                    PRIVATE ECHMETShared
                                                    PRIVATE SysComp)
    add_test(formlixs_analyte_nois formlixs_analyte_nois_exe)

    add_executable(formlixs_analyte_is_exe src/tests/formlixs_analyte_is.cpp)
    target_link_libraries(formlixs_analyte_is_exe PRIVATE LEMNG
                                                  PRIVATE ECHMETShared
                                                  PRIVATE SysComp)
    add_test(formlixs_analyte_is formlixs_analyte_is_exe)

    add_executable(formlixs_analyte_inv_nois_exe src/tests/formlixs_analyte_inv_nois.cpp)
    target_link_libraries(formlixs_analyte_inv_nois_exe PRIVATE LEMNG
                                                        PRIVATE ECHMETShared
                                                        PRIVATE SysComp)
    add_test(formlixs_analyte_inv_nois formlixs_analyte_inv_nois_exe)

    add_executable(formlixs_analyte_inv_is_exe src/tests/formlixs_analyte_inv_is.cpp)
    target_link_libraries(formlixs_analyte_inv_is_exe PRIVATE LEMNG
                                                      PRIVATE ECHMETShared
                                                      PRIVATE SysComp)
    add_test(formlixs_analyte_inv_is formlixs_analyte_inv_is_exe)

    add_executable(formlixs_analyte_sys_nois_exe src/tests/formlixs_analyte_sys_nois.cpp)
    target_link_libraries(formlixs_analyte_sys_nois_exe PRIVATE LEMNG
                                                        PRIVATE ECHMETShared
                                                        PRIVATE SysComp)
    add_test(formlixs_analyte_sys_nois formlixs_analyte_sys_nois_exe)

    add_executable(formlixs_analyte_sys_is_exe src/tests/formlixs_analyte_sys_is.cpp)
    target_link_libraries(formlixs_analyte_sys_is_exe PRIVATE LEMNG
                                                      PRIVATE ECHMETShared
                                                      PRIVATE SysComp)
    add_test(formlixs_analyte_sys_is formlixs_analyte_sys_is_exe)
endif()

install(TARGETS LEMNG
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ECHMET/LEMNG)
