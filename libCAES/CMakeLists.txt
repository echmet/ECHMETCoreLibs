cmake_minimum_required(VERSION 2.8.12)

if (CMAKE_MAJOR_VERSION LESS 3)
    project(libCAES)
else ()
    project(libCAES LANGUAGES CXX)
endif ()

option(DISABLE_AVX512 "Disables use of AVX512 instruction set" OFF)

add_definitions("-D_USE_MATH_DEFINES -DEIGEN_MAX_ALIGN_BYTES=64")

if (DISABLE_AVX512)
    add_definitions("-DECHMET_DISABLE_AVX512")
endif ()

include_directories(${INCLUDE_DIRECTORIES}
                    SYSTEM "${EIGEN_INCLUDE_DIR}")

set(libCAES_SRCS
    src/caes.cpp
    src/caes_extended.cpp
    src/solvercontextimpl.cpp)

if (UNIX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
            set_source_files_properties(src/caes.cpp COMPILE_FLAGS "-Wno-unused-variable")
    endif ()
elseif (WIN32)
    if (MINGW OR MSYS)
        set_source_files_properties(src/caes.cpp COMPILE_FLAGS "-Wno-unused-variable")
    endif ()
endif ()

if (ENABLE_X86_EXTENSIONS)
    set(libCAES_SRCS
        ${libCAES_SRCS}
        src/solverinternal_sse2.cpp
        src/solverinternal_avx.cpp
        src/solverinternal_avx512.cpp
        src/solverinternal_fma3.cpp
        src/chargesummer_sse2.cpp
        src/chargesummer_avx.cpp
        src/chargesummer_avx512.cpp
        src/chargesummer_fma3.cpp
        src/estimator_helpers_sse2.cpp
        src/estimator_helpers_avx.cpp
        src/estimator_helpers_avx512.cpp
        src/estimator_helpers_fma3.cpp
        src/totalequilibrium_avx.cpp
        src/totalequilibrium_avx512.cpp
        src/totalequilibrium_fma3.cpp
        src/internals/newtonraphson_avx.cpp
        src/internals/newtonraphson_avx512.cpp
        src/internals/newtonraphson_fma3.cpp
        src/vecmath/vecmath.cpp
        src/vecmath/vecmath_avx.cpp
        src/vecmath/vecmath_avx512.cpp
        src/vecmath/vecmath_sse2.cpp
        src/vecmath/vecmath_fma3.cpp)

    if (UNIX)
        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
            ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
            set_source_files_properties(src/chargesummer_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/chargesummer_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/chargesummer_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/solverinternal_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/solverinternal_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/solverinternal_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/estimator_helpers_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/estimator_helpers_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/estimator_helpers_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/totalequilibrium_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/totalequilibrium_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/totalequilibrium_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/vecmath/vecmath_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/vecmath/vecmath_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/vecmath/vecmath_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/internals/newtonraphson_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/internals/newtonraphson_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/internals/newtonraphson_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
        endif()
    elseif (WIN32)
        if (MINGW OR MSYS)
            set_source_files_properties(src/chargesummer_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/chargesummer_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/chargesummer_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/solverinternal_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/solverinternal_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/solverinternal_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/estimator_helpers_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/estimator_helpers_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/estimator_helpers_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/totalequilibrium_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/totalequilibrium_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/totalequilibrium_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/vecmath/vecmath_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/vecmath/vecmath_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/vecmath/vecmath_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
            set_source_files_properties(src/internals/newtonraphson_avx.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx")
            set_source_files_properties(src/internals/newtonraphson_avx512.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma -mavx512f")
            set_source_files_properties(src/internals/newtonraphson_fma3.cpp COMPILE_FLAGS "-msse3 -mssse3 -msse4 -mavx -mavx2 -mfma")
        elseif (MSVC)
            set_source_files_properties(src/chargesummer_avx.cpp COMPILE_FLAGS "/arch:AVX")
            set_source_files_properties(src/chargesummer_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
            set_source_files_properties(src/solverinternal_avx.cpp COMPILE_FLAGS "/arch:AVX")
            set_source_files_properties(src/solverinternal_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
            set_source_files_properties(src/estimator_helpers_avx.cpp COMPILE_FLAGS "/arch:AVX")
            set_source_files_properties(src/estimator_helpers_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
            set_source_files_properties(src/totalequilibrium_avx.cpp COMPILE_FLAGS "/arch:AVX")
            set_source_files_properties(src/totalequilibrium_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
            set_source_files_properties(src/vecmath/vecmath_avx.cpp COMPILE_FLAGS "/arch:AVX")
            set_source_files_properties(src/vecmath/vecmath_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
            set_source_files_properties(src/internals/newtonraphson_avx.cpp COMPILE_FLAGS "/arch:AVX")
            set_source_files_properties(src/internals/newtonraphson_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
        endif ()
    endif ()
endif()

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libECHMETShared/internal"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libIonProps/include")

if (!WIN32)
    set(CAES_LINK_MATH
        PRIVATE m)
endif()

if ((BUILD_STATIC_LIBS) OR (EMSCRIPTEN))
    add_library(CAES STATIC ${libCAES_SRCS})
else ()
    add_library(CAES SHARED ${libCAES_SRCS})
endif ()

target_include_directories(CAES PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

if (EMSCRIPTEN)
    set_target_properties(ECHMETShared
                          PROPERTIES VERSION 0.1)
else ()
    set_target_properties(CAES
                          PROPERTIES VERSION 0.3
                                     SOVERSION 0.3
                                     LINK_FLAGS ${DEFAULT_SYMVER_LINK})
endif ()
target_link_libraries(CAES
                      PRIVATE ECHMETShared
                      PRIVATE SysComp
                      PRIVATE IonProps
                      PRIVATE ${ECHMET_HIPREC_LIBS}
                      ${CAES_LINK_MATH})

install(TARGETS CAES
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ECHMET/CoreLibs)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetcaes.h"
              "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetcaes_extended.h"
        DESTINATION include/ECHMET/CoreLibs)
