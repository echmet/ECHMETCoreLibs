cmake_minimum_required(VERSION 2.8.12)

if (CMAKE_MAJOR_VERSION LESS 3)
    project(libCAES)
else ()
    project(libCAES LANGUAGES CXX)
endif ()

add_definitions("-D_USE_MATH_DEFINES -DEIGEN_MAX_ALIGN_BYTES=32")

if (WIN32)
    if (MINGW)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN_INCLUDE_DIR} -Wno-unused-function")
    elseif (MSYS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN_INCLUDE_DIR} -Wno-unused-function")
    else ()
        include_directories(${INCLUDE_DIRECTORIES}
                            "${EIGEN_INCLUDE_DIR}")
    endif ()
elseif (UNIX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN_INCLUDE_DIR} -Wno-unused-function")
    endif()
endif()

set(libCAES_SRCS
    src/caes.cpp
    src/caes_extended.cpp
    src/solvercontextimpl.cpp
    src/solverinternal_sse2.cpp
    src/solverinternal_avx.cpp
    src/solverinternal_fma3.cpp
    src/chargesummer_sse2.cpp
    src/chargesummer_avx.cpp
    src/chargesummer_fma3.cpp
    src/internals/newtonraphson_avx.cpp
    src/internals/newtonraphson_fma3.cpp
    src/vecmath/vecmath.cpp
    src/vecmath/vecmath_avx.cpp
    src/vecmath/vecmath_sse2.cpp
    src/vecmath/vecmath_fma3.cpp)

if (UNIX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        set_source_files_properties(src/chargesummer_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/chargesummer_fma3.cpp COMPILE_FLAGS "-mfma")
        set_source_files_properties(src/solverinternal_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/solverinternal_fma3.cpp COMPILE_FLAGS "-mfma")
        set_source_files_properties(src/vecmath/vecmath_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/vecmath/vecmath_fma3.cpp COMPILE_FLAGS "-mfma")
        set_source_files_properties(src/internals/newtonraphson_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/internals/newtonraphson_fma3.cpp COMPILE_FLAGS "-mfma")
    endif()
elseif (WIN32)
    if (MINGW OR MSYS)
        set_source_files_properties(src/chargesummer_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/chargesummer_fma3.cpp COMPILE_FLAGS "-mfma")
        set_source_files_properties(src/solverinternal_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/solverinternal_fma3.cpp COMPILE_FLAGS "-mfma")
        set_source_files_properties(src/vecmath/vecmath_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/vecmath/vecmath_fma3.cpp COMPILE_FLAGS "-mfma")
        set_source_files_properties(src/internals/newtonraphson_avx.cpp COMPILE_FLAGS "-mavx")
        set_source_files_properties(src/internals/newtonraphson_fma3.cpp COMPILE_FLAGS "-mfma")
    elseif (MSVC)
        set_source_files_properties(src/chargesummer_avx.cpp COMPILE_FLAGS "/arch:AVX")
        set_source_files_properties(src/chargesummer_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
        set_source_files_properties(src/solverinternal_avx.cpp COMPILE_FLAGS "/arch:AVX")
        set_source_files_properties(src/solverinternal_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
        set_source_files_properties(src/vecmath/vecmath_avx.cpp COMPILE_FLAGS "/arch:AVX")
        set_source_files_properties(src/vecmath/vecmath_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
        set_source_files_properties(src/internals/newtonraphson_avx.cpp COMPILE_FLAGS "/arch:AVX")
        set_source_files_properties(src/internals/newtonraphson_fma3.cpp COMPILE_FLAGS "/arch:AVX2")
	endif ()
endif()

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libECHMETShared/internal"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libIonProps/include")


set(CAES_LINK_LIBS
    ECHMETShared
    SysComp
    ${ECHMET_HIPREC_LIBS})
			
if (!WIN32)
    set(CAES_LINK_LIBS
        ${CAES_LINK_LIBS}
        m)
endif()

add_library(CAES SHARED ${libCAES_SRCS})
target_include_directories(CAES PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(CAES
                      PROPERTIES VERSION 0.2
                                 SOVERSION 0.2
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK})
target_link_libraries(CAES
                      PRIVATE ${CAES_LINK_LIBS})

if (WIN32)
    install(TARGETS CAES
            RUNTIME DESTINATION lib
            ARCHIVE DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
else ()
    install(TARGETS CAES
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
endif ()
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetcaes.h"
              "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetcaes_extended.h"
        DESTINATION include/ECHMET/CoreLibs)
if (WIN32)
    if (MINGW OR MSYS)
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libCAES.lib" DESTINATION lib RENAME CAES.lib)
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libCAES.dll" DESTINATION lib RENAME CAES.dll)
    endif ()
endif ()
