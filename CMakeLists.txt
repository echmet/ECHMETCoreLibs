cmake_minimum_required(VERSION 2.8.12)

project(ECHMET_CORE_LIBS)

option(ECHMET_USE_HIGH_PRECISION "Use arbitrary precision data type instead of IEEE754 double everywhere" OFF)
option(ECHMET_HACKS_ENABLED "Enable ad-hoc hacks - do not switch this on" OFF)
option(MANUAL_HIPREC_LIBS_PATH "Use manually specified paths to GMP and MPFR libraries on UNIX" OFF)

if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions("-DECHMET_DEBUG_OUTPUT")
else()
    add_definitions("-DEIGEN_NO_DEBUG")
endif()

if (ECHMET_USE_HIGH_PRECISION)
    set("ECHMET_USE_HIGH_PRECISION" "1")
endif()

if (ECHMET_HACKS_ENABLED)
    add_definitions("-DECHMET_HACKS_MD")
endif()

add_definitions("-DECHMET_DLL_BUILD -DUSE_ECHMET_CONTAINERS")

if (WIN32)
    set("ECHMET_PLATFORM_WIN32" "1")
	
    set(ECHMET_HIPREC_LIBS
        ${GMP_LIBRARY_BIN}
        ${MPFR_LIBRARY_BIN})
    set(ECHMET_HIPREC_INCLUDE_DIR
        ${GMP_INCLUDE_DIR}
        ${MPFR_INCLUDE_DIR})
		
    include_directories(SYSTEM ${ECHMET_HIPREC_INCLUDE_DIR})

    if (MINGW)
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.2)
            message(FATAL_ERROR "MinGW 4.9.2 or newer is required to build ECHMETCoreLibs")
        endif ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -ftree-vectorize -g -O3 -msse2")
        add_definitions("-DECHMET_COMPILER_MINGW")
    elseif (MSYS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -ftree-vectorize -g -O3 -msse2")
        add_definitions("-DECHMET_COMPILER_MSYS")
    elseif (MSVC)
        add_definitions("-DECHMET_COMPILER_MSVC")
        set(DEFAULT_SYMVER_LINK " ")
    else()
        message(FATAL_ERROR "Compiler not recognized or not supported")
    endif()
elseif (UNIX)
    set("ECHMET_PLATFORM_UNIX" "1")

    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND
            CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.2)
            message(FATAL_ERROR "GCC 4.9.2 or newer is required to build ECHMETCoreLibs")
        elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" AND
	        CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.4)
            message(FATAL_ERROR "Clang 3.4 or newer is required to build ECHMETCoreLibs")
        endif ()

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -g -O3")
        add_definitions("-DECHMET_COMPILER_GCC_LIKE ")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    endif ()

    if (MANUAL_HIPREC_LIBS_PATH)
        set(ECHMET_HIPREC_LIBS
            ${GMP_LIBRARY_BIN}
            ${MPFR_LIBRARY_BIN})
        set(ECHMET_HIPREC_INCLUDE_DIR
            ${GMP_INCLUDE_DIR}
            ${MPFR_INCLUDE_DIR})
    
        include_directories(SYSTEM ${ECHMET_HIPREC_INCLUDE_DIR})
    else ()
        set(ECHMET_HIPREC_LIBS
            gmp
            mpfr)
    endif ()
else()
    message(FATAL_ERROR "Platform not recognized or not supported")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/libECHMETShared/echmetcorelibs_config.h.in" "${CMAKE_CURRENT_BINARY_DIR}/echmetcorelibs_config.h")

include_directories(${INCLUDE_DIRECTORIES}
                    ${CMAKE_CURRENT_BINARY_DIR}
                   "${CMAKE_CURRENT_SOURCE_DIR}/libECHMETShared")

add_subdirectory(libECHMETShared)
add_subdirectory(libSysComp)
add_subdirectory(libIonProps)
add_subdirectory(libCAES)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/echmetcorelibs_config.h"
        DESTINATION include/ECHMET/CoreLibs)
