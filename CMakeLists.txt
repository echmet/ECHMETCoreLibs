cmake_minimum_required(VERSION 2.8.12)

project(ECHMET_CORE_LIBS)

option(ECHMET_USE_HIGH_PRECISION "Use arbitrary precision data type instead of IEEE754 double everywhere" OFF)
option(ECHMET_HACKS_ENABLED "Enable ad-hoc hacks - do not switch this on" OFF)
option(MANUAL_HIPREC_LIBS_PATH "Use manually specified paths to GMP and MPFR libraries on UNIX" OFF)
option(ENABLE_DIAGNOSTIC_OUTPUT "Enables diagnostic text output. This output is enabled automatically in debug builds" OFF)
option(USE_CLANG_FPROF "Enables Clang-specific code profiling instrumentation" OFF)
option(BUILD_STATIC_LIBS "Build static libraries" OFF)
option(DISABLE_HANDCRAFTED_SIMD "Disable hand-crafted SIMD code" OFF)

if (WIN32)
    option(ECHMET_FORCE_WINXP "Define target Windows version to Windows XP")
endif ()

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

if (ENABLE_DIAGNOSTIC_OUTPUT)
    add_definitions("-DECHMET_DEBUG_OUTPUT")
endif()

add_definitions("-DECHMET_DLL_BUILD -DUSE_ECHMET_CONTAINERS")

if (WIN32)
    set("ECHMET_PLATFORM_WIN32" "1")
    if (ECHMET_FORCE_WINXP)
        add_definitions("-D_WIN32_WINNT=0x0501")
    endif ()

    link_directories(${LINK_DIRECTORIES}
                     "${LIBGMP_DIR}/bin/"
                     "${LIBMPFR_DIR}/bin/")

    set(ECHMET_HIPREC_INCLUDE_DIR
        "${LIBGMP_DIR}/include"
        "${LIBMPFR_DIR}/include")

    include_directories(SYSTEM ${ECHMET_HIPREC_INCLUDE_DIR})

    if (MSVC)
        set(ECHMET_HIPREC_LIBS
            libgmp-10
            libmpfr-6)
    elseif (ECHMET_CORE_LIBS_XCOMPILE)
        set(ECHMET_HIPREC_LIBS
            libgmp
            libmpfr)
    else ()
        set(ECHMET_HIPREC_LIBS
            gmp-10
            mpfr-6)
    endif ()

    if (MINGW)
        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9.2)
            message(FATAL_ERROR "MinGW 4.9.2 or newer is required to build ECHMETCoreLibs")
        endif ()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -ftree-vectorize -g -O3 -msse2")
        add_definitions("-DECHMET_COMPILER_MINGW")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (MSYS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -ftree-vectorize -g -O3 -msse2")
        add_definitions("-DECHMET_COMPILER_MSYS")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    elseif (MSVC)
        set(CMAKE_CONFIGURATION_TYPES "Debug;Release;RelWithDebInfo;MinSizeRel;RelDebOpti")
        set(CMAKE_CXX_FLAGS_RELDEBOPTI "/MD /Zi /Ox /Ob2 /Oy /GL /Oi /DNDEBUG")
        set(CMAKE_SHARED_LINKER_FLAGS_RELDEBOPTI "/debug /INCREMENTAL:NO /LTCG")
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

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -std=c++14 -Wall -Wextra -pedantic -Wno-unused-function -fvisibility=hidden -ftree-vectorize -g -O3 -msse2")

        if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" AND
            CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 9.1.0 OR
            CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL 9.1.0)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
        endif ()
	if (USE_CLANG_FPROF)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
		set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fprofile-instr-generate")
	endif()

        add_definitions("-DECHMET_COMPILER_GCC_LIKE ")
        set(DEFAULT_SYMVER_LINK "-Wl,--default-symver")
    endif ()

    if (MANUAL_HIPREC_LIBS_PATH)
	    link_directories(${LINK_DIRECTORIES}
	                     "${LIBGMP_DIR}/lib/"
                         "${LIBMPFR_DIR}/lib/")
        set(ECHMET_HIPREC_INCLUDE_DIR
            "${LIBGMP_DIR}/include"
            "${LIBMPFR_DIR}/include")

        include_directories(SYSTEM ${ECHMET_HIPREC_INCLUDE_DIR})
    endif ()

    set(ECHMET_HIPREC_LIBS
        gmp
        mpfr)
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
