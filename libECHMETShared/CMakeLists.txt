cmake_minimum_required(VERSION 2.8.12)

if (CMAKE_MAJOR_VERSION LESS 3)
    project(libECHMETShared)
else ()
    project(libECHMETShared LANGUAGES CXX)
endif ()

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}")

set(ECHMETShared_SRCS
    internal/echmetelems.cpp
    internal/echmetmath.cpp
    internal/initializers.cpp
    internal/phchconsts.cpp
    internal/cpufeatures_common.cpp)

if (DISABLE_HANDCRAFTED_SIMD)
    set(ECHMETShared_SRCS
        ${ECHMETShared_SRCS}
        internal/cpufeatures_generic.cpp)
else ()
    if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" OR
        ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64" OR WIN32)
        set(ECHMETShared_SRCS
            ${ECHMETShared_SRCS}
            internal/cpufeatures_x86.cpp)
    else ()
         set(ECHMETShared_SRCS
            ${ECHMETShared_SRCS}
            internal/cpufeatures_generic.cpp)
    endif ()
endif()

if (!WIN32)
    set(ECHMET_SHARED_LINK_LIBS
        ${ECHMET_SHARED_LINK_LIBS}
        m)
endif()

if (ECHMET_USE_HIGH_PRECISION)
    set(ECHMET_SHARED_LINK_LIBS
        ${ECHMET_SHARED_LINK_LIBS}
        ${ECHMET_HIPREC_LIBS})
endif ()

if (BUILD_STATIC_LIBS)
    add_library(ECHMETShared STATIC ${ECHMETShared_SRCS})
else ()
    add_library(ECHMETShared SHARED ${ECHMETShared_SRCS})
endif ()

target_include_directories(ECHMETShared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(ECHMETShared
                      PROPERTIES VERSION 0.1
                                 SOVERSION 0.1
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK})
target_link_libraries(ECHMETShared
                      PRIVATE ${ECHMET_SHARED_LINK_LIBS}
                      PRIVATE ${ECHMET_HIPREC_LIBS})

if (WIN32)
    install(TARGETS ECHMETShared
            RUNTIME DESTINATION lib
            ARCHIVE DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
else ()
    install(TARGETS ECHMETShared
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
endif ()
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/echmetelems.h"
              "${CMAKE_CURRENT_SOURCE_DIR}/echmetmath.h"
              "${CMAKE_CURRENT_SOURCE_DIR}/echmetmodule.h"
              "${CMAKE_CURRENT_SOURCE_DIR}/echmetphchconsts.h"
        DESTINATION include/ECHMET/CoreLibs)
install(DIRECTORY containers/ DESTINATION include/ECHMET/CoreLibs/containers FILES_MATCHING PATTERN "*.h")
if (WIN32)
    if (MINGW OR MSYS)
        if (ECHMET_CORE_LIBS_XCOMPILE)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.dll.a" DESTINATION lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.dll" DESTINATION lib)
        else ()
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.lib" DESTINATION lib RENAME ECHMETShared.lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.dll" DESTINATION lib RENAME ECHMETShared.dll)
        endif ()
    endif ()
endif ()

