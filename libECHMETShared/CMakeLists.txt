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

if (NOT ENABLE_X86_EXTENSIONS)
    message(STATUS "x86 extensions disabled")
    set(ECHMETShared_SRCS
        ${ECHMETShared_SRCS}
        internal/cpufeatures_generic.cpp)
else ()
    set(ECHMETShared_SRCS
        ${ECHMETShared_SRCS}
        internal/cpufeatures_x86.cpp)
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

if ((BUILD_STATIC_LIBS) OR (EMSCRIPTEN))
    add_library(ECHMETShared STATIC ${ECHMETShared_SRCS})
else ()
    add_library(ECHMETShared SHARED ${ECHMETShared_SRCS})
endif ()

target_include_directories(ECHMETShared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
if (EMSCRIPTEN)
    set_target_properties(ECHMETShared
                          PROPERTIES VERSION 0.1)
else ()
    set_target_properties(ECHMETShared
                          PROPERTIES VERSION 0.1
                                     SOVERSION 0.1
                                     LINK_FLAGS ${DEFAULT_SYMVER_LINK})
endif ()
target_link_libraries(ECHMETShared
                      PRIVATE ${ECHMET_SHARED_LINK_LIBS}
                      PRIVATE ${ECHMET_HIPREC_LIBS})

install(TARGETS ECHMETShared
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ECHMET/CoreLibs)
install(FILES
            "${CMAKE_CURRENT_SOURCE_DIR}/echmetelems.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/echmetmath.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/echmetmodule.h"
            "${CMAKE_CURRENT_SOURCE_DIR}/echmetphchconsts.h"
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ECHMET/CoreLibs)
install(DIRECTORY containers/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ECHMET/CoreLibs/containers FILES_MATCHING PATTERN "*.h")
if (WIN32)
    if (MINGW OR MSYS)
        if (ECHMET_CORE_LIBS_XCOMPILE)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.dll.a" DESTINATION ${CMAKE_INSTALL_LIBDIR})
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.dll" DESTINATION ${CMAKE_INSTALL_BINDIR})
        else ()
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.lib" DESTINATION ${CMAKE_INSTALL_LIBDIR} RENAME ECHMETShared.lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libECHMETShared.dll" DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME ECHMETShared.dll)
        endif ()
    endif ()
endif ()
