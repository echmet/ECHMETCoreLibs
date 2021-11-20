cmake_minimum_required(VERSION 2.8.12)

if (CMAKE_MAJOR_VERSION LESS 3)
    project(libSysComp)
else ()
    project(libSysComp LANGUAGES CXX)
endif ()

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(libSysComp_SRCS
    src/syscomp.cpp
   )

if (!WIN32)
    set(SYSCOMP_LINK_MATH
        "PRIVATE m")
endif ()

if (ECHMET_USE_HIGH_PRECISION)
    set(SYSCOMP_LINK_LIBS
        ${SYSCOMP_LINK_LIBS}
	${ECHMET_HIPREC_LIBS})
endif ()

if (BUILD_STATIC_LIBS)
    add_library(SysComp STATIC ${libSysComp_SRCS})
else ()
    add_library(SysComp SHARED ${libSysComp_SRCS})
endif ()

target_include_directories(SysComp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(SysComp
                      PROPERTIES VERSION 0.1
                                 SOVERSION 0.1
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK}
                      PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetsyscomp.h")
target_link_libraries(SysComp
                      PRIVATE ECHMETShared
                      ${SYSCOMP_LINK_MATH})

if (WIN32)
    install(TARGETS SysComp
            RUNTIME DESTINATION lib
            ARCHIVE DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
else ()
    install(TARGETS SysComp
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
endif ()
if (WIN32)
    if (MINGW OR MSYS)
        if (ECHMET_CORE_LIBS_XCOMPILE)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libSysComp.dll.a" DESTINATION lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libSysComp.dll" DESTINATION lib)
        else ()
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libSysComp.lib" DESTINATION lib RENAME SysComp.lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libSysComp.dll" DESTINATION lib RENAME SysComp.dll)
        endif ()
    endif ()
endif ()
