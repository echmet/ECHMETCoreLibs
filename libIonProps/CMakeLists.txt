cmake_minimum_required(VERSION 2.8.12)

if (CMAKE_MAJOR_VERSION LESS 3)
    project(libIonProps)
else ()
    project(libIonProps LANGUAGES CXX)
endif ()

option(IONPROPS_DISABLE_COMPLEX_ONSFUO "Do not correct mobilities of complexes by Onsager-Fuoss law" ON)

add_definitions("-D_USE_MATH_DEFINES")

if (IONPROPS_DISABLE_COMPLEX_ONSFUO)
    add_definitions("-DIONPROPS_DISABLE_COMPLEX_ONSFUO")
endif ()

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
    endif ()
endif ()

set(libIonProps_SRCS
    src/ionprops.cpp)

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libSysComp/include")

set(IONPROPS_LINK_LIBS
    ECHMETShared
    SysComp)
			
if (!WIN32)
    set(IONPROPS_LINK_LIBS
        ${IONPROPS_LINK_LIBS}
        m)
endif ()

set(IONPROPS_LINK_LIBS
    "${IONPROPS_LINK_LIBS}"
    ${ECHMET_HIPREC_LIBS})

if (BUILD_STATIC_LIBS)
    add_library(IonProps STATIC ${libIonProps_SRCS})
else ()
    add_library(IonProps SHARED ${libIonProps_SRCS})
endif ()

target_include_directories(IonProps PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(IonProps
                      PROPERTIES VERSION 0.2
                                 SOVERSION 0.2
                                 LINK_FLAGS ${DEFAULT_SYMVER_LINK}
                      PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetionprops.h")
target_link_libraries(IonProps
                      PRIVATE ${IONPROPS_LINK_LIBS})

if (WIN32)
    install(TARGETS IonProps
            RUNTIME DESTINATION lib
            ARCHIVE DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
else ()
    install(TARGETS IonProps
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
endif ()
if (WIN32)
    if (MINGW OR MSYS)
        if (ECHMET_CORE_LIBS_XCOMPILE)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libIonProps.dll.a" DESTINATION lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libIonProps.dll" DESTINATION lib)
        else ()
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libIonProps.lib" DESTINATION lib RENAME IonProps.lib)
            install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libIonProps.dll" DESTINATION lib RENAME IonProps.dll)
        endif ()
    endif ()
endif ()
