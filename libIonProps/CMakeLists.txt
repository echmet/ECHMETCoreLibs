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

include_directories(${INCLUDE_DIRECTORIES}
                    SYSTEM "${EIGEN_INCLUDE_DIR}")

set(libIonProps_SRCS
    src/ionprops.cpp)

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libSysComp/include")

if (!WIN32)
    set(IONPROPS_LINK_MATH
        ${IONPROPS_LINK_MATH}
        PRIVATE m)
endif ()


if ((BUILD_STATIC_LIBS) OR (EMSCRIPTEN))
    add_library(IonProps STATIC ${libIonProps_SRCS})
else ()
    add_library(IonProps SHARED ${libIonProps_SRCS})
endif ()

target_include_directories(IonProps PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
if (EMSCRIPTEN)
    set_target_properties(IonProps
                          PROPERTIES VERSION 0.2
                                     SOVERSION 0.2
                          PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetionprops.h")
else ()
    set_target_properties(IonProps
                          PROPERTIES VERSION 0.2
                                     SOVERSION 0.2
                                     LINK_FLAGS ${DEFAULT_SYMVER_LINK}
                          PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetionprops.h")
endif ()
target_link_libraries(IonProps
                      PRIVATE ECHMETShared
                      PRIVATE SysComp
                      PRIVATE ${ECHMET_HIPREC_LIBS}
                      ${IONPROPS_LINK_MATH})

install(TARGETS IonProps
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ECHMET/CoreLibs)
