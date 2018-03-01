cmake_minimum_required(VERSION 2.8.12)

project(libCAES LANGUAGES CXX)

add_definitions("-D_USE_MATH_DEFINES")

if (WIN32)
    if (MINGW)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN_INCLUDE_DIR} -Wno-unused-function")
    elseif (MSYS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN_INCLUDE_DIR} -Wno-unused-function")
    endif()
elseif (UNIX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
	${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${EIGEN_INCLUDE_DIR} -Wno-unused-function")
    endif()
endif()

set(libCAES_SRCS
    src/caes.cpp
    src/caes_extended.cpp)

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../libECHMETShared/internal"
		    "${CMAKE_CURRENT_SOURCE_DIR}/../libIonProps/include"
		    "${EIGEN_INCLUDE_DIR}")

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
                                 SOVERSION 0.2)
target_link_libraries(CAES
                      PRIVATE ${CAES_LINK_LIBS})

if (WIN32)
    install(TARGETS CAES
            RUNTIME DESTINATION lib
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
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libCAES.lib" DESTINATION lib RENAME CAES.lib)
	install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libCAES.dll" DESTINATION lib RENAME CAES.dll)
endif ()
