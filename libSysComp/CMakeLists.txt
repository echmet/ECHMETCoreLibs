cmake_minimum_required(VERSION 2.8.12)

project(libSysComp LANGUAGES CXX)

include_directories(${INCLUDE_DIRECTORIES}
                    "${CMAKE_CURRENT_SOURCE_DIR}/include")

set(libSysComp_SRCS
    src/syscomp.cpp
   )

set(SYSCOMP_LINK_LIBS
    ECHMETShared)
			
if (!WIN32)
    set(SYSCOMP_LINK_LIBS
        ${SYSCOMP_LINK_LIBS}
        m)
endif ()

if (ECHMET_USE_HIGH_PRECISION)
    set(SYSCOMP_LINK_LIBS
        ${SYSCOMP_LINK_LIBS}
	${ECHMET_HIPREC_LIBS})
endif ()

add_library(SysComp SHARED ${libSysComp_SRCS})
target_include_directories(SysComp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(SysComp
                      PROPERTIES VERSION 0.1
                                 SOVERSION 0.1
                      PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/include/echmetsyscomp.h")
target_link_libraries(SysComp
                      PRIVATE ${SYSCOMP_LINK_LIBS})

if (WIN32)
    install(TARGETS SysComp
            RUNTIME DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
else ()
    install(TARGETS SysComp
            LIBRARY DESTINATION lib
            PUBLIC_HEADER DESTINATION include/ECHMET/CoreLibs)
endif ()
