#! /bin/sh
source ./ref_tool_glob.sh

clang -c jsonloader/constituents_json_ldr.c \
	-I${LIBJANSSON_INCLUDE}

clang++ -std=c++11 -Wall -Wextra -pedantic -g -O0 \
	-c json_input_processor.cpp jsonloader/inputreader.cpp \
	-I${ECL_INCLUDE} \
	-DECHMET_COMPILER_GCC_LIKE

clang++ -std=c++98 -Wall -Wextra -pedantic -g -O0 \
	ref_tool.cpp \
	json_input_processor.o inputreader.o \
	constituents_json_ldr.o /home/madcat/Devel/ECHMET/jansson-bin/lib/libjansson.a \
	-o ref_tool \
	-I../include \
	-I${ECL_INCLUDE} \
	-DECHMET_COMPILER_GCC_LIKE \
	-L${ECL_BIN} \
	-L../build \
	-lECHMETShared -lSysComp -lCAES -lIonProps
	
