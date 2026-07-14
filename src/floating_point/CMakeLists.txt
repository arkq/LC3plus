#******************************************************************************
#                        ETSI TS 103 634 V1.7.1                               *
#              Low Complexity Communication Codec Plus (LC3plus)              *
#                                                                             *
# Copyright licence is solely granted through ETSI Intellectual Property      *
# Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
# estoppel or otherwise.                                                      *
#*****************************************************************************/

cmake_minimum_required(VERSION 3.13)

if(WINDOWS)
    set(CMAKE_SYSTEM_NAME Windows)
    set(CMAKE_C_COMPILER i686-w64-mingw32-gcc)
endif()

project(LC3plus_FL C)

# ------------------ Options ------------------
set(NAME_LC3 "LC3plus" CACHE STRING "Output binary name")
set(LIB_LC3  "LC3plus" CACHE STRING "Shared library base name (-> libLC3plus.so)")

option(NO_POST_REL_CHANGES_TEST "Disable post-release non-bitexact changes" OFF)
option(SHORT_PLC_FADEOUT        "Enable short PLC fadeout tuning" OFF)
option(GCOV                     "Build with coverage instrumentation" OFF)
option(AFL                      "Build for AFL fuzzing" OFF)

set(OPTIM "0" CACHE STRING "Optimization level 0..3")
set(CLANG "0" CACHE STRING "0=gcc, 1=clang+msan, 2=clang+asan, 3=clang+ubsan")
set(WMCFR "" CACHE STRING "Enable WMC instrumentation with given frame size (5/10/...). Empty = off.")

if(CLANG STREQUAL "1" OR CLANG STREQUAL "2" OR CLANG STREQUAL "3")
    set(CMAKE_C_COMPILER clang)
    set(OPTIM "2")
endif()

# ------------------ Compile flags ------------------
set(WARN_FLAGS
    -pedantic -Wcast-qual -Wall -W -Wextra -Wno-long-long
    -Wpointer-arith -Wstrict-prototypes -Wmissing-prototypes
    -Werror-implicit-function-declaration)

if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    add_compile_options(-std=c99 -fPIC ${WARN_FLAGS})
endif()

add_link_options(-g)

if(NOT DEFINED DEBUG OR NOT DEBUG STREQUAL "0")
    add_compile_options(-g3)
    add_link_options(-g3)
endif()

if(SHORT_PLC_FADEOUT)
    add_compile_definitions(PLC_TUNING_SHORT_FADEOUT)
endif()

if(NO_POST_REL_CHANGES_TEST)
    add_compile_definitions(NO_POST_REL_CHANGES)
endif()

if(WMCFR)
    add_compile_definitions(WMOPS WMC_FRAME_SIZE_MS=${WMCFR})
endif()

if(GCOV)
    add_compile_options(-fprofile-arcs -ftest-coverage)
    add_link_options(-fprofile-arcs -ftest-coverage)
endif()

if(CLANG STREQUAL "1")
    add_compile_options(-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer)
    add_link_options(-fsanitize=memory -fsanitize-memory-track-origins -fno-omit-frame-pointer)
elseif(CLANG STREQUAL "2")
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address -fno-omit-frame-pointer)
elseif(CLANG STREQUAL "3")
    add_compile_options(-fsanitize=undefined)
    add_link_options(-fsanitize=undefined)
endif()

add_compile_options(-O${OPTIM})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# ------------------ Source lists ------------------
# fft/*.c is #included from lc3plus_fft.c, not compiled standalone.
file(GLOB ALL_SRCS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.c")

set(LC3PLUS_SRCS ${ALL_SRCS})

# Shared lib excludes the two driver entry points
set(LIB_SRCS ${ALL_SRCS})
list(FILTER LIB_SRCS EXCLUDE REGEX "/ccConvert\\.c$")
list(FILTER LIB_SRCS EXCLUDE REGEX "/codec_exe\\.c$")

# ------------------ Targets ------------------
add_executable(${NAME_LC3} ${LC3PLUS_SRCS})
target_link_libraries(${NAME_LC3} PRIVATE m)
set_target_properties(${NAME_LC3} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_library(LC3plus_shared SHARED EXCLUDE_FROM_ALL ${LIB_SRCS})
target_link_libraries(LC3plus_shared PRIVATE m)
set_target_properties(LC3plus_shared PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_NAME ${LIB_LC3})

add_custom_target(show_help
    COMMAND ${CMAKE_COMMAND} -E echo "Targets: ${NAME_LC3} (default), ${LIB_LC3}"
    COMMAND ${CMAKE_COMMAND} -E echo "Configure with -D<OPT>=<VAL>:"
    COMMAND ${CMAKE_COMMAND} -E echo "  NO_POST_REL_CHANGES_TEST [ON,OFF]"
    COMMAND ${CMAKE_COMMAND} -E echo "  SHORT_PLC_FADEOUT        [ON,OFF]"
    COMMAND ${CMAKE_COMMAND} -E echo "  GCOV                     [ON,OFF]"
    COMMAND ${CMAKE_COMMAND} -E echo "  WINDOWS                  [ON,OFF]"
    COMMAND ${CMAKE_COMMAND} -E echo "  OPTIM                    [0,1,2,3]"
    COMMAND ${CMAKE_COMMAND} -E echo "  CLANG                    [0,1,2,3]"
    COMMAND ${CMAKE_COMMAND} -E echo "  NAME_LC3                 binary name (default LC3plus)")

