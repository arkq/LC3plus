#/******************************************************************************
#*                        ETSI TS 103 634 V1.3.1                               *
#*              Low Complexity Communication Codec Plus (LC3plus)              *
#*                                                                             *
#* Copyright licence is solely granted through ETSI Intellectual Property      *
#* Rights Policy, 3rd April 2019. No patent licence is granted by implication, *
#* estoppel or otherwise.                                                      *
#******************************************************************************/

#########################################################################################################
# How to clean all cmake cache files? Use:
# find . -iwholename '*cmake*' -not -name CMakeLists.txt -delete
#########################################################################################################

cmake_minimum_required(VERSION 3.13)
project(LC3plus_fec_control_unit_lib VERSION 1.0.1 DESCRIPTION "Forward error correction control unit library for LC3plus codec")

set(OPTION_SETUP "1" CACHE STRING "Byte slots setup")

if ((OPTION_SETUP GREATER 3) OR (OPTION_SETUP LESS 0))
    message( FATAL_ERROR "Setup variable must be in range [0:3]")
endif()

add_library(LC3plus_fec_control_unit STATIC
    fec_control_unit_config.c
    fec_control_unit_ferEstim.c
    fec_control_unit_ferReductionEstim.c
    fec_control_unit_meanEstim.c
    fec_control_unit.c
)

set_target_properties(LC3plus_fec_control_unit PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(LC3plus_fec_control_unit PROPERTIES PUBLIC_HEADER fec_control_unit.h)

target_compile_definitions(LC3plus_fec_control_unit PRIVATE SETUP=${OPTION_SETUP})

if(MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 /WX")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-long-long -pedantic -Werror -g -O3")
endif()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSETUP=${OPTION_SETUP} -I ../../src/fixed_point/")

MESSAGE("")
MESSAGE(STATUS "--------------------------------------------------------------")
MESSAGE(STATUS "Compiler Options:")
MESSAGE(STATUS "-----------------")
MESSAGE(STATUS "C compiler: " ${CMAKE_C_COMPILER})
MESSAGE(STATUS "C flags   : " ${CMAKE_C_FLAGS})
MESSAGE(STATUS "--------------------------------------------------------------")
MESSAGE(STATUS "Build Options:")
MESSAGE(STATUS "Byte slots setup: " ${OPTION_SETUP})
MESSAGE(STATUS "Control with -DOPTION_SETUP=VAR")
MESSAGE(STATUS "The following setups are currently supported:")
MESSAGE(STATUS "   0 : NB   with 40 byte slots")
MESSAGE(STATUS "   1 : WB   with 40 byte slots (default)")
MESSAGE(STATUS "   2 : SSWB with 40 byte slots")
MESSAGE(STATUS "   3 : SWB  with 80 byte slots")
MESSAGE(STATUS "--------------------------------------------------------------")
MESSAGE("")
