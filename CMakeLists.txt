cmake_minimum_required(VERSION 3.0)

project("onf" C)

# The version number.
set(ONF_VERSION "v0.1.0")

# Modules
include(CheckIncludeFile)
include(GNUInstallDirs)

set(LINKER_FLAGS "-lm -lz")

# Check and see if OpenMP is available
find_package(OpenMP)
if (OPENMP_FOUND)
    option(OPENMP "OpenMP support" ON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else ()
    option(OPENMP "OpenMP support" OFF)
endif ()

set(CMAKE_BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES src/array.c src/file.c src/onf.c src/onf_math.c)
file(GLOB VENDOR_INCLUDES ${PROJECT_SOURCE_DIR}/vendor/*)
file(GLOB VENDOR_SOURCES ${PROJECT_SOURCE_DIR}/vendor/*/*.c)

include_directories(${VENDOR_INCLUDES})

add_executable(build ${SOURCE_DIR}/build.c ${SOURCES} ${VENDOR_SOURCES})
target_link_libraries(build PRIVATE ${LINKER_FLAGS})

add_executable(predict ${SOURCE_DIR}/predict.c ${SOURCES} ${VENDOR_SOURCES})
target_link_libraries(predict PRIVATE ${LINKER_FLAGS})

# Find the rya library
find_library(RYA_FILE_LIB "rya_file")
find_library(RYA_FORMAT_LIB "rya_format")

find_path(RYA_FILE_H_PATH "rya_file.h")
find_path(RYA_FORMAT_H_PATH "rya_format.h")

if (RYA_FILE_LIB AND RYA_FILE_H_PATH AND RYA_FORMAT_LIB AND RYA_FORMAT_H_PATH)
    target_link_libraries(build PRIVATE ${RYA_FILE_LIB})
    target_link_libraries(build PRIVATE ${RYA_FORMAT_LIB})

    target_link_libraries(predict PRIVATE ${RYA_FILE_LIB})
    target_link_libraries(predict PRIVATE ${RYA_FORMAT_LIB})

    include_directories(${RYA_FILE_H_PATH} ${RYA_FORMAT_H_PATH})
else ()
    message(FATAL_ERROR "rya library not found")
endif ()


configure_file("${SOURCE_DIR}/build.h.in" "${SOURCE_DIR}/build.h")
