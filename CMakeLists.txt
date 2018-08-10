cmake_minimum_required(VERSION 3.0)

project("onf" C)

# The version number.
set(ONF_VERSION "v0.1.0")

# Check and see if OpenMP is available
find_package(OpenMP)
if (OPENMP_FOUND)
    option(OPENMP "OpenMP support" ON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else()
    option(OPENMP "OpenMP support" OFF)
endif()


set(CMAKE_BINARY_DIR ${CMAKE_BINARY_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES ${SOURCE_DIR}/*.c)
file(GLOB VENDOR_INCLUDES ${PROJECT_SOURCE_DIR}/vendor/*)
file(GLOB VENDOR_SOURCES ${PROJECT_SOURCE_DIR}/vendor/*/*.c)

include_directories(${VENDOR_INCLUDES})

add_executable(main ${SOURCES} ${VENDOR_SOURCES})
target_link_libraries(main PRIVATE -lm -lz)

configure_file("${SOURCE_DIR}/main.h.in" "${SOURCE_DIR}/main.h")
