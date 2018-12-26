cmake_minimum_required(VERSION 2.8)
project(Whitebox_crypto_AES)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
option(NTL_NO_DEFAULT_PATH "Using just given paths when searching for NTL" OFF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCE_FILES
    base.h
    base.cpp
    BGEAttack.cpp
    BGEAttack.h
    GenericAES.cpp
    GenericAES.h
    LinearAffineEq.cpp
    LinearAffineEq.h
    LinearAffineEq_test.cpp
    md5.c
    md5.h
    MixingBijections.cpp
    MixingBijections.h
    NTLUtils.cpp
    NTLUtils.h
    WBAES.cpp
    WBAES.h
    WBAESGenerator.cpp
    WBAESGenerator.h
    RingBuffer.h
    RingBuffer.cpp
    EncTools.cpp
    EncTools.h
    InputObject.cpp
    InputObject.h
    InputObjectBuffer.cpp
    InputObjectBuffer.h
    InputObjectIstream.cpp
    InputObjectIstream.h
    InputObjectOstream.cpp
    InputObjectOstream.h
    InputObjectIOstream.cpp
    InputObjectIOstream.h
)

set(TEST_SOURCE_FILES
        tests/WBAes.cpp
        tests/tester.cpp
        tests/Commons.cpp
        tests/Commons.h
        tests/RingBufferTest.cpp
        tests/EncToolsTest.cpp
)

add_executable(main ${SOURCE_FILES} main.cpp)
add_executable(testing ${SOURCE_FILES} testing.cpp testing.h)
add_executable(gtesting ${SOURCE_FILES} ${TEST_SOURCE_FILES})

# NTL
if (NTL_NO_DEFAULT_PATH)
    set(NTL_MODIF NO_DEFAULT_PATH)
    message("Not using default paths for NTL search, ${NTL_MODIF}")
endif (NTL_NO_DEFAULT_PATH)
message("HOME: $ENV{HOME}")
find_path(NTL_INCLUDE_PATH NAMES NTL/mat_GF2.h PATHS $ENV{HOME}/ntl $ENV{HOME}/ntl/include /usr/local /opt/local /usr ${NTL_MODIF})
message("NTL_INCLUDE_PATH: ${NTL_INCLUDE_PATH}")

find_library(NTL_LIB ntl PATHS ${CMAKE_CURRENT_SOURCE_DIR} $ENV{HOME}/ntl $ENV{HOME}/ntl/lib /usr/local /opt/local /usr ${NTL_MODIF})
if(NOT NTL_LIB)
    message(FATAL_ERROR "ntl library not found.  Rerun cmake with -DCMAKE_PREFIX_PATH=\"<path to lib1>;<path to lib2>\"")
endif()
message("NTL_LIB: ${NTL_LIB}")
include_directories(${NTL_INCLUDE_PATH})

# BOOST
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS program_options serialization iostreams random)
include_directories(${Boost_INCLUDE_DIRS})

# We need thread support
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package( Threads )
#find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Download and install GoogleTest
ExternalProject_Add(
        gtest
        URL https://github.com/google/googletest/archive/master.zip
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
        # Disable install step
        INSTALL_COMMAND ""
)

# Get GTest source and binary directories from CMake project
ExternalProject_Get_Property(gtest source_dir binary_dir)

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set libgtest properties
set_target_properties(libgtest PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/lib/libgtest.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
        )

# Create a libgmock target to be used as a dependency by test programs
add_library(libgmock IMPORTED STATIC GLOBAL tests/tester.cpp)
add_dependencies(libgmock gtest)

# Set libgmock properties
set_target_properties(libgmock PROPERTIES
        "IMPORTED_LOCATION" "${binary_dir}/lib/libgmock.a"
        "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
        )

# I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
include_directories(
        "${source_dir}/googletest/include"
        "${source_dir}/googlemock/include")

# Linking
if(NOT MSVC)
    set(PThreadLib -pthread)
endif()

target_link_libraries(main ${NTL_LIB} ${Boost_LIBRARIES} ${PThreadLib})
target_link_libraries(testing ${NTL_LIB} ${Boost_LIBRARIES} ${PThreadLib})
target_link_libraries(gtesting ${NTL_LIB} ${Boost_LIBRARIES} libgtest libgmock ${CMAKE_THREAD_LIBS_INIT} ${PThreadLib})

