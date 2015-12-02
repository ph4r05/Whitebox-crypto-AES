cmake_minimum_required(VERSION 2.8)
project(Whitebox_crypto_AES)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#find_package(Ntl)
#include_directories(${Ntl_INCLUDE_DIRS})

set(SOURCE_FILES
    base.h
    BGEAttack.cpp
    BGEAttack.h
    BGEAttack_test.cpp
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
    WBAESGenerator.h)

add_executable(main ${SOURCE_FILES} main.cpp)
add_executable(testing ${SOURCE_FILES} testing.cpp testing.h)

# NTL
find_library(NTL_LIB ntl)
if(NOT NTL_LIB)
    message(FATAL_ERROR "ntl library not found.  Rerun cmake with -DCMAKE_PREFIX_PATH=\"<path to lib1>;<path to lib2>\"")
endif()
message(NTL_LIB)

link_directories("/usr/local/lib")
include_directories("/usr/local/include")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS program_options serialization iostreams)
include_directories(${Boost_INCLUDE_DIRS})

target_link_libraries(main ${NTL_LIB} ${Boost_LIBRARIES})
target_link_libraries(testing ${NTL_LIB} ${Boost_LIBRARIES})
