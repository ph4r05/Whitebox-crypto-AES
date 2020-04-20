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

# CMakeLists-gtest.txt.in
# Download and unpack googletest at configure time
configure_file(CMakeLists-gtest.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
)
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download"
)

# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds the following targets:
# gtest, gtest_main, gmock and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build"
)

if(CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include"
                        "${gmock_SOURCE_DIR}/include"
    )
endif()

# Linking
if(NOT MSVC)
    set(PThreadLib -pthread)
endif()

target_link_libraries(main ${NTL_LIB} ${Boost_LIBRARIES} ${PThreadLib})
target_link_libraries(testing ${NTL_LIB} ${Boost_LIBRARIES} ${PThreadLib})
target_link_libraries(gtesting ${NTL_LIB} ${Boost_LIBRARIES} gtest_main ${CMAKE_THREAD_LIBS_INIT} ${PThreadLib})
add_test(NAME WBTests COMMAND gtesting)

