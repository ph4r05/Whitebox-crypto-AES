file(GLOB SRCS *.cpp)

add_executable(testwb ${SRCS})

target_link_libraries(testwb
        Whitebox_crypto_AES
        libgtest
        libgmock
        )

install(TARGETS testfoo DESTINATION bin)
