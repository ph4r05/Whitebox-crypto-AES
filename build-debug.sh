#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p build-debug \
    && cd build-debug \
    && cmake -DCMAKE_BUILD_TYPE=Debug "${DIR}" \
    && make
