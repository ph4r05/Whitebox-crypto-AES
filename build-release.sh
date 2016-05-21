#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir -p build-release \
    && cd build-release \
    && cmake -DCMAKE_BUILD_TYPE=Release "${DIR}" \
    && make
