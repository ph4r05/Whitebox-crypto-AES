//
// Created by Dusan Klinec on 02.01.18.
//

#ifndef WHITEBOX_CRYPTO_AES_INPUTOBJECTOFSTREAM_H
#define WHITEBOX_CRYPTO_AES_INPUTOBJECTOFSTREAM_H

#include "InputObject.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

template<typename T>
class InputObjectOstream  : public InputObject<T>{
public:
    explicit InputObjectOstream(std::ostream *os) : os(os) {}

    ssize_t write(const T *buffer, size_t maxSize) override {
        os->write((const char*)buffer, maxSize);
        return maxSize;
    }

    ssize_t read(T *buffer, size_t maxSize) override {
        throw std::invalid_argument("Ostream does not support read");
    }

    bool isGood() override {
        return os->good();
    }

    bool eof() override {
        return os->eof();
    }

    bool isFull() override {
        return !os->good();
    }

protected:
    std::ostream * os;
};


#endif //WHITEBOX_CRYPTO_AES_INPUTOBJECTOFSTREAM_H
