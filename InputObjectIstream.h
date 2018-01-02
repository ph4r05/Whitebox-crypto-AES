//
// Created by Dusan Klinec on 02.01.18.
//

#ifndef WHITEBOX_CRYPTO_AES_INPUTOBJECTIFSTREAM_H
#define WHITEBOX_CRYPTO_AES_INPUTOBJECTIFSTREAM_H

#include "InputObject.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

template<typename T>
class InputObjectIstream : public InputObject<T>{
public:
    explicit InputObjectIstream(std::istream *is) : is(is) {}

    ssize_t write(const T *buffer, size_t maxSize) override {
        throw std::invalid_argument("Istream does not support write");
    }

    ssize_t read(T *buffer, size_t maxSize) override {
        is->read((char*)buffer, maxSize);
        return is->gcount();
    }

    bool eof() override {
        return is->eof();
    }

    bool isFull() override {
        throw std::invalid_argument("Istream does not support write");
    }

    bool isGood() override {
        return is->good();
    }

protected:
    std::istream * is;
};


#endif //WHITEBOX_CRYPTO_AES_INPUTOBJECTIFSTREAM_H
