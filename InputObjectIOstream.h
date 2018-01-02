//
// Created by Dusan Klinec on 02.01.18.
//

#ifndef WHITEBOX_CRYPTO_AES_INPUTOBJECTIOSTREAM_H
#define WHITEBOX_CRYPTO_AES_INPUTOBJECTIOSTREAM_H

#include "InputObject.h"
#include <fstream>
#include <iostream>
#include <stdexcept>

template<typename T>
class InputObjectIOstream : public InputObject<T>{
public:
    explicit InputObjectIOstream(std::iostream *ios) : ios(ios) {}

    ssize_t write(const T *buffer, size_t maxSize) override {
        ios->write((const char*)buffer, maxSize);
        return maxSize;
    }

    ssize_t read(T *buffer, size_t maxSize) override {
        ios->read((char*)buffer, maxSize);
        return ios->gcount();
    }

    bool isGood() override {
        return ios->good();
    }

    bool eof() override {
        return ios->eof();
    }

    bool isFull() override {
        return !ios->good();
    }

protected:
    std::iostream * ios;
};


#endif //WHITEBOX_CRYPTO_AES_INPUTOBJECTIOSTREAM_H
