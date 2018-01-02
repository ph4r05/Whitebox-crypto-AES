//
// Created by Dusan Klinec on 02.01.18.
//

#ifndef WHITEBOX_CRYPTO_AES_INPUTOBJECT_H
#define WHITEBOX_CRYPTO_AES_INPUTOBJECT_H

#include <ios>

/**
 * Input output object wrapper
 * @tparam T
 */
template<typename T>
class InputObject {
public:
    virtual ssize_t write(const T * buffer, size_t maxSize)=0;
    virtual ssize_t read(T * buffer, size_t maxSize)=0;
    virtual bool isGood()=0;
    virtual bool eof()=0;
    virtual bool isFull()=0;
};


#endif //WHITEBOX_CRYPTO_AES_INPUTOBJECT_H
