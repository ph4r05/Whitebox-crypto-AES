//
// Created by Dusan Klinec on 02.01.18.
//

#ifndef WHITEBOX_CRYPTO_AES_INPUTOBJECTBUFFER_H
#define WHITEBOX_CRYPTO_AES_INPUTOBJECTBUFFER_H

#include "InputObject.h"
#include <memory>

template<typename T>
class InputObjectBuffer : public InputObject<T>{
public:
    explicit InputObjectBuffer(size_t size);

    ssize_t write(const T *buffer, size_t maxSize) override;
    ssize_t read(T *buffer, size_t maxSize) override;
    bool eof() override;
    bool isFull() override;
    bool isGood() override;

    T * getBuffer(){
        return this->buffer.get();
    }

    size_t getPos(){
        return this->pos;
    }

    void clear(){
        start = 0;
        pos = 0;
    }

protected:
    size_t size;
    size_t start=0;
    size_t pos=0;
    std::unique_ptr<T[]> buffer = nullptr;
};

template<typename T>
InputObjectBuffer<T>::InputObjectBuffer(size_t size) {
    this->size = size;
    this->pos = 0;
    this->buffer = std::unique_ptr<T[]>(new T[size]);
}

template<typename T>
ssize_t InputObjectBuffer<T>::write(const T *buffer, size_t maxSize) {
    ssize_t toCopy = std::min(maxSize, size-pos);
    if (toCopy > 0) {
        memcpy(this->buffer.get() + pos, buffer, (size_t) toCopy);
        pos += toCopy;
        return toCopy;
    }

    return -1;
}

template<typename T>
ssize_t InputObjectBuffer<T>::read(T *buffer, size_t maxSize) {
    ssize_t toCopy = std::min(maxSize, pos-start);
    if (toCopy > 0) {
        memcpy(buffer, this->buffer.get() + start, (size_t) toCopy);
        start += toCopy;
        return toCopy;
    }

    return -1;
}

template<typename T>
bool InputObjectBuffer<T>::eof() {
    return start >= pos;
}

template<typename T>
bool InputObjectBuffer<T>::isFull() {
    return pos >= size;
}

template<typename T>
bool InputObjectBuffer<T>::isGood() {
    return !eof();
}

#endif //WHITEBOX_CRYPTO_AES_INPUTOBJECTBUFFER_H
