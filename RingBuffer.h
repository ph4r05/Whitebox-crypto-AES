//
// Created by Dusan Klinec on 31.12.17.
//

#ifndef WHITEBOX_CRYPTO_AES_RINGBUFFER_H
#define WHITEBOX_CRYPTO_AES_RINGBUFFER_H

#include <iostream>
#include <fstream>
#include <iosfwd>
#include <ios>
#include <memory>
#include "InputObject.h"
#include <boost/circular_buffer.hpp>

template<typename T>
class RingBuffer {
public:
    typedef RingBuffer<T> this_type;
    typedef std::pair<T*, size_t> array_range;
    typedef std::pair<array_range, array_range> array_ranges;

    explicit RingBuffer(size_t size);

    /**
     *
     * @return true if empty
     */
    bool isEmpty();

    /**
     *
     * @return true if full
     */
    bool isFull();

    /**
     *
     * @return Number of bytes stored in the ring buffer
     */
    size_t getBytesAvailable();

    /**
     *
     * @return Number of the free space bytes in the ring buffer
     */
    size_t getSpaceAvailable();

    /**
     * Read bytes from the ring buffer to the provided buffer.
     *
     * @param buffer buffer to write ring buffer data to
     * @param maxLength the maximal length of the destination buffer
     * @return number of bytes written to the buffer
     */
    ssize_t read(T * buffer, size_t maxLength);
    ssize_t read(InputObject<T> * buffer, ssize_t maxLength);

    /**
     * Write data to the circular buffer from the provided buffer
     *
     * @param buffer buffer to read data from
     * @param maxLength maximal number of available bytes in the read buffer
     * @return number of bytes read from the buffer
     */
    ssize_t write(T const *buffer, size_t maxLength);
    ssize_t write(InputObject<T> *buffer, size_t maxLength);

    /**
    * If ring buffer is empty, resets writing target to the beginning of the buffer
    * and returns internal buffer array.
    *
    * If buffer is not empty, NULL is returned.
    */
    T * resetBufferIfEmpty();

    /**
    * Return current native reading buffer part, with number of bytes available in this reading section.
    * This is reading ring buffer, thus for complete read there may be 2 calls needed.
    */
    T * getReadingBuffer(size_t * bytesAvailable);

    /**
    * Return length of the contiguous buffer for writing.
    * Since we are using ring buffer it may wrap around.
    */
    ssize_t getContiguousWriteBufferLen();

    /**
    * Return length of the contiguous buffer for reading.
    * Since we are using ring buffer it may wrap around.
    */
    ssize_t getContiguousReadBufferLen();

    /**
    * Returns contiguous buffer for writing even if ring buffer is not empty.
    */
    T * getContiguousWriteBuffer();

    /**
    * Returns contiguous buffer for reading even if ring buffer is not empty.
    */
    T * getContiguousReadBuffer();

    /**
    * Designed to be used in conjunction with resetBufferIfEmpty.
    * This should be called when some data is written natively to the buffer using resetBufferIfEmpty
    * and direct pointer to the memory.
    */
    bool setBytesWritten(ssize_t bytesWritten);

    /**
    * Designed to be used in conjunction with getReadingBuffer.
    * This should be called whe some data is read from native pointer using getReadingBuffer
    * to advance internal pointers.
    */
    bool setBytesRead(ssize_t bytesRead);

    /**
     * Returns 2 consecutive arrays for writing to the ring buffer
     * @return
     */
    array_ranges array_free();

    /**
     * Returns 2 consecutive arrays for reading from the ring buffer
     * @return
     */
    array_ranges array_read();

private:
    size_t _buffSize = 0;
    std::unique_ptr<T[]> _buffer = nullptr;
    size_t _fill  = 0;
    size_t _use   = 0;
    size_t _count = 0;
};


template<typename T>
RingBuffer<T>::RingBuffer(size_t size) {
    this->_buffer = std::unique_ptr<T[]>(new T[size]);
    this->_buffSize = size;
    this->_fill = 0;
    this->_use = 0;
    this->_count = 0;
}

template<typename T>
bool RingBuffer<T>::isEmpty() {
    return _count == 0;
}

template<typename T>
bool RingBuffer<T>::isFull() {
    return _count == _buffSize;
}

template<typename T>
size_t RingBuffer<T>::getBytesAvailable() {
    return _count;
}

template<typename T>
size_t RingBuffer<T>::getSpaceAvailable() {
    return _buffSize - _count;
}

template<typename T>
ssize_t RingBuffer<T>::read(T *buffer, size_t maxLength) {
    if (_count == 0){
        return 0;
    }

    ssize_t read = 0;
    array_ranges ars = this->array_read();

    // Number of bytes that can be read in this call.
    ssize_t bytesToRead = maxLength < 0 ? _count : std::min(_count, static_cast<size_t>(maxLength));

    if (ars.first.second > 0 && bytesToRead > 0){
        ssize_t toRead = std::min(bytesToRead, (ssize_t)ars.first.second);
        memcpy(buffer, ars.first.first, (size_t)toRead);
        read += toRead;
        bytesToRead -= toRead;
    }

    if (ars.second.second > 0 && bytesToRead > 0){
        ssize_t toRead = std::min(bytesToRead, (ssize_t)ars.first.second);
        memcpy(buffer+read, ars.second.first, (size_t)toRead);
        read += toRead;
    }
    setBytesRead(read);

    return read;
}

template<typename T>
ssize_t RingBuffer<T>::read(InputObject<T> * buffer, ssize_t maxLength) {
    if (_count == 0){
        return 0;
    }

    ssize_t read = 0;
    array_ranges ars = this->array_read();

    // Number of bytes that can be read in this call.
    ssize_t bytesToRead = maxLength < 0 ? _count : std::min(_count, static_cast<size_t>(maxLength));

    if (ars.first.second > 0 && bytesToRead > 0 && buffer->isGood()){
        ssize_t toRead = std::min(bytesToRead, (ssize_t)ars.first.second);
        buffer->write((const T*)ars.first.first, (size_t)toRead);
        read += toRead;
        bytesToRead -= toRead;
    }

    if (ars.second.second > 0 && bytesToRead > 0 && buffer->isGood()){
        ssize_t toRead = std::min(bytesToRead, (ssize_t)ars.first.second);
        buffer->write((const T*)ars.second.first, (size_t)toRead);
        read += toRead;
    }
    setBytesRead(read);

    return read;
}

template<typename T>
ssize_t RingBuffer<T>::write(T const *buffer, size_t maxLength) {
    if (isFull()){
        return 0;
    }

    ssize_t written = 0;
    array_ranges ars = this->array_free();

    // Number of bytes that will be written in this call.
    ssize_t bytesToWrite = maxLength < 0 ? getSpaceAvailable() : std::min(getSpaceAvailable(), maxLength);

    if (ars.first.second > 0 && bytesToWrite > 0){
        ssize_t toWrite = std::min(bytesToWrite, (ssize_t)ars.first.second);
        memcpy(ars.first.first, buffer + written, (size_t)toWrite);
        written += toWrite;
        bytesToWrite -= toWrite;
    }

    if (ars.second.second > 0 && bytesToWrite > 0){
        ssize_t toWrite = std::min(bytesToWrite, (ssize_t)ars.second.second);
        memcpy(ars.second.first, buffer + written, (size_t)bytesToWrite);
        written += toWrite;
        bytesToWrite -= toWrite;
    }

    setBytesWritten(written);
    return written;
}

template<typename T>
ssize_t RingBuffer<T>::write(InputObject<T> *buffer, size_t maxLength) {
    if (isFull()){
        return 0;
    }

    ssize_t written = 0;
    array_ranges ars = this->array_free();

    // Number of bytes that will be written in this call.
    ssize_t bytesToWrite = maxLength < 0 ? getSpaceAvailable() : std::min(getSpaceAvailable(), maxLength);

    if (ars.first.second > 0 && bytesToWrite > 0 && buffer->isGood()){
        ssize_t toWrite = std::min(bytesToWrite, (ssize_t)ars.first.second);
        ssize_t actuallyWritten = buffer->read((T*)ars.first.first, (size_t)toWrite);
        written += actuallyWritten;
        bytesToWrite -= actuallyWritten;
    }

    if (ars.second.second > 0 && bytesToWrite > 0 && buffer->isGood()){
        ssize_t toWrite = std::min(bytesToWrite, (ssize_t)ars.second.second);
        ssize_t actuallyWritten = buffer->read((T*)ars.second.first, (size_t)toWrite);
        written += actuallyWritten;
        bytesToWrite -= actuallyWritten;
    }

    setBytesWritten(written);
    return written;
}

template<typename T>
T *RingBuffer<T>::resetBufferIfEmpty() {
    if (isEmpty()){
        return NULL;
    }

    _fill = 0;
    _use  = 0;
    return _buffer.get();
}

template<typename T>
T *RingBuffer<T>::getReadingBuffer(size_t *bytesAvailable) {
    if (isEmpty()){
        return nullptr;
    }

    if (bytesAvailable != nullptr){
        *bytesAvailable = (size_t) getContiguousReadBufferLen();
    }

    return getContiguousReadBuffer();
}

template<typename T>
ssize_t RingBuffer<T>::getContiguousWriteBufferLen() {
    // Minimum from (space available in buffer, space to the right boundary).
    // First parameter for cases like   |    FeeeU    |, e=empty to use.
    // Second parameter for cases like  |    U   Feeee|, e=empty to use.
    return std::min(_buffSize - _count, _buffSize - _fill);
}

template<typename T>
ssize_t RingBuffer<T>::getContiguousReadBufferLen() {
    // Minimum from (bytes available in buffer, bytes to the right boundary).
    // First parameter for cases like   |    UxxxF    |, e=bytes to use.
    // Second parameter for cases like  |    F   Uxxxx|, e=bytes to use.
    return std::min(_count, _buffSize - _use);
}

template<typename T>
T *RingBuffer<T>::getContiguousWriteBuffer() {
    return _buffer.get() + _fill;
}

template<typename T>
T *RingBuffer<T>::getContiguousReadBuffer() {
    return _buffer.get() + _use;
}

template<typename T>
bool RingBuffer<T>::setBytesWritten(ssize_t bytesWritten) {
    if (_count + bytesWritten > _buffSize){
        return false;
    }

    _fill   = (_fill + bytesWritten) % _buffSize; // Should never need to modulo.
    _count += bytesWritten;
    assert((_count != _buffSize || _use == _fill) && "Cyclic buffer invariant failed");
    return true;
}

template<typename T>
bool RingBuffer<T>::setBytesRead(ssize_t bytesRead) {
    if (_count < bytesRead){
        return false;
    }

    _use    = (_use + bytesRead) % _buffSize; // Should never need to modulo.
    _count -= bytesRead;
    assert((_count != 0 || _use == _fill) && "Cyclic buffer invariant failed");
    return true;
}

template<typename T>
typename RingBuffer<T>::array_ranges RingBuffer<T>::array_read() {
    RingBuffer::array_range r1 = RingBuffer::array_range(nullptr, 0);
    RingBuffer::array_range r2 = RingBuffer::array_range(nullptr, 0);

    if (_count == 0){
        return RingBuffer::array_ranges(r1, r2);
    }

    T * bytes = _buffer.get();

    // Number of bytes that can be read in this call.
    ssize_t bytesToRead = _count;

    // Number of bytes that can be copied right after _fill position without need to modulo.
    const ssize_t bytesToReadRight = std::min(bytesToRead, static_cast<ssize_t>(_buffSize - _use));

    // Phase 1 - memcpy up to the right boundary.
    size_t nuse = _use;
    if (bytesToReadRight > 0){
        r1 = RingBuffer::array_range(bytes + _use, bytesToReadRight);

        nuse         = (_use + bytesToReadRight) % _buffSize;  // Should be 0.
        bytesToRead -= bytesToReadRight;
    }

    // Phase 2 - write rest of bytes after modular rotation.
    if (bytesToRead > 0){
        r2 = RingBuffer::array_range(bytes + nuse, bytesToRead);
    }

    return RingBuffer::array_ranges(r1, r2);
}

template<typename T>
typename RingBuffer<T>::array_ranges RingBuffer<T>::array_free() {
    RingBuffer::array_range r1 = RingBuffer::array_range(nullptr, 0);
    RingBuffer::array_range r2 = RingBuffer::array_range(nullptr, 0);

    if (isFull()){
        return RingBuffer::array_ranges(r1, r2);
    }

    T * bytes = _buffer.get();

    // Free space in ring buffer.
    const ssize_t freeSpace = _buffSize - _count;

    // Number of bytes that will be written in this call.
    ssize_t bytesToWrite = freeSpace;

    // Number of bytes that can be copied right after _fill position without need to modulo.
    const ssize_t bytesToWriteRight = std::min(bytesToWrite, static_cast<ssize_t>(_buffSize - _fill));

    // Phase 1 - memcpy to the right boundary.
    size_t nfill = _fill;
    if (bytesToWriteRight > 0) {
        r1 =RingBuffer::array_range(bytes + _fill, bytesToWriteRight);

        bytesToWrite -= bytesToWriteRight;
        nfill = (_fill + bytesToWriteRight) % _buffSize;  // Should be 0, if we reached the right boundary.
    }

    // Phase 2 - write rest of bytes after modular rotation.
    if (bytesToWrite > 0) {
        r2 = RingBuffer::array_range(bytes + nfill, bytesToWrite);
    }

    return RingBuffer::array_ranges(r1, r2);
}

//template<typename T>
//class RingBuffer : public boost::circular_buffer<T> {
//public:
//    using base_t = boost::circular_buffer<T>; // The original parent type
//    typename base_t::array_range array_free_one();
//    typename base_t::array_range array_free_two();
//};
//
//template<typename T>
//typename boost::circular_buffer<T>::array_range RingBuffer<T>::array_free_one() {
//    typename base_t::array_range ar1 = this->array_one();
//    typename base_t::array_range ar2 = this->array_two();
//    typename base_t::iterator it = this->begin();
//    return boost::circular_buffer<T, std::__1::allocator<T>>::array_range();
//}
//
//template<typename T>
//typename boost::circular_buffer<T>::array_range RingBuffer<T>::array_free_two() {
//    return boost::circular_buffer<T, std::__1::allocator<T>>::array_range();
//}
#endif //WHITEBOX_CRYPTO_AES_RINGBUFFER_H
