//
// Created by Dusan Klinec on 02.01.18.
//

#include <gtest/gtest.h>
#include <sstream>
#include "../RingBuffer.h"
#include "../InputObject.h"
#include "../InputObjectBuffer.h"
#include "Commons.h"

using namespace std;

class RingBufferTest : public ::testing::Test {
protected:
    void SetUp() override {
        size = 16;
        buffer = unique_ptr<RingBuffer<BYTE>>(new RingBuffer<BYTE>(size));
        byteBuffer1 = unique_ptr<BYTE[]>(new BYTE[size]);
        byteBuffer2 = unique_ptr<BYTE[]>(new BYTE[size]);
        bbuf1 = byteBuffer1.get();
        bbuf2 = byteBuffer2.get();

        memset((void *) bbuf1, 0xff, size);
        memset((void *) bbuf2, 0xaa, size);
    }

    void TearDown() override {

    }

    size_t size;
    unique_ptr<RingBuffer<BYTE>> buffer;
    unique_ptr<BYTE[]> byteBuffer1;
    unique_ptr<BYTE[]> byteBuffer2;
    BYTE * bbuf1;
    BYTE * bbuf2;
};


TEST_F(RingBufferTest, BaseAssertions)
{
    EXPECT_EQ(buffer->isEmpty(), true);
    EXPECT_EQ(buffer->isFull(), false);
    EXPECT_EQ(buffer->getBytesAvailable(), 0);
    EXPECT_EQ(buffer->getSpaceAvailable(), size);
}

TEST_F(RingBufferTest, BasicReadWrite)
{
    EXPECT_EQ(buffer->read(bbuf1, size), 0);
    EXPECT_EQ(buffer->write(bbuf1, 8), 8);
    EXPECT_EQ(buffer->getBytesAvailable(), 8);
    EXPECT_EQ(buffer->getContiguousReadBufferLen(), 8);
    EXPECT_TRUE(memcmp(bbuf1, buffer->getContiguousReadBuffer(), 8) == 0);
    EXPECT_EQ(buffer->write(bbuf2, 16), 8);
    EXPECT_TRUE(buffer->isFull());
    EXPECT_FALSE(buffer->isEmpty());

    EXPECT_EQ(buffer->read(bbuf2, 3), 3);
    EXPECT_EQ(buffer->read(bbuf2, 3), 3);
    EXPECT_EQ(buffer->getBytesAvailable(), size - 6);
    EXPECT_EQ(buffer->write(bbuf1, 16), 6);
    EXPECT_TRUE(memcmp(bbuf1, buffer->getContiguousReadBuffer(), 2) == 0);
    EXPECT_TRUE(*(buffer->getContiguousReadBuffer()+1) == 0xff);
    EXPECT_TRUE(*(buffer->getContiguousReadBuffer()+2) == 0xaa);
    EXPECT_EQ(buffer->read(bbuf2, size), size);
    EXPECT_FALSE(buffer->isFull());
    EXPECT_TRUE(buffer->isEmpty());
}

TEST_F(RingBufferTest, BasicStreamReadWrite)
{
    InputObjectBuffer<BYTE> iob(1024);
    iob.write(reinterpret_cast<const BYTE *>(bbuf1), 16);
    iob.write(reinterpret_cast<const BYTE *>(bbuf1), 16);

    EXPECT_EQ(buffer->write(&iob, 10), 10);
    EXPECT_EQ(buffer->write(&iob, 10), size - 10);
    EXPECT_TRUE(buffer->isFull());
    EXPECT_FALSE(buffer->isEmpty());

    EXPECT_EQ(buffer->read(&iob, 1), 1);
    EXPECT_EQ(buffer->read(&iob, 99), size - 1);
}
