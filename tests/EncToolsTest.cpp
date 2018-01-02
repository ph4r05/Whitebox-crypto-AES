//
// Created by Dusan Klinec on 02.01.18.
//

#include <gtest/gtest.h>
#include <sstream>
#include <iostream>
#include "../base.h"
#include "../EncTools.h"
#include "../InputObject.h"
#include "../InputObjectBuffer.h"
#include "Commons.h"

using namespace std;

class EncToolsTest : public ::testing::Test {

protected:
    void SetUp() override {
        io1 = std::unique_ptr<InputObjectBuffer<BYTE>>(new InputObjectBuffer<BYTE>(size + 256));
        io2 = std::unique_ptr<InputObjectBuffer<BYTE>>(new InputObjectBuffer<BYTE>(size + 256));
        io3 = std::unique_ptr<InputObjectBuffer<BYTE>>(new InputObjectBuffer<BYTE>(size + 256));
        char buff[16];

        for(int i = 0; i < size; i++){
            buff[0] = (char)('a' + (i % 26));
            io1->write((const BYTE*)buff, 1);
        }

        initDefaultModulus();
        wbaes = unique_ptr<WBAES>(new WBAES());
        generator = unique_ptr<WBAESGenerator>(new WBAESGenerator());
        coding = unique_ptr<ExtEncoding>(new ExtEncoding);
        coding2 = unique_ptr<ExtEncoding>(new ExtEncoding);

        generator->generateExtEncoding(coding.get(), WBAESGEN_EXTGEN_ID);
        generator->generateExtEncoding(coding2.get(), 0);
        generator->generateTables(const_cast<BYTE *>(test_key), KEY_SIZE_16, wbaes.get(), coding.get(), true);
        generator->generateTables(const_cast<BYTE *>(test_key), KEY_SIZE_16, wbaes.get(), coding.get(), false);
    }

    void TearDown() override {

    }

    const size_t size = 8192;
    std::unique_ptr<InputObjectBuffer<BYTE>> io1;
    std::unique_ptr<InputObjectBuffer<BYTE>> io2;
    std::unique_ptr<InputObjectBuffer<BYTE>> io3;

    unique_ptr<WBAES> wbaes;
    unique_ptr<WBAESGenerator> generator;
    unique_ptr<ExtEncoding> coding;
    unique_ptr<ExtEncoding> coding2;
};

TEST_F(EncToolsTest, EncDecBase)
{
    EncTools::processData(false, wbaes.get(), generator.get(), io1.get(), io2.get(), coding.get(), false, false, nullptr, nullptr, nullptr);
    EXPECT_EQ(io2->getPos(), size);

    EncTools::processData(true, wbaes.get(), generator.get(), io2.get(), io3.get(), coding.get(), false, false, nullptr, nullptr, nullptr);
    EXPECT_TRUE(memcmp(io1.get()->getBuffer(), io3.get()->getBuffer(), size) == 0);
}

TEST_F(EncToolsTest, EncDecBaseInvalidSize)
{
    io1->write((BYTE*)"ZZ", 2);
    EXPECT_THROW(
            EncTools::processData(false, wbaes.get(), generator.get(), io1.get(), io2.get(), coding.get(), false, false, nullptr, nullptr, nullptr),
            std::invalid_argument);
}

TEST_F(EncToolsTest, EncDecBasePadding)
{
    EncTools::processData(false, wbaes.get(), generator.get(), io1.get(), io2.get(), coding.get(), true, false, nullptr, nullptr, nullptr);
    EXPECT_EQ(io2->getPos(), size + 16);

    EncTools::processData(true, wbaes.get(), generator.get(), io2.get(), io3.get(), coding.get(), true, false, nullptr, nullptr, nullptr);
    EXPECT_TRUE(memcmp(io1.get()->getBuffer(), io3.get()->getBuffer(), size) == 0);
}

TEST_F(EncToolsTest, EncDecPadding)
{
    io1->write((BYTE*)"XX", 2);
    EncTools::processData(false, wbaes.get(), generator.get(), io1.get(), io2.get(), coding.get(), true, false, nullptr, nullptr, nullptr);
    EXPECT_EQ(io2->getPos(), size + 16);

    EncTools::processData(true, wbaes.get(), generator.get(), io2.get(), io3.get(), coding.get(), true, false, nullptr, nullptr, nullptr);
    EXPECT_TRUE(memcmp(io1.get()->getBuffer(), io3.get()->getBuffer(), size + 2) == 0);
}

TEST_F(EncToolsTest, EncDecPaddingCbcZeroIv)
{
    io1->write((BYTE*)"XX", 2);
    EncTools::processData(false, wbaes.get(), generator.get(), io1.get(), io2.get(), coding.get(), true, true, nullptr, nullptr, nullptr);
    EXPECT_EQ(io2->getPos(), size + 16);

    EncTools::processData(true, wbaes.get(), generator.get(), io2.get(), io3.get(), coding.get(), true, true, nullptr, nullptr, nullptr);
    EXPECT_TRUE(memcmp(io1.get()->getBuffer(), io3.get()->getBuffer(), size + 2) == 0);
}

TEST_F(EncToolsTest, EncDecPaddingCbcIv)
{
    BYTE iv[N_BYTES];
    for(int i = 0; i < N_BYTES; i++){
        iv[i] = (BYTE) i;
    }

    io1->write((BYTE*)"XX", 2);
    EncTools::processData(false, wbaes.get(), generator.get(), io1.get(), io2.get(), coding.get(), true, true, iv, nullptr, nullptr);
    EXPECT_EQ(io2->getPos(), size + 16);

    EncTools::processData(true, wbaes.get(), generator.get(), io2.get(), io3.get(), coding.get(), true, true, iv, nullptr, nullptr);
    EXPECT_TRUE(memcmp(io1.get()->getBuffer(), io3.get()->getBuffer(), size + 2) == 0);
}

