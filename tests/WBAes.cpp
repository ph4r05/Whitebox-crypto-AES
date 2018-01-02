//
// Created by Dusan Klinec on 28.12.17.
//

#include <gtest/gtest.h>
#include "../WBAES.h"
#include "../WBAESGenerator.h"
#include "../EncTools.h"
#include "Commons.h"

using namespace std;

TEST(WBAES, GenerateKey)
{
    initDefaultModulus();

    shared_ptr<WBAES> wbaes(new WBAES());
    WBAESGenerator generator;
    ExtEncoding coding;
    generator.generateExtEncoding(&coding, 0);
    generator.generateTables(const_cast<BYTE *>(test_key), KEY_SIZE_16, wbaes.get(), &coding, true);
    generator.generateTables(const_cast<BYTE *>(test_key), KEY_SIZE_16, wbaes.get(), &coding, false);

    W128b plain{}, cipher{}, state{};
    arr_to_W128b(const_cast<BYTE *>(test_plaintext), 0, plain);
    arr_to_W128b(const_cast<BYTE *>(test_plaintext), 0, state);
    arr_to_W128b(const_cast<BYTE *>(test_ciphertext), 0, cipher);

    generator.applyExternalEnc(state, &coding, true);
    wbaes->encrypt(state);
    generator.applyExternalEnc(state, &coding, false);
    EXPECT_TRUE(compare_W128b(state, cipher));

    generator.applyExternalEnc(state, &coding, true);
    wbaes->decrypt(state);
    generator.applyExternalEnc(state, &coding, false);
    EXPECT_TRUE(compare_W128b(state, plain));

    std::stringstream inout;
    generator.save(inout, wbaes.get(), &coding);
}

TEST(WBAES, EncDec)
{
    initDefaultModulus();
    std::stringstream inout;

    // Generating part
    // Different scope to prevent unwanted access after generation.
    {
        shared_ptr<WBAES> wbaes(new WBAES());
        WBAESGenerator generator;
        ExtEncoding coding;

        generator.generateExtEncoding(&coding, 0);
        generator.generateTables(const_cast<BYTE *>(test_key), KEY_SIZE_16, wbaes.get(), &coding, true);
        generator.generateTables(const_cast<BYTE *>(test_key), KEY_SIZE_16, wbaes.get(), &coding, false);
        generator.save(inout, wbaes.get(), &coding);
    }

    // Load phase
    shared_ptr<WBAES> wbaes2(new WBAES());
    WBAESGenerator generator2;
    ExtEncoding coding2;

    generator2.load(inout, wbaes2.get(), &coding2);

    W128b plain{}, cipher{}, state{};
    arr_to_W128b(const_cast<BYTE *>(test_plaintext), 0, plain);
    arr_to_W128b(const_cast<BYTE *>(test_plaintext), 0, state);
    arr_to_W128b(const_cast<BYTE *>(test_ciphertext), 0, cipher);

    generator2.applyExternalEnc(state, &coding2, true);
    wbaes2->encrypt(state);
    generator2.applyExternalEnc(state, &coding2, false);
    EXPECT_TRUE(compare_W128b(state, cipher));

    generator2.applyExternalEnc(state, &coding2, true);
    wbaes2->decrypt(state);
    generator2.applyExternalEnc(state, &coding2, false);
    EXPECT_TRUE(compare_W128b(state, plain));
}

TEST(WBAES, TestVectors)
{
    initDefaultModulus();

    WBAESGenerator generator;
    shared_ptr<WBAES> wbaes(new WBAES());

    // Test WB AES with test vectors.
    // This test also demonstrates usage of external encodings by wrapping AES
    int errors = generator.testWithVectors(false, wbaes.get());
    EXPECT_EQ(errors, 0);

    // Test WB AES with test vectors - no external encoding.
    errors = generator.testWithVectors(false, wbaes.get(), WBAESGEN_EXTGEN_ID);
    EXPECT_EQ(errors, 0);
}

