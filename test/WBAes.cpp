#include <gtest/gtest.h>
#include <WBAES.h>
#include <WBAESGenerator.h>

using namespace std;

TEST(WBAES, GenerateKey)
{
    // A given AES key.
    BYTE key[16] = {0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6, 0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c };

    GF2X gf2x = GF2XFromLong(0x11B, 9);
    GF2E::init(gf2x);

    shared_ptr<WBAES> wbaes(new WBAES());
    WBAESGenerator generator;
    ExtEncoding coding;
    generator.generateExtEncoding(&coding, 0);
    generator.generateTables(key, KEY_SIZE_16, wbaes.get(), &coding, true);
    generator.generateTables(key, KEY_SIZE_16, wbaes.get(), &coding, false);

    BYTE plaintext[16] = { 0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34 };
    BYTE ciphertext[16] = { 0x39, 0x25, 0x84, 0x1d, 0x02, 0xdc, 0x09, 0xfb, 0xdc, 0x11, 0x85, 0x97, 0x19, 0x6a, 0x0b, 0x32 };

    W128b plain, cipher, state;
    arr_to_W128b(plaintext, 0, plain);
    arr_to_W128b(plaintext, 0, state);
    arr_to_W128b(ciphertext, 0, cipher);

    generator.applyExternalEnc(state, &coding, true);
    wbaes->encrypt(state);
    generator.applyExternalEnc(state, &coding, false);
    EXPECT_TRUE(compare_W128b(state, cipher));

    generator.applyExternalEnc(state, &coding, true);
    wbaes->decrypt(state);
    generator.applyExternalEnc(state, &coding, false);
    EXPECT_TRUE(compare_W128b(state, plain));

    std::ostringstream out;
    generator.save("X:/key.bin", wbaes.get(), &coding);
}

TEST(WBAES, EncDec)
{
    // test key file can be well-accessed.
    ifstream fl("X:/key.bin");
    fl.seekg(0, ios::end);
    size_t len = fl.tellg();
    EXPECT_TRUE(len == 2663044);
    char *ret = new char[len];
    fl.seekg(0, ios::beg);
    fl.read(ret, len);
    fl.close();

    BYTE plaintext[16] = { 0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d, 0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34 };
    BYTE ciphertext[16] = { 0x39, 0x25, 0x84, 0x1d, 0x02, 0xdc, 0x09, 0xfb, 0xdc, 0x11, 0x85, 0x97, 0x19, 0x6a, 0x0b, 0x32 };

    GF2X gf2x = GF2XFromLong(0x11B, 9);
    GF2E::init(gf2x);

    shared_ptr<WBAES> wbaes(new WBAES());
    WBAESGenerator generator;
    ExtEncoding coding;

    try
    {
        generator.load("X:/key.bin", wbaes.get(), &coding);
    }
    catch (const boost::archive::archive_exception& e)
    {
        printf(e.what());

    }

    W128b plain, cipher, state;
    arr_to_W128b(plaintext, 0, plain);
    arr_to_W128b(plaintext, 0, state);
    arr_to_W128b(ciphertext, 0, cipher);

    generator.applyExternalEnc(state, &coding, true);
    wbaes->encrypt(state);
    generator.applyExternalEnc(state, &coding, false);
    EXPECT_TRUE(compare_W128b(state, cipher));

    generator.applyExternalEnc(state, &coding, true);
    wbaes->decrypt(state);
    generator.applyExternalEnc(state, &coding, false);
    EXPECT_TRUE(compare_W128b(state, plain));
}
