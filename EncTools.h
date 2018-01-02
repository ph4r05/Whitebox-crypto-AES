//
// Created by Dusan Klinec on 30.12.17.
//

#ifndef WHITEBOX_CRYPTO_AES_ENCTOOLS_H
#define WHITEBOX_CRYPTO_AES_ENCTOOLS_H


#include "WBAESGenerator.h"
#include "InputObject.h"

class EncTools {
public:
    /**
     *
     * @param decrypt true/false
     * @param inf [in] input stream to process
     * @param out [in, optional] output stream to return result to
     * @param coding [in, optional] external coding
     * @param cacc [in,out, optional] time_t time accumulator
     * @param pacc [in,out, optional] clock_t time accumulator
     */
    static void processData(bool decrypt, WBAES * wbaes, WBAESGenerator * generator,
                            InputObject<BYTE> * inf, InputObject<BYTE> * out,
                            ExtEncoding * coding = nullptr, bool padding = false, bool cbc = false, BYTE * iv = nullptr,
                            time_t *cacc = nullptr, clock_t * pacc = nullptr);
};


#endif //WHITEBOX_CRYPTO_AES_ENCTOOLS_H
