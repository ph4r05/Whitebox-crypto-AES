//
// Created by Dusan Klinec on 30.12.17.
//

#include "EncTools.h"
#include "GenericAES.h"
#include "WBAESGenerator.h"

void EncTools::processData(bool decrypt, istream &inf, ofstream * out, ExtEncoding * coding,
                           time_t *cacc, clock_t * pacc)
{
    GenericAES defAES;
    defAES.init(0x11B, 0x03);

    WBAESGenerator generator;
    auto * genAES = new WBAES;

    // read the file
    const int buffSize       = 4096;
    const long int iters     = buffSize / N_BYTES;
    unsigned long long blockCount = 0;
    auto * memblock          = new char[buffSize];
    char blockbuff[N_BYTES];

    // time measurement of just the cipher operation
    time_t cstart=0, cend=0;
    if (cacc){
        *cacc = 0;
    }

    clock_t pstart=0, pend=0;
    if (pacc){
        *pacc = 0;
    }

    // measure the time here
    do {
        streamsize bRead;

        // read data from the file to the buffer
        inf.read(memblock, buffSize);
        bRead = inf.gcount();
        if (inf.bad()) {
            break;
        }

        // here we have data in the buffer - lets encrypt them
        W128b state{};
        long int iter2comp = min(iters, (long int) ceil((float)bRead / N_BYTES));

        for(int k = 0; k < iter2comp; k++, blockCount++){
            arr_to_W128b(memblock, k * 16UL, state);

            if (cacc) {
                time(&cstart);
            }
            if(pacc) {
                pstart = clock();
            }

            // encryption
            if (coding) {
                generator.applyExternalEnc(state, coding, true);
            }

            if (decrypt){
                genAES->decrypt(state);
            } else {
                genAES->encrypt(state);
            }

            if (coding) {
                generator.applyExternalEnc(state, coding, false);
            }

            if (pacc) {
                pend = clock();
                *pacc += (pend - pstart);
            }

            if (cacc) {
                time(&cend);
                *cacc += (cend - cstart);
            }

            // if wanted, store to file
            if (out){
                W128b_to_arr(blockbuff, 0, state);
                out->write(blockbuff, N_BYTES);
            }
        }

        if (inf.eof()){
            break;
        }

    } while(true);

    if (out){
        out->flush();
    }

    // free allocated memory
    delete genAES;
    delete[] memblock;
}
