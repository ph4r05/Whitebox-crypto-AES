//
// Created by Dusan Klinec on 30.12.17.
//

#include <stdexcept>
#include "EncTools.h"
#include "GenericAES.h"
#include "WBAESGenerator.h"
#include "RingBuffer.h"

void EncTools::processData(bool decrypt, WBAES * wbaes, WBAESGenerator * generator,
                           istream * inf, ostream * out, ExtEncoding * coding, bool padding,
                           time_t *cacc, clock_t * pacc)
{
    // read the file
    const int buffSize       = 4096;
    const long int iters     = buffSize / N_BYTES;
    const unique_ptr<RingBuffer<BYTE>> buffer = unique_ptr<RingBuffer<BYTE>>(new RingBuffer<BYTE>(buffSize + N_BYTES));

    unsigned long long blockCount = 0;
    char blockbuff[N_BYTES];
    bool paddingOk = true;

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
    bool eof = false;
    do {
        // read data from the file to the buffer
        streamsize bRead = buffer->write(inf, buffSize);
        streamsize bufferSize = buffer->getBytesAvailable();
        eof = inf->eof();
        if (inf->bad()) {
            break;
        }

        // here we have data in the buffer - lets encrypt them
        W128b state{};
        auto blocks = (float)bufferSize / N_BYTES;
        auto blocks_rounded = (long int) (eof ? ceil(blocks) : floor(blocks));
        auto iter2comp = min(iters, blocks_rounded);

        // Padding not enabled but input data is not block aligned - exception.
        if (eof && (!padding || decrypt) && (long int)blocks != blocks_rounded){
            throw std::invalid_argument("Error: Padding not enabled, input data not block aligned");
        }

        // Add PKCS5 padding bytes to the buffer so we are block aligned
        if (eof && padding && !decrypt){
            auto missingBytes = blocks_rounded * N_BYTES - bufferSize;
            assert(missingBytes > 0 && missingBytes <= 16 && "Padding size is invalid");
            memset(blockbuff, (char)missingBytes, (size_t)missingBytes);
            buffer->read((BYTE*)blockbuff, (size_t)missingBytes);
        }

        for(int k = 0; k < iter2comp; k++, blockCount++){
            buffer->read((BYTE*)blockbuff, 16);
            arr_to_W128b(blockbuff, k * 16UL, state);

            if (cacc) {
                time(&cstart);
            }

            if(pacc) {
                pstart = clock();
            }

            // encryption
            if (coding && generator) {
                generator->applyExternalEnc(state, coding, true);
            }

            if (decrypt){
                wbaes->decrypt(state);
            } else {
                wbaes->encrypt(state);
            }

            if (coding && generator) {
                generator->applyExternalEnc(state, coding, false);
            }

            if (pacc) {
                pend = clock();
                *pacc += (pend - pstart);
            }

            if (cacc) {
                time(&cend);
                *cacc += (cend - cstart);
            }

            // result of the cipher operation
            if (out){
                W128b_to_arr(blockbuff, 0, state);
                ssize_t writeBytes = N_BYTES;

                // If is the last decryption block with the padding - remove the padding.
                if (decrypt && padding && eof && k + 1 < iter2comp){
                    char paddingVal = blockbuff[N_BYTES - 1];
                    if (paddingVal <= 0 || paddingVal > N_BYTES){
                        paddingOk = false;
                    }

                    for(int px = N_BYTES - paddingVal; paddingOk && px < N_BYTES; ++px){
                        paddingOk &= blockbuff[px] == paddingVal;
                    }

                    if (paddingOk){
                        writeBytes -= paddingVal;
                    }
                }

                if (writeBytes > 0) {
                    out->write(blockbuff, (size_t)writeBytes);
                }
            }
        }

    } while(!eof);

    if (out){
        out->flush();
    }

    if (!paddingOk){
        throw std::invalid_argument("Padding is not OK");
    }
}
