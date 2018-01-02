//
// Created by Dusan Klinec on 30.12.17.
//

#include <stdexcept>
#include "EncTools.h"
#include "GenericAES.h"
#include "WBAESGenerator.h"
#include "RingBuffer.h"

void EncTools::processData(bool decrypt, WBAES * wbaes, WBAESGenerator * generator,
                           InputObject<BYTE> * inf, InputObject<BYTE> * out,
                           ExtEncoding * coding, bool padding, bool cbc, BYTE * iv,
                           time_t *cacc, clock_t * pacc)
{
    // read the file
    const int buffSize       = 4096;
    const long int iters     = buffSize / N_BYTES;
    const unique_ptr<RingBuffer<BYTE>> buffer = unique_ptr<RingBuffer<BYTE>>(new RingBuffer<BYTE>(buffSize + N_BYTES));

    unsigned long long blockCount = 0;
    char blockbuff[N_BYTES];
    char blockBuffOut[N_BYTES];
    char prevBlock[N_BYTES] = {0};
    bool paddingOk = true;

    // IV cbc decryption init
    if (cbc){
        iv ? memcpy(prevBlock, iv, N_BYTES) : memset(prevBlock, 0, N_BYTES);
    }

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
        if (!eof && !inf->isGood()) {
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
            auto missingBytes = (blocks_rounded * N_BYTES - bufferSize);
            auto paddingByte = missingBytes == 0 ? 16 : missingBytes;
            assert(paddingByte > 0 && paddingByte <= 16 && "Padding size is invalid");
            memset(blockbuff, (char)paddingByte, (size_t)paddingByte);
            buffer->write((BYTE*)blockbuff, (size_t)paddingByte);
            iter2comp += missingBytes == 0 ? 1 : 0;
        }

        for(int k = 0; k < iter2comp; k++, blockCount++){
            buffer->read((BYTE*)blockbuff, 16);

            // CBC xor for encryption
            if (cbc && !decrypt){
                EncTools::xorIv((BYTE*)blockbuff, (BYTE*)prevBlock);
            }

            arr_to_W128b(blockbuff, 0, state);

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
            W128b_to_arr(blockBuffOut, 0, state);
            ssize_t writeBytes = N_BYTES;

            // Decrypt CBC
            if (cbc && decrypt){
                EncTools::xorIv((BYTE*)blockBuffOut, (BYTE*)prevBlock);
                EncTools::copyBlock((BYTE*)prevBlock, (BYTE*)blockbuff);
            }
            if (cbc && !decrypt){
                EncTools::copyBlock((BYTE*)prevBlock, (BYTE*)blockBuffOut);
            }

            // If is the last decryption block with the padding - remove the padding.
            if (decrypt && padding && eof && k + 1 >= iter2comp){
                char paddingVal = blockBuffOut[N_BYTES - 1];
                if (paddingVal <= 0 || paddingVal > N_BYTES){
                    paddingOk = false;
                }

                for(int px = N_BYTES - paddingVal; paddingOk && px < N_BYTES; ++px){
                    paddingOk &= blockBuffOut[px] == paddingVal;
                }

                if (paddingOk){
                    writeBytes -= paddingVal;
                }
            }

            if (out && writeBytes > 0) {
                out->write((BYTE*)blockBuffOut, (size_t)writeBytes);
            }
        }

    } while(!eof);

    if (!paddingOk){
        throw std::invalid_argument("Padding is not OK");
    }
}
