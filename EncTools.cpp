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
    unsigned int paddingByte=0;

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

        if (eof && coding && coding->flags != WBAESGEN_EXTGEN_ID && (long int)blocks != blocks_rounded){
            throw std::invalid_argument("Error: External IO encodings used, input has to be block aligned (processed by external IO encoding)");
        }

        // PKCS5 padding bytes computation so the input is block aligned
        if (eof && padding && !decrypt && coding->flags == WBAESGEN_EXTGEN_ID){
            auto missingBytes = (blocks_rounded * N_BYTES - bufferSize);
            paddingByte = static_cast<unsigned int>(missingBytes == 0 ? N_BYTES : missingBytes);
            assert(paddingByte > 0 && paddingByte <= N_BYTES && "Padding size is invalid");
            iter2comp += missingBytes == 0 ? 1 : 0;
        }

        for(int k = 0; k < iter2comp; k++, blockCount++){
            // Read 1 AES block from the ring buffer.
            ssize_t bytesRead = buffer->read((BYTE*)blockbuff, N_BYTES);

            // Pad to 16 bytes before IO encodings - may happen for the last block
            if (bytesRead >= 0 && bytesRead != N_BYTES){
                memset(blockbuff + bytesRead, 0, static_cast<size_t>(N_BYTES - bytesRead));
            }

            // strip IO encodings before processing further.
            if (coding && generator && coding->flags != WBAESGEN_EXTGEN_ID) {
                generator->applyExternalEnc((BYTE*)blockbuff, coding, true);
            }

            // Encryption padding
            if (eof && padding && !decrypt && k + 1 >= iter2comp) {
                memset(blockbuff + 16 - paddingByte, (char) paddingByte, (size_t) paddingByte);
            }

            // CBC xor for encryption
            if (cbc && !decrypt){
                EncTools::xorIv((BYTE*)blockbuff, (BYTE*)prevBlock);
            }

            // Timing - core cipher
            if (cacc) {
                time(&cstart);
            }

            if(pacc) {
                pstart = clock();
            }

            // Cipher main operation
            arr_to_W128b(blockbuff, 0, state);
            if (decrypt){
                wbaes->decrypt(state);
            } else {
                wbaes->encrypt(state);
            }
            W128b_to_arr(blockBuffOut, 0, state);

            // Timing - core cipher
            if (pacc) {
                pend = clock();
                *pacc += (pend - pstart);
            }

            if (cacc) {
                time(&cend);
                *cacc += (cend - cstart);
            }

            // result of the cipher operation
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
            // If IO encodings are used padding cannot be stripped
            if (decrypt && padding && eof && k + 1 >= iter2comp && coding->flags == WBAESGEN_EXTGEN_ID){
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

            // add IO encodings before output.
            if (coding && generator && coding->flags != WBAESGEN_EXTGEN_ID) {
                generator->applyExternalEnc((BYTE*)blockBuffOut, coding, false);
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
