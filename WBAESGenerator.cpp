/*
 * WBAESGenerator.cpp
 *
 *  Created on: Mar 10, 2013
 *      Author: ph4r05
 */

#include "WBAESGenerator.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

// Enable NTL library here
NTL_CLIENT

using namespace std;
using namespace NTL;

int WBAESGenerator::shiftRowsLBijection[N_BYTES] = {
		0, 13, 10, 7,
		4,  1, 14, 11,
		8,  5,  2, 15,
		12, 9,  6,  3
};

// Shift rows selector
int WBAESGenerator::shiftRows[N_BYTES] = {
		 0,  5, 10, 15,
		 4,  9, 14,  3,
		 8, 13,  2,  7,
		12,  1,  6, 11
};

WBAESGenerator::WBAESGenerator() {
	// TODO Auto-generated constructor stub

}

WBAESGenerator::~WBAESGenerator() {
	// TODO Auto-generated destructor stub
}



void WBAESGenerator::encGenerateCodingMap(WBACR_AES_CODING_MAP& map, int *codingCount){
	int r,i,cIdx=0;
	for(r=0; r<N_ROUNDS; r++){
		// iterate over strips/MC cols
		for(i=0; i<4; i++){
			//
			// Allocation part, OUTPUT direction creates/defines new mapping
			//
			// Allocate new coding for T2, boxes, output direction (input is always set by others)
			ALLOCW08x32Coding(map.eT2[r][i][0], cIdx);
			ALLOCW08x32Coding(map.eT2[r][i][1], cIdx);
			ALLOCW08x32Coding(map.eT2[r][i][2], cIdx);
			ALLOCW08x32Coding(map.eT2[r][i][3], cIdx);

			// Allocate new coding for T3, boxes, output direction (input is always set by others)
			ALLOCW08x32Coding(map.eT3[r][i][0], cIdx);
			ALLOCW08x32Coding(map.eT3[r][i][1], cIdx);
			ALLOCW08x32Coding(map.eT3[r][i][2], cIdx);
			ALLOCW08x32Coding(map.eT3[r][i][3], cIdx);

			// Allocate new coding for XOR boxes, layer 1,2
			ALLOCXORCoding(map.eXOR1[r][i], 0, cIdx);
			ALLOCXORCoding(map.eXOR1[r][i], 8, cIdx);
			ALLOCXORCoding(map.eXOR1[r][i], 16, cIdx);

			// Allocate new coding for XOR boxes, layer 3,4
			ALLOCXORCoding(map.eXOR2[r][i], 0, cIdx);
			ALLOCXORCoding(map.eXOR2[r][i], 8, cIdx);
			ALLOCXORCoding(map.eXOR2[r][i], 16, cIdx);

			//
			// Connecting part - connecting allocated codings together
			//

			// Connect T2 boxes to XOR input boxes
			CONNECT_W08x32_TO_XOR_H(map.eT2[r][i][0], map.eXOR1[r][i], 0);	// HIGH part of XOR1_0 table IN connect to OUT of T2_0 table
			CONNECT_W08x32_TO_XOR_L(map.eT2[r][i][1], map.eXOR1[r][i], 0);  // LOW  part of XOR1_0 table IN connect to OUT of T2_1 table
			CONNECT_W08x32_TO_XOR_H(map.eT2[r][i][2], map.eXOR1[r][i], 8);  // HIGH part of XOR1_1 table IN connect to OUT of T2_2 table
			CONNECT_W08x32_TO_XOR_L(map.eT2[r][i][3], map.eXOR1[r][i], 8);  // LOW  part of XOR1_1 table IN connect to OUT of T2_3 table

			// Connect XOR layer 1 to XOR layer 2
			CONNECT_XOR_TO_XOR_H(map.eXOR1[r][i], 0, map.eXOR1[r][i], 16);  // HIGH part of XOR1_2 table IN connect to OUT of XOR1_0
			CONNECT_XOR_TO_XOR_L(map.eXOR1[r][i], 8, map.eXOR1[r][i], 16);  // LOW  part of XOR1_2 table IN connect to OUT of XOR1_1

			// Connect result XOR layer 2 to B boxes (T3)
			CONNECT_XOR_TO_W08x32(map.eXOR1[r][i], 16, map.eT3[r][i][0]);
			CONNECT_XOR_TO_W08x32(map.eXOR1[r][i], 18, map.eT3[r][i][1]);
			CONNECT_XOR_TO_W08x32(map.eXOR1[r][i], 20, map.eT3[r][i][2]);
			CONNECT_XOR_TO_W08x32(map.eXOR1[r][i], 22, map.eT3[r][i][3]);

			// Connect B boxes to XOR
			CONNECT_W08x32_TO_XOR_H(map.eT3[r][i][0], map.eXOR2[r][i], 0);	// HIGH part of XOR2_0 table IN connect to OUT of T3_0 table
			CONNECT_W08x32_TO_XOR_L(map.eT3[r][i][1], map.eXOR2[r][i], 0);  // LOW  part of XOR2_0 table IN connect to OUT of T3_1 table
			CONNECT_W08x32_TO_XOR_H(map.eT3[r][i][2], map.eXOR2[r][i], 8);  // HIGH part of XOR2_1 table IN connect to OUT of T3_2 table
			CONNECT_W08x32_TO_XOR_L(map.eT3[r][i][3], map.eXOR2[r][i], 8);  // LOW  part of XOR2_1 table IN connect to OUT of T3_3 table

			// Connect XOR layer 3 to XOR layer 4
			CONNECT_XOR_TO_XOR_H(map.eXOR2[r][i], 0, map.eXOR2[r][i], 16);  // HIGH part of XOR1_2 table IN connect to OUT of XOR1_0
			CONNECT_XOR_TO_XOR_L(map.eXOR2[r][i], 8, map.eXOR2[r][i], 16);  // LOW  part of XOR1_2 table IN connect to OUT of XOR1_1

			// in last round - we are done with 4x4 bijections
			if (r==N_ROUNDS-1) break;
			int newIdx;

			// Connect result XOR layer 4 to T2 boxes in next round unless this is the last one
			newIdx = WBAESGenerator::shiftRowsLBijection[4*i+0];
			CONNECT_XOR_TO_W08x32(map.eXOR2[r][i], 16, map.eT2[r+1][ newIdx % 4 ][ newIdx / 4]);
			newIdx = WBAESGenerator::shiftRowsLBijection[4*i+1];
			CONNECT_XOR_TO_W08x32(map.eXOR2[r][i], 18, map.eT2[r+1][ newIdx % 4 ][ newIdx / 4]);
			newIdx = WBAESGenerator::shiftRowsLBijection[4*i+2];
			CONNECT_XOR_TO_W08x32(map.eXOR2[r][i], 20, map.eT2[r+1][ newIdx % 4 ][ newIdx / 4]);
			newIdx = WBAESGenerator::shiftRowsLBijection[4*i+3];
			CONNECT_XOR_TO_W08x32(map.eXOR2[r][i], 22, map.eT2[r+1][ newIdx % 4 ][ newIdx / 4]);
		}
	}

	*codingCount = cIdx+1;
}

int WBAESGenerator::generateMixingBijections(MB08x08_TABLE ** L08x08[MB_CNT_08x08_PER_ROUND], int L08x08rounds, MB32x32_TABLE ** MB32x32[MB_CNT_32x32_PER_ROUND], int MB32x32rounds){
	int r,i;

	// Generate all required 8x8 mixing bijections.
	for(r=0; r<L08x08rounds; r++){
		for(i=0; i<MB_CNT_08x08_PER_ROUND; i++){
			generateMixingBijection(L08x08[r][i]->mb, 8, 4);
			L08x08[r][i]->inv = inv(L08x08[r][i]->mb);
		}
	}

	// Generate all required 32x32 mixing bijections.
	for(r=0; r<MB32x32rounds; r++){
		for(i=0; i<MB_CNT_32x32_PER_ROUND; i++){
			generateMixingBijection(MB32x32[r][i]->mb, 32, 4);
			MB32x32[r][i]->inv = inv(MB32x32[r][i]->mb);
		}
	}
	return 0;
}

int WBAESGenerator::generateMixingBijections(){
	return generateMixingBijections((MB08x08_TABLE***) &this->MB_L08x08, MB_CNT_08x08_ROUNDS, (MB32x32_TABLE***) &this->MB_MB32x32, MB_CNT_32x32_ROUNDS);
}

void WBAESGenerator::encGenerateTables(BYTE *key, enum keySize ksize){
	WBACR_AES_CODING_MAP 		codingMap;
	int 						codingCount;
	int							i,j,r,b,k;
	WBAES						genAES;

	// Initialize IO coding map (networked fashion of mappings)
	encGenerateCodingMap(codingMap, &codingCount);

	// Preparing all 4Bits internal encoding/decoding bijections
	CODING4X4_TABLE* pCoding04x04 = new CODING4X4_TABLE[codingCount+1];
	this->generate4X4Bijections(pCoding04x04, codingCount+1);

	// Input/Output coding
	CODING8X8_TABLE* pCoding08x08 = new CODING8X8_TABLE[16];
	this->generate8X8Bijections(pCoding08x08, 16);

	// Generate mixing bijections
	MB08x08_TABLE eMB_L08x08 [MB_CNT_08x08_ROUNDS][MB_CNT_08x08_PER_ROUND];
	MB32x32_TABLE eMB_MB32x32[MB_CNT_32x32_ROUNDS][MB_CNT_32x32_PER_ROUND];
	this->generateMixingBijections((MB08x08_TABLE***) &eMB_L08x08, MB_CNT_08x08_ROUNDS, (MB32x32_TABLE***) &eMB_MB32x32, MB_CNT_32x32_ROUNDS);

	// Generate random dual AES instances, key schedule
	vec_GF2E vecRoundKey[N_ROUNDS][N_SECTIONS];
	vec_GF2E vecKey[N_ROUNDS][N_SECTIONS];
	for(i=0; i<N_ROUNDS * N_SECTIONS; i++){
		int rndPolynomial = 0;//rand() % AES_IRRED_POLYNOMIALS;
		int rndGenerator = 0;//rand() % AES_GENERATORS;
		this->AESCipher[i].init(rndPolynomial, rndGenerator);

		// convert BYTE[] to key
		vecKey[i/N_SECTIONS][i%N_SECTIONS].SetLength(ksize);
		for(j=0; j<ksize; j++){
			vecKey[i/N_SECTIONS][i%N_SECTIONS].put(j, GF2EFromLong(key[j], 8));
		}

		// Prepare key schedule from vector representation of encryption key
		this->AESCipher[i].expandKey(vecRoundKey[i/N_SECTIONS][i%N_SECTIONS], vecKey[i/N_SECTIONS][i%N_SECTIONS], ksize);
	}

	// Connect 8x8 mapping to input/output
	// TODO: finish this...

	// 0..9 rounds, include MixColumns
	for(r=0; r<N_ROUNDS-1; r++){
		// Iterate by mix cols/sections/dual AES-es
		for(i=0; i<=N_SECTIONS; i++){
			// Restore modulus for current AES for computation in GF2E.
			this->AESCipher[i].restoreModulus();

			//
			// Build L lookup table from L_k stripes using shiftRowsLBijection (Lr_k is just simplification for indexes)
			// Now we are determining Lbox that will be used in next round.
			// Also pre-compute lookup tables by matrix multiplication
			//
			mat_GF2 * Lr_k[4];
			BYTE	  Lr_k_table[4][256];
			for(i=0; i<=N_SECTIONS; i++){
				for(j=0; j<N_SECTIONS; j++){
					Lr_k[j] = &(eMB_L08x08[r][this->shiftRowsLBijection[i*N_SECTIONS + j]].mb);
					for(b=0; b<256; b++){
						mat_GF2 tmpMat(INIT_SIZE, 8, 1);
						BYTE_to_matGF2(b, tmpMat, 0, 0);
						// multiply with 8x8 mixing bijection to obtain transformed value
						tmpMat = *(Lr_k[j]) * tmpMat;
						// convert back to byte value
						Lr_k_table[j][b] = matGF2_to_BYTE(tmpMat,0,0);
					}
				}
			}

			//
			// T table construction (Type2)
			//
			for(j=0; j<N_SECTIONS; j++){
				// Build tables - for each byte
				for(b=0; b<256; b++){
					W32b 		mapResult;
					mat_GF2 	mPreMB;
					mat_GF2E 	mcres;
					int	 		bb = b;

					// Decode input with IO coding
					bb = iocoding_encode08x08(b, codingMap.eT2[r][i][j].IC, true, pCoding04x04, pCoding08x08);
					// Mixing bijection - matrix multiplication in GF2, but only in case we are not in first round.
					// This uses inversion of transformation used in previous round in T3 box, so in first round
					// no L multiplication is done.
					if(r>0){
						mat_GF2 tmpMat(INIT_SIZE, 8, 1);
						BYTE_to_matGF2(bb, tmpMat, 0, 0);
						tmpMat = eMB_L08x08[r-1][i*N_SECTIONS].inv * tmpMat;
						bb = matGF2_to_BYTE(tmpMat, 0, 0);
					}

					// TODO: build T_i box by composing with round key

					// SBox transformation with dedicated AES for this round and section
					bb = this->AESCipher[r*4 + i].ByteSub(bb);
					GF2E tmpE = GF2EFromLong(bb, 8);

					//
					// MixColumn, Mixing bijection part
					//

					// Build [0 tmpE 0 0]^T stripe where tmpE is in j-th position
					mat_GF2E zj(INIT_SIZE, 4, 1);
					zj.put(j,0, tmpE);
					// Multiply with MC matrix from our AES dedicated for this round.
					mcres = this->AESCipher[r*4 + i].mixColMat * zj;
					// Apply 32x32 Mixing bijection, mPreMB is initialized to mat_GF2 with 32x1 dimensions,
					// GF2E values are encoded to binary column vectors
					mat_GF2E_to_mat_GF2_col(mPreMB, mcres, AES_FIELD_DIM);
					mPreMB = eMB_MB32x32[r][i].mb * mPreMB;
					// Convert transformed vector back to values
					matGF2_to_W32b(mPreMB, 0, 0, mapResult);
					// Encode mapResult with out encoding
					iocoding_encode32x32(mapResult, mapResult, codingMap.eT2[r][i][j], false, pCoding04x04, pCoding08x08);
					// Store result value to lookup table
					genAES.eTab2[r][i*4+j][b] = mapResult;
				}
			}


			//
			// B table construction (Type3) - just mixing bijections and L strip
			//
			for(j=0; j<N_SECTIONS; j++){
				// Build tables - for each byte
				for(b=0; b<256; b++){
					W32b mapResult;
					int	 bb = b;
					// Decode with IO encoding
					bb = iocoding_encode08x08(b, codingMap.eT3[r][i][j].IC, true, pCoding04x04, pCoding08x08);
					// Transform bb to matrix, to perform mixing bijection operation (matrix multiplication)
					mat_GF2 tmpMat(INIT_SIZE, 32, 1);
					// builds binary matrix [0 0 bb 0], if j==2
					BYTE_to_matGF2(bb, tmpMat, j*8, 0);
					// Build MB multiplication result
					tmpMat = eMB_MB32x32[r][j].inv * tmpMat;
					// Encode using L mixing bijection (another matrix multiplication)
					// Map bytes from result via L bijections
					mapResult.B[0] = Lr_k_table[0][matGF2_to_BYTE(tmpMat, 8*0, 0)];
					mapResult.B[1] = Lr_k_table[1][matGF2_to_BYTE(tmpMat, 8*1, 0)];
					mapResult.B[2] = Lr_k_table[2][matGF2_to_BYTE(tmpMat, 8*2, 0)];
					mapResult.B[3] = Lr_k_table[3][matGF2_to_BYTE(tmpMat, 8*3, 0)];
					// Encode mapResult with out encoding
					iocoding_encode32x32(mapResult, mapResult, codingMap.eT3[r][i][j], false, pCoding04x04, pCoding08x08);
					// Store result value to lookup table
					genAES.eTab3[r][i*4+j][b] = mapResult;
				}
			}

			//
			// XOR boxes
			//
			for(j=0; j<6; j++){
				// every master XOR table consists of 8 small XOR tables
				for(k=0; k<8; k++){
					// Build tables - for each byte
					for(b=0; b<256; b++){
						int	 bb = b;
						CODING * xorCoding = j > 2 ? &(codingMap.eXOR2[r][i][(j-3)*8+k]) : &(codingMap.eXOR1[r][i][j*8+k]);
						bb = iocoding_encode08x08(b, xorCoding->IC, true, pCoding04x04, pCoding08x08);
						bb = HI(bb) ^ LO(bb);
						bb =  iocoding_encode08x08(b, xorCoding->OC, false, pCoding04x04, pCoding08x08);
						//
						//            ________________________ ROUND
						//           |   _____________________ Section with same AES structure/MixCol stripe
						//           |  |   __________________ Master XOR table in section (6 in total, 3 up, 3 down)
						//           |  |  |   _______________ Slave XOR table in master table, 8 in total
						//           |  |  |  |   ____________ Mapping for input value
						//           |  |  |  |  |
						genAES.eXTab[r][i][j][k][b] = bb;
					}
				}
			}
		}
	}


	delete[] pCoding04x04;
}

int WBAESGenerator::generate4X4Bijections(CODING4X4_TABLE * tbl, size_t size){
	unsigned int i,c;
	for(i=0; i<size; i++){
		c |= generate4X4Bijection(&tbl[i].coding, &tbl[i].invCoding);
	}

	return c;
}

int WBAESGenerator::generate8X8Bijections(CODING8X8_TABLE * tbl, size_t size){
	unsigned int i,c;
	for(i=0; i<size; i++){
		c |= generate8X8Bijection(&tbl[i].coding, &tbl[i].invCoding);
	}

	return c;
}

int WBAESGenerator::generate4X4Bijection(BIJECT4X4 *biject, BIJECT4X4 *invBiject){
	return generateRandomBijection((unsigned char*)biject, (unsigned char*)invBiject, 16, 1);
}

int WBAESGenerator::generate8X8Bijection(BIJECT8X8 *biject, BIJECT8X8 *invBiject){
	return generateRandomBijection((unsigned char*)biject, (unsigned char*)invBiject, 256, 1);
}

