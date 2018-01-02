/*
 * WBAESGenerator.cpp
 *
 *  Created on: Mar 10, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
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

int WBAESGenerator::shiftRowsLBijectionInv[N_BYTES] = {
		 0,  5, 10, 15,
		 4,  9, 14,  3,
		 8, 13,  2,  7,
		12,  1,  6, 11
};

// Shift rows selector
int WBAESGenerator::shiftRows[N_BYTES] = {
		 0,  1,  2,  3,
		 5,  6,  7,  4,
		10, 11,  8,  9,
		15, 12, 13, 14
};

int WBAESGenerator::shiftRowsInv[N_BYTES] = {
		0,  1,  2,  3,
		 7,  4,  5,  6,
		10, 11,  8,  9,
		13, 14, 15, 12
};

WBAESGenerator::WBAESGenerator() {
	useDualAESARelationsIdentity=false;
	useDualAESIdentity=false;
	useDualAESSimpeAlternate=false;
	useIO04x04Identity=false;
	useIO08x08Identity=true;
	useMB08x08Identity=false;
	useMB32x32Identity=false;
	pCoding04x04 = NULL;
	pCoding08x08 = NULL;
	codingMap    = NULL;
}

WBAESGenerator::~WBAESGenerator() {
	;
}

void WBAESGenerator::generateCodingMap(WBACR_AES_CODING_MAP* map, int *codingCount, bool encrypt){
	int r,i,cIdx=0;

	assert(map!=NULL);
	assert(codingCount!=NULL);

	// Encryption/Decryption dependent operation and tables
	int (&shiftOp)[N_BYTES]                           = encrypt ? (this->shiftRowsLBijection) : (this->shiftRowsLBijectionInv);
	W08x128Coding (&edT1)[2][N_BYTES]                 = encrypt ? (map->eT1)                  : (map->dT1);
	W08x32Coding  (&edT2)[N_ROUNDS][N_SECTIONS][4]    = encrypt ? (map->eT2)                  : (map->dT2);
	W08x32Coding  (&edT3)[N_ROUNDS][N_SECTIONS][4]    = encrypt ? (map->eT3)                  : (map->dT3);
	CODING        (&edXOR1)[N_ROUNDS][N_SECTIONS][24] = encrypt ? (map->eXOR1)                : (map->dXOR1);
	CODING        (&edXOR2)[N_ROUNDS][N_SECTIONS][24] = encrypt ? (map->eXOR2)                : (map->dXOR2);
	CODING        (&edXOR3)[2][XTB_CNT_T1]            = encrypt ? (map->eXOR3)                : (map->dXOR3);

	// At first allocate new bijections for T1 output tables
	// Allocate encodings for XOR cascade summing output of T1 boxes
	for(r=0;r<2;r++){
		// T1 out
		for(i=0; i<N_BYTES; i++){
			ALLOCW08x128Coding(edT1[r][i], cIdx);
		}

		// XOR table cascade for T1 out sum, 8,4,2,1 = 15 XOR tables
		// Caution! Last 128-bit XOR table from T1[1] is output from whole cipher -> no allocation for this
		for(i=0; i<XTB_CNT_T1; i+=32){
			if (r==1 && i==(XTB_CNT_T1-32)) continue; 	// not for output XOR table

			ALLOCXOR128Coding(edXOR3[r], i, cIdx);
		}

		// Now connect output of T1 table (edT1[0]) for first round of XOR cascade
		// then connect whole XOR cascade together
		for(i=0; i<N_BYTES; i+=2){
			// first numerical argument is XOR table offset
			// = (i/2 * 32); offset in global XOR3 table. It takes 32 XOR tables to SUM 2xT1.
			// Thus in first level we have 8 pairs of T1, we need 8*32 = 256 XOR tables
			int xtbId = i*16;

			CONNECT_W08x32_TO_XOR_H_EX(edT1[r][i+0], edXOR3[r], xtbId+0,  0);  // HIGH part of XOR3_0 table IN connect to OUT of T1_0 table,  0..31 bit
			CONNECT_W08x32_TO_XOR_L_EX(edT1[r][i+1], edXOR3[r], xtbId+0,  0);  // LOW  part of XOR3_0 table IN connect to OUT of T1_1 table,  0..31 bit

			CONNECT_W08x32_TO_XOR_H_EX(edT1[r][i+0], edXOR3[r], xtbId+8,  4);  // HIGH part of XOR3_0 table IN connect to OUT of T1_0 table, 32..63 bit
			CONNECT_W08x32_TO_XOR_L_EX(edT1[r][i+1], edXOR3[r], xtbId+8,  4);  // LOW  part of XOR3_0 table IN connect to OUT of T1_1 table, 32..63 bit

			CONNECT_W08x32_TO_XOR_H_EX(edT1[r][i+0], edXOR3[r], xtbId+16,  8); // HIGH part of XOR3_0 table IN connect to OUT of T1_0 table, 64..95 bit
			CONNECT_W08x32_TO_XOR_L_EX(edT1[r][i+1], edXOR3[r], xtbId+16,  8); // LOW  part of XOR3_0 table IN connect to OUT of T1_1 table, 64..95 bit

			CONNECT_W08x32_TO_XOR_H_EX(edT1[r][i+0], edXOR3[r], xtbId+24, 12); // HIGH part of XOR3_0 table IN connect to OUT of T1_0 table, 96..127 bit
			CONNECT_W08x32_TO_XOR_L_EX(edT1[r][i+1], edXOR3[r], xtbId+24, 12); // LOW  part of XOR3_0 table IN connect to OUT of T1_1 table, 96..127 bit
		}

		// Now connect XOR cascade 8->4, 4->2, 2->1
		//
		// They are starting at offsets           00   32   64   96  128   160    192    224
		// We have 8*32 XOR tables, they sum T1: [01] [23] [45] [67] [89] [1011] [1213] [1415]
		// The task is to connect                  \   /     \   /     \   /        \    /
		//                                         [0123]    [4567]   [891011]    [12131415]
		// On indexes                               256       288       320          352
		//                                             \       /          \          /
		//                                              \     /            \        /
		//                                             [01234567]       [89101112131415]
		// On indexes                                     384                 416
		//                                                  \                 /
		//                                                   \               /
        //                                                [0123456789101112131415]
		// On index                                                 448
		//
		for(i=0; i<8; i+=2){
			// index of XOR tables we can use on 2. level
			int xtbId = i*16+256;
			CONNECT_XOR_TO_XOR_128_H(edXOR3[r], (i+0)*32, edXOR3[r], xtbId);
			CONNECT_XOR_TO_XOR_128_L(edXOR3[r], (i+1)*32, edXOR3[r], xtbId);
		}

		// 4->2
		CONNECT_XOR_TO_XOR_128_H(edXOR3[r], 256, edXOR3[r], 384);
		CONNECT_XOR_TO_XOR_128_L(edXOR3[r], 288, edXOR3[r], 384);
		CONNECT_XOR_TO_XOR_128_H(edXOR3[r], 320, edXOR3[r], 416);
		CONNECT_XOR_TO_XOR_128_L(edXOR3[r], 352, edXOR3[r], 416);
		// 2->1
		CONNECT_XOR_TO_XOR_128_H(edXOR3[r], 384, edXOR3[r], 448);
		CONNECT_XOR_TO_XOR_128_L(edXOR3[r], 416, edXOR3[r], 448);
	}

	// Now connect XOR3 tables form R=0 (sums T1 input table) to input of T2 tables
	// Result is stored in last XOR table starting on 448 offset, result is stored in LOW value
	// Note that ShiftRows is done here, every Sbox uses result of ShiftRows operation on its input
	//
	// 128-bit XOR has output indexed by rows, same as state.
	//
	// Connects last XOR:
	// 00 01 02 03 | 04 05 06 07 | 08 09 10 11 | 12 13 14 15  -- classical numbering (according to enc. routine)
	// 00 04 08 12 | 13 01 05 09 | 10 14 06 02 | 07 11 15 03  -- T2 boxes, indexed by column first (processed by cols in enc routine)
	//
	for(i=0; i<N_BYTES; i++){
		int newIdx = shiftOp[idxTranspose(i)];
		CONNECT_XOR_TO_W08x32(edXOR3[0], 448+i*2, edT2[0][newIdx/4][newIdx%4]);	// veryfied, OK
	}

	//
	// In the last round there is only T1 table, with defined output mapping by user (external)
	// so it is not allocated here. There are no XOR tables and T3 tables in 10. round.
	//
	// Thus encode only round 1..9.
	// Last round 9 output coding from XOR2 master table
	// is connected to T1[1] input coding in round 10.
	for(r=0; r<(N_ROUNDS-1); r++){
		// iterate over strips/MC cols
		for(i=0; i<4; i++){
			//
			// Allocation part, OUTPUT direction creates/defines new mapping
			//
			// Allocate new coding for T2, boxes, output direction (input is always set by others)
			ALLOCW08x32Coding(edT2[r][i][0], cIdx);
			ALLOCW08x32Coding(edT2[r][i][1], cIdx);
			ALLOCW08x32Coding(edT2[r][i][2], cIdx);
			ALLOCW08x32Coding(edT2[r][i][3], cIdx);

			// Allocate new coding for T3, boxes, output direction (input is always set by others)
			ALLOCW08x32Coding(edT3[r][i][0], cIdx);
			ALLOCW08x32Coding(edT3[r][i][1], cIdx);
			ALLOCW08x32Coding(edT3[r][i][2], cIdx);
			ALLOCW08x32Coding(edT3[r][i][3], cIdx);

			// Allocate new coding for XOR boxes, layer 1,2
			ALLOCXORCoding(edXOR1[r][i], 0, cIdx);
			ALLOCXORCoding(edXOR1[r][i], 8, cIdx);
			ALLOCXORCoding(edXOR1[r][i], 16, cIdx);

			// Allocate new coding for XOR boxes, layer 3,4
			ALLOCXORCoding(edXOR2[r][i], 0, cIdx);
			ALLOCXORCoding(edXOR2[r][i], 8, cIdx);
			ALLOCXORCoding(edXOR2[r][i], 16, cIdx);

			//
			// Connecting part - connecting allocated codings together
			//

			// Connect T2 boxes to XOR input boxes
			CONNECT_W08x32_TO_XOR_H(edT2[r][i][0], edXOR1[r][i], 0);  // HIGH part of XOR1_0 table IN connect to OUT of T2_0 table
			CONNECT_W08x32_TO_XOR_L(edT2[r][i][1], edXOR1[r][i], 0);  // LOW  part of XOR1_0 table IN connect to OUT of T2_1 table
			CONNECT_W08x32_TO_XOR_H(edT2[r][i][2], edXOR1[r][i], 8);  // HIGH part of XOR1_1 table IN connect to OUT of T2_2 table
			CONNECT_W08x32_TO_XOR_L(edT2[r][i][3], edXOR1[r][i], 8);  // LOW  part of XOR1_1 table IN connect to OUT of T2_3 table

			// Connect XOR layer 1 to XOR layer 2
			CONNECT_XOR_TO_XOR_H(edXOR1[r][i], 0, edXOR1[r][i], 16);  // HIGH part of XOR1_2 table IN connect to OUT of XOR1_0
			CONNECT_XOR_TO_XOR_L(edXOR1[r][i], 8, edXOR1[r][i], 16);  // LOW  part of XOR1_2 table IN connect to OUT of XOR1_1

			// Connect result XOR layer 2 to B boxes (T3)
			CONNECT_XOR_TO_W08x32(edXOR1[r][i], 16, edT3[r][i][0]);
			CONNECT_XOR_TO_W08x32(edXOR1[r][i], 18, edT3[r][i][1]);
			CONNECT_XOR_TO_W08x32(edXOR1[r][i], 20, edT3[r][i][2]);
			CONNECT_XOR_TO_W08x32(edXOR1[r][i], 22, edT3[r][i][3]);

			// Connect B boxes to XOR
			CONNECT_W08x32_TO_XOR_H(edT3[r][i][0], edXOR2[r][i], 0);	// HIGH part of XOR2_0 table IN connect to OUT of T3_0 table
			CONNECT_W08x32_TO_XOR_L(edT3[r][i][1], edXOR2[r][i], 0);  // LOW  part of XOR2_0 table IN connect to OUT of T3_1 table
			CONNECT_W08x32_TO_XOR_H(edT3[r][i][2], edXOR2[r][i], 8);  // HIGH part of XOR2_1 table IN connect to OUT of T3_2 table
			CONNECT_W08x32_TO_XOR_L(edT3[r][i][3], edXOR2[r][i], 8);  // LOW  part of XOR2_1 table IN connect to OUT of T3_3 table

			// Connect XOR layer 3 to XOR layer 4
			CONNECT_XOR_TO_XOR_H(edXOR2[r][i], 0, edXOR2[r][i], 16);  // HIGH part of XOR1_2 table IN connect to OUT of XOR1_0
			CONNECT_XOR_TO_XOR_L(edXOR2[r][i], 8, edXOR2[r][i], 16);  // LOW  part of XOR1_2 table IN connect to OUT of XOR1_1

			int newIdx;

			if (r<(N_ROUNDS-2)){
				// Connect result XOR layer 4 to T2 boxes in next round
				newIdx = shiftOp[4*i+0]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 16, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
				newIdx = shiftOp[4*i+1]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 18, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
				newIdx = shiftOp[4*i+2]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 20, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
				newIdx = shiftOp[4*i+3]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 22, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
			} else {
				// Connect result XOR layer 4 to T1 boxes in last round; r==8
				newIdx = shiftOp[4*i+0]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 16, edT1[1][newIdx]);
				newIdx = shiftOp[4*i+1]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 18, edT1[1][newIdx]);
				newIdx = shiftOp[4*i+2]; CONNECT_XOR_TO_W08x32(edXOR2[r][i], 20, edT1[1][newIdx]);
				newIdx = shiftOp[4*i+3];CONNECT_XOR_TO_W08x32(edXOR2[r][i], 22, edT1[1][newIdx]);
			}
		}
	}

	*codingCount = cIdx+1;
}

int WBAESGenerator::generateMixingBijections(
		MB08x08_TABLE L08x08[][MB_CNT_08x08_PER_ROUND], int L08x08rounds,
		MB32x32_TABLE MB32x32[][MB_CNT_32x32_PER_ROUND], int MB32x32rounds,
		bool MB08x08Identity, bool MB32x32Identity){
	int r,i;

	// Generate all required 8x8 mixing bijections.
	for(r=0; r<L08x08rounds; r++){
		for(i=0; i<MB_CNT_08x08_PER_ROUND; i++){
			if (!MB08x08Identity){
				generateMixingBijection(L08x08[r][i].mb, 8, 4);
				L08x08[r][i].inv = inv(L08x08[r][i].mb);
			} else {
				ident(L08x08[r][i].mb, 8);
				ident(L08x08[r][i].inv, 8);
			}
		}
	}

	// Generate all required 32x32 mixing bijections.
	for(r=0; r<MB32x32rounds; r++){
		for(i=0; i<MB_CNT_32x32_PER_ROUND; i++){
			if (!MB32x32Identity){
				generateMixingBijection(MB32x32[r][i].mb, 32, 4);
				MB32x32[r][i].inv = inv(MB32x32[r][i].mb);
			} else {
				ident(MB32x32[r][i].mb, 32);
				ident(MB32x32[r][i].inv, 32);
			}
		}
	}
	return 0;
}

int WBAESGenerator::generateMixingBijections(bool identity){
	return generateMixingBijections(this->MB_L08x08, MB_CNT_08x08_ROUNDS, this->MB_MB32x32, MB_CNT_32x32_ROUNDS, identity);
}

void WBAESGenerator::generateExtEncoding(ExtEncoding * extc, int flags){
	int k;

	// generate 8x8 bijections at first
	for(k=0; k<2; k++){
		bool identity = (k==0 && (flags & WBAESGEN_EXTGEN_fCID) > 0) || (k==1 && (flags & WBAESGEN_EXTGEN_lCID) > 0);
		this->generate4X4Bijections(&(extc->lfC[k][0]), 2*N_BYTES, identity);
	}

	// generate mixing bijection
	for(k=0; k<2; k++){
		bool identity = (k==0 && (flags & WBAESGEN_EXTGEN_IDMID) > 0) || (k==1 && (flags & WBAESGEN_EXTGEN_ODMID) > 0);
		if (!identity){
			generateMixingBijection(extc->IODM[k].mb, 128, 4);
			extc->IODM[k].inv = inv(extc->IODM[k].mb);
		} else {
			ident(extc->IODM[k].mb,  128);
			ident(extc->IODM[k].inv, 128);
		}
	}
	extc->flags = flags;
}

void WBAESGenerator::generateT1Tables(WBAES * genAES, ExtEncoding * extc, bool encrypt){
	// To initialize T1[1] map, coding map is needed, since it takes input from last round, for this we need key material
	// to add S-box to T1[1], so it is not done here...
	int i,j,b;

	// Encryption/Decryption dependent operation and tables
	AES_TB_TYPE1  (&genAES_edTab1)[2][N_BYTES]  = encrypt ? genAES->eTab1                : genAES->dTab1;
	W08x128Coding (&codingMap_edT1)[2][N_BYTES] = encrypt ? codingMap->eT1               : codingMap->dT1;

	// At first initialize T1[0]
	for(i=0; i<N_BYTES; i++){
		// i-th T1 table, indexed by cols

		// Build tables - for each byte
		for(b=0; b<256; b++){
			W128b mapResult;
			int	  bb = b;
			// Decode with IO encoding
			bb = HILO(extc->lfC[0][2*i+0].invCoding[HI(b)], extc->lfC[0][2*i+1].invCoding[LO(b)]);
			// Transform bb to matrix, to perform mixing bijection operation (matrix multiplication)
			mat_GF2 tmpMat(INIT_SIZE, 128, 1);
			// builds binary matrix [0 0 bb 0 0 0 0 0 0 0 0 0 0 0 0 0], if i==2
			BYTE_to_matGF2(bb, tmpMat, idxTranspose(i)*8, 0);
			// Build MB multiplication result
			tmpMat = extc->IODM[0].inv * tmpMat;
			// Encode 128-bit wide output to map result
			for(j=0; j<16; j++){
				mapResult.B[j] = matGF2_to_BYTE(tmpMat, 8*j, 0);
			}
			// Encode mapResult with out encoding of T1 table
			iocoding_encode128x128(mapResult, mapResult, codingMap_edT1[0][i], false, pCoding04x04, pCoding08x08);
			// Store result value to lookup table
			W128CP(genAES_edTab1[0][i][b], mapResult);
		}
	}
}

void WBAESGenerator::generateTables(BYTE *key, enum keySize ksize, WBAES * genAES, ExtEncoding* extc, bool encrypt){
	int 						codingCount;
	int							i,j,r,b,k;
	GenericAES					defaultAES;

	// Initialize IO coding map (networked fashion of mappings)
	this->codingMap = new WBACR_AES_CODING_MAP;
	generateCodingMap(codingMap, &codingCount, encrypt);

	// Preparing all 4Bits internal encoding/decoding bijections
	this->pCoding04x04 = new CODING4X4_TABLE[codingCount+1];
	this->generate4X4Bijections(pCoding04x04, codingCount+1, useIO04x04Identity);

	// Generate mixing bijections
	MB08x08_TABLE eMB_L08x08 [MB_CNT_08x08_ROUNDS][MB_CNT_08x08_PER_ROUND];
	MB32x32_TABLE eMB_MB32x32[MB_CNT_32x32_ROUNDS][MB_CNT_32x32_PER_ROUND];
	this->generateMixingBijections(eMB_L08x08, MB_CNT_08x08_ROUNDS, eMB_MB32x32, MB_CNT_32x32_ROUNDS, useMB08x08Identity, useMB32x32Identity);

	// Encryption/Decryption dependent functions and tables
	int (&nextTbox)[N_BYTES]     = encrypt ? (this->shiftRowsLBijection) : (this->shiftRowsLBijectionInv);
	int (&shiftRowsOp)[N_BYTES]  = encrypt ? (this->shiftRows)    		 : (this->shiftRowsInv);
	W32XTB        (&genAES_edXTab)[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS]  = encrypt ? genAES->eXTab   : genAES->dXTab;
	W32XTB        (&genAES_edXTabEx)[2][15][4]                          = encrypt ? genAES->eXTabEx : genAES->dXTabEx;
	AES_TB_TYPE1  (&genAES_edTab1)[2][N_BYTES]			 		        = encrypt ? genAES->eTab1   : genAES->dTab1;
	AES_TB_TYPE2  (&genAES_edTab2)[N_ROUNDS][N_BYTES]			 		= encrypt ? genAES->eTab2   : genAES->dTab2;
	AES_TB_TYPE3  (&genAES_edTab3)[N_ROUNDS][N_BYTES]			 		= encrypt ? genAES->eTab3   : genAES->dTab3;
	W08x128Coding (&codingMap_edT1)[2][N_BYTES]					        = encrypt ? codingMap->eT1  : codingMap->dT1;
	W08x32Coding  (&codingMap_edT2)[N_ROUNDS][N_SECTIONS][4]	        = encrypt ? codingMap->eT2  : codingMap->dT2;
	W08x32Coding  (&codingMap_edT3)[N_ROUNDS][N_SECTIONS][4]            = encrypt ? codingMap->eT3  : codingMap->dT3;
	CODING        (&codingMap_edXOR1)[N_ROUNDS][N_SECTIONS][24]         = encrypt ? codingMap->eXOR1: codingMap->dXOR1;
	CODING        (&codingMap_edXOR2)[N_ROUNDS][N_SECTIONS][24]         = encrypt ? codingMap->eXOR2: codingMap->dXOR2;
	CODING        (&codingMap_edXOR3)[2][XTB_CNT_T1]                    = encrypt ? codingMap->eXOR3: codingMap->dXOR3;

	// Init T1[0] tables - for the first round
	this->generateT1Tables(genAES, extc, encrypt);

#ifdef AES_BGE_ATTACK
	// If there are 8x8 output bijections, just generate identities
	GF256_func_t (&edOutputBijection)[N_ROUNDS][N_BYTES] = encrypt ? (genAES->eOutputBijection) : (genAES->dOutputBijection);
	for(r=0;r<N_ROUNDS;r++){
		for(i=0;i<N_BYTES;i++){
			for(k=0;k<GF256;k++){
				edOutputBijection[r][i][k]=k;
			}
		}
	}
#endif


	// A1, A2 relations
	int genA[N_ROUNDS * N_SECTIONS];						// constant a for A1A2 relations
	int genI[N_ROUNDS * N_SECTIONS];						// exponent for A1A2 relations
	vec_GF2E genA1[N_ROUNDS * N_SECTIONS];					// A1 relation
	vec_GF2E genA2[N_ROUNDS * N_SECTIONS];					// A2 relation

	// Generate random dual AES instances, key schedule
	vec_GF2E vecRoundKey[N_ROUNDS][N_SECTIONS];
	vec_GF2E vecKey[N_ROUNDS][N_SECTIONS];
	vec_GF2E defaultKey;									// key for default AES in GF2E representation
	vec_GF2E expandedKey, backupKey;						// expanded key for default AES
	defaultAES.initFromIndex(0,0);							// index 0,0 represents 0x11D and 0x03, so default AES
	BYTEArr_to_vec_GF2E(key, ksize, defaultKey);			// convert BYTE key to GF2E key
	defaultAES.expandKey(expandedKey, defaultKey, ksize);	// key schedule for default AES
	backupKey = expandedKey;								// backup default AES expanded key for test routine
	for(i=0; i<N_ROUNDS * N_SECTIONS; i++){
		int rndPolynomial = useDualAESIdentity ? 0 : phrand() % AES_IRRED_POLYNOMIALS;
		int rndGenerator  = useDualAESIdentity ? 0 : phrand() % AES_GENERATORS;
		genA[i]			  = useDualAESARelationsIdentity ? 1 : (phrand() % 255) + 1;
		genI[i]			  = useDualAESARelationsIdentity ? 0 : phrand() % 8;
		if (useDualAESSimpeAlternate && !useDualAESIdentity){
			rndPolynomial = (i)%2 == 0 ? 0: AES_IRRED_POLYNOMIALS-1;
			rndGenerator  = (i)%2 == 0 ? 0: AES_GENERATORS-1;
		}

		this->AESCipher[i].initFromIndex(rndPolynomial, rndGenerator);

		// convert BYTE[] to key
		BYTEArr_to_vec_GF2E(key, ksize, vecKey[i/N_SECTIONS][i%N_SECTIONS]);
		this->AESCipher[i].applyT(vecKey[i/N_SECTIONS][i%N_SECTIONS]);

		// Prepare key schedule from vector representation of encryption key
		this->AESCipher[i].expandKey(vecRoundKey[i/N_SECTIONS][i%N_SECTIONS], vecKey[i/N_SECTIONS][i%N_SECTIONS], ksize);

		// generate A1 A2 relations
		this->AESCipher[i].generateA1A2Relations(genA1[i], genA2[i], genA[i], genI[i]);
		if (encrypt){
			if (this->AESCipher[i].testA1A2Relations(genA1[i], genA2[i]) != 0) cout << "Error in A1A2 generator" << endl;
		} else {
			if (this->AESCipher[i].testA1A2Relations(genA2[i], genA1[i], false) != 0) cout << "Error in A1A2 generator" << endl;
		}

		if (this->AESCipher[i].testA1XorLinearity(genA1[i])!=0) cout << "Error in A1 linearity!" << endl;
	}

	//
	// Build L lookup table from L_k stripes using shiftRowsLBijection (Lr_k is just simplification for indexes)
	// Now we are determining Lbox that will be used in next round.
	// Also pre-compute lookup tables by matrix multiplication
	//
	for(r=0; r<N_ROUNDS; r++){


		// Iterate by mix cols/sections/dual AES-es
		for(i=0; i<N_SECTIONS; i++){
			// Restore modulus for current AES for computation in GF2E.
			this->AESCipher[r*4 + i].restoreModulus();

			mat_GF2 * Lr_k[4];
			BYTE	  Lr_k_table[4][256];

			for(j=0; r<(N_ROUNDS-1) && j<N_SECTIONS; j++){
				Lr_k[j] = &(eMB_L08x08[r][nextTbox[i*N_SECTIONS + j]].mb);
				for(b=0; b<256; b++){
					mat_GF2 tmpMat(INIT_SIZE, 8, 1);
					BYTE_to_matGF2(b, tmpMat, 0, 0);

					// multiply with 8x8 mixing bijection to obtain transformed value
					tmpMat = *(Lr_k[j]) * tmpMat;
					// convert back to byte value
					Lr_k_table[j][b] = matGF2_to_BYTE(tmpMat,0,0);
				}
			}

			//
			// T table construction (Type2, if r=last one, then T1)
			//
			for(j=0; j<N_SECTIONS; j++){
//cout << "T["<<r<<"]["<<i<<"]["<<j<<"] key = 16*" << r << " + " << ((int)shiftRowsOp[ j*4 + i ]) <<  " = " << (vecRoundKey[r][i][16*r + shiftRowsOp[ j*4 + i ]]) << endl;

				// Build tables - for each byte
				for(b=0; b<256; b++){
					GF2E		tmpGF2E;
					W32b 		mapResult;
					mat_GF2 	mPreMB;
					mat_GF2E 	mcres;
					int	 		bb = b;

					// In the first round we apply codings from T1 tables.
					// Decode input with IO coding
					// For the last round, INPUT coding is for T1 box, otherwise for T2 box
					if (r < (N_ROUNDS-1)){
						bb = iocoding_encode08x08(bb, codingMap_edT2[r][i][j].IC, true, pCoding04x04, pCoding08x08);
					} else {
						bb = iocoding_encode08x08(bb, codingMap_edT1[1][i*4+j].IC, true, pCoding04x04, pCoding08x08);
					}

					// Dual AES mapping.
					// Tapply(curAES, TapplyInv(prevAES, state)) for rounds > 0
					// Tapply(curAES, state) for round == 0 	(no inversion could be done, first round)
					tmpGF2E = GF2EFromLong(bb, AES_FIELD_DIM);

					//
					// Mixing bijection - removes effect induced in previous round (inversion here)
					// Note: for DualAES, data from prev round comes here in prev Dual AES encoding, with applied bijection
					// on them. Reversal = apply inverse of mixing bijection, undo prev Dual AES, do cur Dual AES
					// Scheme: Tapply_cur( TapplyInv_prev( L^{-1}_{r-1}(x) ) )
					//
					// Implementation: matrix multiplication in GF2.
					// Inversion to transformation used in previous round in T3 box (so skip this in first round).
					if(r>0){
						mat_GF2 tmpMat = colVector(tmpGF2E, AES_FIELD_DIM);
						tmpMat = eMB_L08x08[r-1][i*N_SECTIONS + j].inv * tmpMat;
						colVector(tmpGF2E, tmpMat, 0);
					}

					//
					// Encryption scenario (decryption is analogous):
					// Applying inverse transformation is a little bit tricky here, illustration follows.
					// We know that indexes to boxes T0 T1 T2 T3 from state array will take ShiftRows() into account:
					//
					// Each quartet corresponds to "section"/column (i-idx) encoded with separate Dual AES.
					// For example 0,5,10,15 state bytes are encoded with particular AES (different from quartet 01,06,11,12),
					// feed to T0,1,2,3 boxes and stored to 00,04,08,12 output state bytes.
					//
					// In next round we again take 0,5,10,15 from state array, but it previous round it was: 0,9,2,11
					// each byte encoded by different dual AES: 0,1,2,3 respectively.
					//
					// Each column in matrix below is encoded with separate dual AES. (matrix = ShiftRows(stateArray))
					//
					//             +------------------------- ShiftRows() again - next round
					//             |              +---------- In this matrix, we take columns again,
					//             |              |           so StateByte_00,06,08,14 -> Tbox_00,01,02,03.
					//             |              |
					// 00 01 02 03 | 00 01 02 03  | Encoded    0  1  2  3
					// 05 06 07 04 | 06 07 04 05  | by DUAL    1  2  3  0
					// 10 11 08 09 | 08 09 10 11  | AES from   2  3  0  1
					// 15 12 13 14 | 14 15 12 13  | prev. rnd  3  0  1  2
					//
					// Thus inverse transformation is ((4*r-1) + ((i+j) % 4)), i=section.
					// Revert Dual AES representation from previous round to default AES.
					//
					// Decryption: principle is the same, so just summary:
					//
					//             +------------------------- ShiftRowsInv() again - next round
					//             |              +---------- In this matrix, we take columns again,
					//             |              |           so StateByte_00,06,08,14 -> Tbox_00,01,02,03.
					//             |              |
					// 00 01 02 03 | 00 01 02 03  | Encoded    0  1  2  3
					// 07 04 05 06 | 06 07 04 05  | by DUAL    3  0  1  2
					// 10 11 08 09 | 08 09 10 11  | AES from   2  3  0  1
					// 13 14 15 12 | 14 15 12 13  | prev. rnd  1  2  3  0
					//
					if(r>0){
						if (encrypt)
							this->AESCipher[(4*(r-1)) + POS_MOD(i+j, 4)].applyTinv(tmpGF2E);
						else
							this->AESCipher[(4*(r-1)) + POS_MOD(i-j, 4)].applyTinv(tmpGF2E);
					}

					// Now apply new transformation - convert to cur Dual AES representation from default AES
					this->AESCipher[4*r + i].applyT(tmpGF2E);

					//
					// Encryption scenario:
					// Build T_i box by composing with round key
					//
					// White box implementation:
					// shiftRows(state)
					// addRoundKey(state, shiftRows(ApplyT(K_{r-1}))) when indexing rounds from 1 and key from 0
					//   K_{r-1} is AES key for default AES,
					//   apply = linear transformation (multiplication by matrix T from dual AES) for changing default AES to dual AES.
					//
					// Rewritten to form:
					// shiftRows(state)
					// addRoundKey(state, ApplyT(shiftRows(K_r)))
					//
					// K_{r}  [x][y] = vecRoundKey[r][i] [16*(r)   + x*4 + y]
					// in this round we want to work with AES from same dual AES, thus we are choosing
					// vecRoundKey[r][i]. Also we have to take effect of ShiftRows() into account, thus apply
					// ShiftRows() transformation on key indexes.
					//
					// Implementation in one section (i) corresponds to one column (0,5,10,15) are indexes taken
					// for computation in one section in WBAES. Inside section (column) we are iterating over
					// rows (j). Key is serialized by rows.
					if (encrypt){
						if (r==0){
							// In first round during encryption, there is no A1 applied in previous round, so just apply here
							tmpGF2E = genA1[4*r+i][getLong(tmpGF2E)];
						}

						GF2E tmpKey = vecRoundKey[r][i][16*r + idxTranspose(shiftRowsOp[ j*4 + i ])];
						tmpGF2E    += genA1[4*r+i][getLong(tmpKey)];
					} else {
						if(r==0) {
							// Decryption & first round => add k_10 to state.
							// Same logic applies here
							// AddRoundKey(State, k_10)  | -> InvShiftRows(State)
							// InvShiftRows(State)       | -> AddRoundKey(State, InvShiftRows(k_10))
							tmpGF2E += vecRoundKey[r][i][16*N_ROUNDS + idxTranspose(shiftRowsOp[ j*4 + i ])];
							tmpGF2E = genA2[4*r+i][getLong(tmpGF2E)];
						}
					}


					// SBox transformation with dedicated AES for this round and section
					// Encryption: ByteSub
					// Decryption: ByteSubInv
					GF2E tmpE = encrypt ?
							  this->AESCipher[r*4 + i].ByteSub(tmpGF2E)
							: this->AESCipher[r*4 + i].ByteSubInv(tmpGF2E);

					// Dual AES with A1 A2 relations, after Sbox apply A2 relation
					tmpE = encrypt ? genA2[4*r+i][getLong(tmpE)] : genA1[4*r+i][getLong(tmpE)];

					// Decryption case:
					// T(x) = Sbox(x) + k
					if (!encrypt){
						tmpE += vecRoundKey[r][i][16*(N_ROUNDS-r-1) + idxTranspose(j*4 + i)];
					}

					// If we are in last round we also have to add k_10, not affected by ShiftRows()
					// And more importantly, build T1
					if (r==N_ROUNDS-1){
						// Adding last encryption key (k_10) by special way is performed only in encryption
						if (encrypt) {
							tmpE += vecRoundKey[r][i][16*(r+1) + idxTranspose(j*4 + i)];
						}

						// revert last dual AES transformation here
						this->AESCipher[4*r + i].applyTinv(tmpE);

						// Now we use output encoding G and quit, no MixColumn or Mixing bijections here.
						W128b mapResult128;
						bb = getLong(tmpE);

						// Transform bb to matrix, to perform mixing bijection operation (matrix multiplication)
						mat_GF2 tmpMat2(INIT_SIZE, 128, 1);
						// builds binary matrix [0 0 bb 0 0 0 0 0 0 0 0 0 0 0 0 0], if curByte==2
						BYTE_to_matGF2(bb, tmpMat2, (i*N_SECTIONS + j)*8, 0);
						// Build MB multiplication result
						tmpMat2 = extc->IODM[1].mb * tmpMat2;
						// Encode 128-bit wide output to map result
						for(int jj=0; jj<16; jj++){
							mapResult128.B[jj] = matGF2_to_BYTE(tmpMat2, jj*8, 0);
						}
						// Encode mapResult with out encoding of T1 table
						iocoding_encode128x128(mapResult128, mapResult128, codingMap_edT1[1][(i*N_SECTIONS + j)], false, pCoding04x04, pCoding08x08);
						// Store result value to lookup table
						W128CP(genAES_edTab1[1][(i*N_SECTIONS + j)][b], mapResult128);
						continue;
					}

					//
					// MixColumn, Mixing bijection part
					//	only in case 1..9 round

					// Build [0 tmpE 0 0]^T stripe where tmpE is in j-th position
					mat_GF2E zj(INIT_SIZE, 4, 1);
					zj.put(j,0, tmpE);

					// Multiply with MC matrix from our AES dedicated for this round, only in 1..9 rounds (not in last round)
					if (encrypt){
						mcres = this->AESCipher[r*4 + i].mixColMat * zj;
					} else {
						mcres = this->AESCipher[r*4 + i].mixColInvMat * zj;
					}
					// Encryption:
					// Dual AES, apply A1 for next round here.
					//
					// We have one 4x1 GF2E matrix. In each next j-iteration we will have
					// different 4x1 matrices, XORed with each other afterwards.
					// XOR is performed by rows, thus all elements in one row are added together,
					// so they have to have same dual AES encodings => we have to apply
					// different transformation on each element.
					//
					// Every resulting element after XOR is passed to different T2 boxes.
					//
					//  Cur. round |  Next round |              |
					// 00 01 02 03 | 00 01 02 03 | AES encoding | 00 01 02 03
					// 05 06 07 04 | 06 07 04 05 | in next      | 03 00 01 02
					// 10 11 08 09 | 08 09 10 11 | round        | 02 03 00 01
					// 15 12 13 14 | 14 15 12 13 |              | 01 02 03 00
					//
					// One i iteration corresponds to one column above. One i=0 iteration should look like this:
					// Every A in next diagram is A_I = A^1_{r+1, I} - simplified syntax
					//
					// | A_0 (02 T(x)) |   | A_0 (03 T(x)) |   | A_0 (01 T(x)) |   | A_0 (01 T(x)) |
					// | A_3 (01 T(x)) | + | A_3 (02 T(x)) | + | A_3 (03 T(x)) | + | A_3 (01 T(x)) |
					// | A_2 (02 T(x)) |   | A_2 (01 T(x)) |   | A_2 (02 T(x)) |   | A_2 (03 T(x)) |
					// | A_1 (03 T(x)) |   | A_1 (01 T(x)) |   | A_1 (01 T(x)) |   | A_1 (02 T(x)) |
					//
					int tmpi;
					for(tmpi=0; tmpi<4; tmpi++){
						if (encrypt){
							this->AESCipher[ 4* r     + i                 ].applyTinv(mcres[tmpi]);
							this->AESCipher[(4*(r+1)) + POS_MOD(i-tmpi, 4)].applyT(   mcres[tmpi]);
							applyLookupTable(genA1[(4*(r+1)) + POS_MOD(i-tmpi, 4)],   mcres[tmpi]);
							this->AESCipher[(4*(r+1)) + POS_MOD(i-tmpi, 4)].applyTinv(mcres[tmpi]);
							this->AESCipher[ 4* r     + i                 ].applyT(   mcres[tmpi]);
						} else {
							this->AESCipher[ 4* r     + i                 ].applyTinv(mcres[tmpi]);
							this->AESCipher[(4*(r+1)) + POS_MOD(i+tmpi, 4)].applyT(   mcres[tmpi]);
							applyLookupTable(genA2[(4*(r+1)) + POS_MOD(i+tmpi, 4)],   mcres[tmpi]);

							//
							// Compensate affine part of A2 relation
							//
							// A2 is not linear in decryption case, but affine.
							// We have here 4 elements (entering XOR), so from 3 of them
							// we have to subtract affine constant = A2[0].
							// Af(a1+a2+a3+a4) = A*a1 + A*a2 + A*a3 + A*a4 + c
							//                 = Af(a1) + Af(a2)+Af(0) + Af(a3)+Af(0) + Af(a4)+Af(0)
							if (j!=0) {
								mcres[tmpi][0] += genA2[(4*(r+1)) + POS_MOD(i+tmpi, 4)][0];
							}

							this->AESCipher[(4*(r+1)) + POS_MOD(i+tmpi, 4)].applyTinv(mcres[tmpi]);
							this->AESCipher[ 4* r     + i                 ].applyT(   mcres[tmpi]);
						}
					}

					// Apply 32x32 Mixing bijection, mPreMB is initialized to mat_GF2 with 32x1 dimensions,
					// GF2E values are encoded to binary column vectors
					mat_GF2E_to_mat_GF2_col(mPreMB, mcres, AES_FIELD_DIM);
					mPreMB = eMB_MB32x32[r][i].mb * mPreMB;

					//
					// TESTING - multiply by inversion


					// Convert transformed vector back to values
					mapResult.l = 0;
					matGF2_to_W32b(mPreMB, 0, 0, mapResult);

					// Encode mapResult with out encoding
					iocoding_encode32x32(mapResult, mapResult, codingMap_edT2[r][i][j], false, pCoding04x04, pCoding08x08);
					// Store result value to lookup table
					genAES_edTab2[r][i*4+j][b] = mapResult;
				}
			}

			// In final round there are no more XOR and T3 boxes
			if (r==N_ROUNDS-1){
				continue;
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
					bb = iocoding_encode08x08(b, codingMap_edT3[r][i][j].IC, true, pCoding04x04, pCoding08x08);
					// Transform bb to matrix, to perform mixing bijection operation (matrix multiplication)
					mat_GF2 tmpMat(INIT_SIZE, 32, 1);
					// builds binary matrix [0 0 bb 0], if j==2
					BYTE_to_matGF2(bb, tmpMat, j*8, 0);
					// Build MB multiplication result
					tmpMat = eMB_MB32x32[r][i].inv * tmpMat;
					// Encode using L mixing bijection (another matrix multiplication)
					// Map bytes from result via L bijections
					mapResult.l = 0;
					mapResult.B[0] = Lr_k_table[0][matGF2_to_BYTE(tmpMat, 8*0, 0)];
					mapResult.B[1] = Lr_k_table[1][matGF2_to_BYTE(tmpMat, 8*1, 0)];
					mapResult.B[2] = Lr_k_table[2][matGF2_to_BYTE(tmpMat, 8*2, 0)];
					mapResult.B[3] = Lr_k_table[3][matGF2_to_BYTE(tmpMat, 8*3, 0)];
					// Encode mapResult with out encoding
					iocoding_encode32x32(mapResult, mapResult, codingMap_edT3[r][i][j], false, pCoding04x04, pCoding08x08);
					// Store result value to lookup table
					genAES_edTab3[r][i*4+j][b] = mapResult;
					// cout << "T3["<<r<<"]["<<i<<"]["<<j<<"]["<<b<<"] = "; dumpW32b(mapResult);
				}
			}

			//
			// XOR boxes
			//
			for(j=0; j<6; j++){
				// every master XOR table consists of 8 small XOR tables
				for(k=0; k<8; k++){
					CODING & xorCoding = j > 2 ? codingMap_edXOR2[r][i][(j-3)*8+k] : codingMap_edXOR1[r][i][j*8+k];
					//
					//                                            ________________________ ROUND
					//                                           |   _____________________ Section with same AES structure/MixCol stripe
					//                                           |  |   __________________ Master XOR table in section (6 in total, 3 up, 3 down)
					//                                           |  |  |   _______________ Slave XOR table in master table, 8 in total
					//                                           |  |  |  |
					generateXorTable(&xorCoding, &(genAES_edXTab[r][i][j][k]));
				}
			}
		}
	}

	//
	// XOR boxes in T1 cascade
	//
	for(r=0; r<2; r++){
		for(i=0; i<15; i++){
			for(j=0; j<4; j++){
				// every master XOR table consists of 8 small XOR tables
				for(k=0; k<8; k++){
					CODING & xorCoding = codingMap_edXOR3[r][32*i+8*j+k];
					//
					//                                              ________________________ ROUND
					//                                             |   _____________________ 0..14 8,4,2,1
					//                                             |  |   __________________ 0..4  (128-bits)
					//                                             |  |  |   _______________ 0..8  (32-bit XOR table)
					//                                             |  |  |  |
					generateXorTable(&xorCoding, &(genAES_edXTabEx[r][i][j][k]));

					//
					// Last XOR table
					//
					if (r==1 && i==14) {
						for(int b=0; b<256; b++){
							genAES_edXTabEx[r][i][j][k][b] = extc->lfC[1][8*j+k].coding[genAES_edXTabEx[r][i][j][k][b]];
						}
					}
				}
			}
		}
	}

	delete[] this->pCoding04x04;
	this->pCoding04x04 = NULL;

	delete[] this->codingMap;
	this->codingMap = NULL;
}

void WBAESGenerator::generateXorTable(CODING * xorCoding, XTB * xtb){
	for(int b=0; b<256; b++){
		int	bb = b;
		bb = iocoding_encode08x08(bb, xorCoding->IC, true, pCoding04x04, pCoding08x08);
		bb = HI(bb) ^ LO(bb);
		bb = iocoding_encode08x08(bb, xorCoding->OC, false, pCoding04x04, pCoding08x08);
		(*xtb)[b] = bb;
	}
}

int WBAESGenerator::generate4X4Bijections(CODING4X4_TABLE * tbl, size_t size, bool identity){
	unsigned long int i=0,c=0;
	for(i=0; i<size; i++){
		// HINT: if you are debugging IO problems, try to turn on and off some bijections,
		// you can very easily localize the problem.

		//if (i>=0x3c0) identity=true;
		c |= generate4X4Bijection(&tbl[i].coding, &tbl[i].invCoding, identity);
	}

	return c;
}

int WBAESGenerator::generate8X8Bijections(CODING8X8_TABLE * tbl, size_t size, bool identity){
	unsigned int i=0,c=0;
	for(i=0; i<size; i++){
		c |= generate8X8Bijection(&tbl[i].coding, &tbl[i].invCoding, identity);
	}

	return c;
}

int WBAESGenerator::generate4X4Bijection(BIJECT4X4 *biject, BIJECT4X4 *invBiject, bool identity){
	if (!identity){
		return generateRandomBijection((unsigned char*)biject, (unsigned char*)invBiject, 16, 1);
	} else {
		int i;
		for(i=0; i<16; i++){
			(*biject)[i]    = i;
			(*invBiject)[i] = i;
		}; return 0;
	}
}

int WBAESGenerator::generate8X8Bijection(BIJECT8X8 *biject, BIJECT8X8 *invBiject, bool identity){
	if (!identity){
		return generateRandomBijection((unsigned char*)biject, (unsigned char*)invBiject, 256, 1);
	} else {
		int i;
		for(i=0; i<256; i++){
			(*biject)[i]    = i;
			(*invBiject)[i] = i;
		}; return 0;
	}
}

int WBAESGenerator::testWithVectors(bool coutOutput, WBAES * genAES, int extCodingFlags){
	// generate table implementation for given key
	ExtEncoding extc;

	//
	// Demonstrate also use of external encodings in practice
	generateExtEncoding(&extc, extCodingFlags);
	if (coutOutput){
		cout << "Generating table implementation for testvector key: " << endl;
		dumpVectorT(GenericAES::testVect128_key, 16);
	}

	//this->useIO04x04Identity=true;
	//this->useIO08x08Identity=true;
	//this->useDualAESARelationsIdentity=true;
	//this->useDualAESIdentity=true;
	//this->useMB08x08Identity=true;
	//this->useMB32x32Identity=true;

	generateTables(GenericAES::testVect128_key, KEY_SIZE_16, genAES, &extc, true);
	generateTables(GenericAES::testVect128_key, KEY_SIZE_16, genAES, &extc, false);

	//genAES.dumpEachRound=true;
	return this->testComputedVectors(coutOutput, genAES, &extc);
}

void WBAESGenerator::applyExternalEnc(W128b& state, ExtEncoding * extc, bool input){
	assert(extc!=NULL);
	if (input){
		// If input -> at first apply linear transformation 128 x 128, then bijection
		// Now we use output encoding G and quit, no MixColumn or Mixing bijections here.

		//
		// Mixing bijection 128x128
		//
		mat_GF2 tmpMat2(INIT_SIZE, 128, 1);
		for(int jj=0; jj<16; jj++){
			BYTE_to_matGF2(state.B[jj], tmpMat2, jj*8, 0);
		}
		tmpMat2 = extc->IODM[0].mb * tmpMat2;

		for(int jj=0; jj<16; jj++){
			state.B[jj] = matGF2_to_BYTE(tmpMat2, jj*8, 0);
		}

		//
		// IO bijection
		//
		for(int jj=0; jj<16; jj++){
			int tt = idxTranspose(jj);
			state.B[jj] = HILO(extc->lfC[0][2*tt+0].coding[HI(state.B[jj])], extc->lfC[0][2*tt+1].coding[LO(state.B[jj])]);
		}
	} else {
		// Output -> decode bijections

		//
		// IO bijection
		//
		for(int jj=0; jj<16; jj++){
			int tt = idxTranspose(jj);
			state.B[jj] = HILO(extc->lfC[1][2*tt+0].invCoding[HI(state.B[jj])], extc->lfC[1][2*tt+1].invCoding[LO(state.B[jj])]);
		}

		//
		// Mixing bijection 128x128
		//
		mat_GF2 tmpMat2(INIT_SIZE, 128, 1);
		for(int jj=0; jj<16; jj++){
			BYTE_to_matGF2(state.B[jj], tmpMat2, idxTranspose(jj)*8, 0);
		}
		tmpMat2 = extc->IODM[1].inv * tmpMat2;

		for(int jj=0; jj<16; jj++){
			state.B[jj] = matGF2_to_BYTE(tmpMat2, idxTranspose(jj)*8, 0);
		}
	}
}

void WBAESGenerator::applyExternalEnc(BYTE * state, ExtEncoding * extc, bool input, size_t numBlocks){
    if (state == nullptr || extc == nullptr){
        return;
    }

    W128b aesState{};

	for(int idx = 0; idx < numBlocks; ++idx) {
		arr_to_W128b(state, static_cast<size_t>(idx * N_BYTES), aesState);
		applyExternalEnc(aesState, extc, input);
		W128b_to_arr((char *) state, static_cast<size_t>(idx * N_BYTES), aesState);
	}
}

int WBAESGenerator::testComputedVectors(bool coutOutput, WBAES * genAES, ExtEncoding * extc){
	int i, err=0;

	// see [http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf]
	if (coutOutput){
		cout << "Testing Dual Whitebox AES generator implementation on test vectors..." << endl;
	}

	for(i=0; i<AES_TESTVECTORS; i++){
		W128b plain{}, cipher{}, state{};
		arr_to_W128b(GenericAES::testVect128_plain[i], 0, plain);
		arr_to_W128b(GenericAES::testVect128_plain[i], 0, state);
		arr_to_W128b(GenericAES::testVect128_cipher[i], 0, cipher);

		// encryption
		applyExternalEnc(state, extc, true);
		genAES->encrypt(state);
		applyExternalEnc(state, extc, false);
		if (coutOutput){
			cout << "Testvector index: " << i << endl;
			cout << "=====================" << endl;
			cout << "Testvector plaintext: " << endl;
			dumpW128b(plain); cout << endl;

			cout << "Testvector ciphertext: " << endl;
			dumpW128b(cipher); cout << endl;

			cout << "Enc(plaintext_test): " << endl;
			dumpW128b(state); cout << endl;
		}

		if (compare_W128b(state, cipher)){
			if (coutOutput) cout << "[  OK  ]    Enc(plaintext) == ciphertext_test" << endl;
		} else {
			err++;
			if (coutOutput) cout << "[ ERROR ]:  Enc(plaintext) != ciphertext_test" << endl;
#			ifdef FAIL
            FAIL() << "[ ERROR ]:  Enc(plaintext) != ciphertext_test";
#			endif
		}

		applyExternalEnc(state, extc, true);
		genAES->decrypt(state);
		applyExternalEnc(state, extc, false);
		if (coutOutput){
			cout << "Dec(Enc(plaintext_test)): " << endl;
			dumpW128b(state); cout << endl;
		}

		if (compare_W128b(state, plain)){
			if (coutOutput) cout << "[  OK  ]    Dec(Enc(plaintext)) == plaintext_test" << endl;
		} else {
			err++;
			if (coutOutput) cout << "[ ERROR ]:  Dec(Enc(plaintext)) != plaintext_test" << endl;
#			ifdef FAIL
			FAIL() << "[ ERROR ]:  Dec(Enc(plaintext)) != plaintext_test";
#			endif
		}

		if (coutOutput){
			cout << "==========================================================================================" << endl;
		}
	}

	return err;
}

int WBAESGenerator::save(const char * filename, WBAES * aes, ExtEncoding * extCoding){
#ifdef WBAES_BOOST_SERIALIZATION
	std::ofstream ofs(filename);
	int code = save(ofs, aes, extCoding);
	ofs.close();

	return code;
#else
	cerr << "WBAESGenerator::save: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

int WBAESGenerator::load(const char * filename, WBAES * aes, ExtEncoding * extCoding){
#ifdef WBAES_BOOST_SERIALIZATION
	std::ifstream ifs(filename);
	int code = load(ifs, aes, extCoding);
	ifs.close();
	return code;
#else
	cerr << "WBAESGenerator::load: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

int WBAESGenerator::save(ostream& out, WBAES * aes, ExtEncoding * extCoding){
#ifdef WBAES_BOOST_SERIALIZATION
	boost::archive::binary_oarchive oa(out);
	if (aes) {
		aes->save(oa);
	}
	if (extCoding) {
		oa << *extCoding;
	}

	return 0;
#else
	cerr << "WBAESGenerator::save: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

int WBAESGenerator::load(istream& ins, WBAES * aes, ExtEncoding * extCoding){
#ifdef WBAES_BOOST_SERIALIZATION
	boost::archive::binary_iarchive ia(ins);
	if (aes) {
		aes->load(ia);
	}
	if (extCoding){
		ia >> *extCoding;
	}

	return 0;
#else
	cerr << "WBAESGenerator::load: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}
