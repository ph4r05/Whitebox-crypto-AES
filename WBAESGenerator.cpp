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
	useDualAESIdentity=false;
	useIO04x04Identity=false;
	useIO08x08Identity=true;
	useMB08x08Identity=false;
	useMB32x32Identity=false;
}

WBAESGenerator::~WBAESGenerator() {
	;
}

void WBAESGenerator::generateCodingMap(WBACR_AES_CODING_MAP& map, int *codingCount, bool encrypt){
	int r,i,cIdx=0;

	// Encryption/Decryption dependent operation and tables
	int (&shiftOp)[N_BYTES] 							 = encrypt ? (this->shiftRowsLBijection) : (this->shiftRowsLBijectionInv);
	W08x32Coding  (&edT2)[N_ROUNDS][N_SECTIONS][4]		 = encrypt ? (map.eT2)  		  : (map.dT2);
	W08x32Coding  (&edT3)[N_ROUNDS][N_SECTIONS][4]       = encrypt ? (map.eT3)   		  : (map.dT3);
	CODING        (&edXOR1)[N_ROUNDS][N_SECTIONS][24]    = encrypt ? (map.eXOR1) 		  : (map.dXOR1);
	CODING        (&edXOR2)[N_ROUNDS][N_SECTIONS][24]    = encrypt ? (map.eXOR2) 		  : (map.dXOR2);

	//
	// In the last round there is no 04x04 bijective output mapping.
	// There are no XOR tables and T3 tables. Output from T2 tables is final and
	// encoded by output 128bit encoding G (not allocated & connected here).
	// Thus encode only round 1..9. Last round 9 output coding from XOR2 master table
	// is connected to T2 input coding in round 10 to revert last bijection.
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
			CONNECT_W08x32_TO_XOR_H(edT2[r][i][0], edXOR1[r][i], 0);	// HIGH part of XOR1_0 table IN connect to OUT of T2_0 table
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
			// Connect result XOR layer 4 to T2 boxes in next round
			newIdx = shiftOp[4*i+0];
			CONNECT_XOR_TO_W08x32(edXOR2[r][i], 16, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
			newIdx = shiftOp[4*i+1];
			CONNECT_XOR_TO_W08x32(edXOR2[r][i], 18, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
			newIdx = shiftOp[4*i+2];
			CONNECT_XOR_TO_W08x32(edXOR2[r][i], 20, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
			newIdx = shiftOp[4*i+3];
			CONNECT_XOR_TO_W08x32(edXOR2[r][i], 22, edT2[r+1][ newIdx / 4 ][ newIdx % 4]);
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

void WBAESGenerator::generateIO128Coding(CODING8X8_TABLE (&coding)[N_BYTES], bool identity){
	this->generate8X8Bijections((CODING8X8_TABLE*) &coding, N_BYTES, identity);
}

void WBAESGenerator::generateTables(BYTE *key, enum keySize ksize, WBAES& genAES, CODING8X8_TABLE* pCoding08x08, bool encrypt){
	WBACR_AES_CODING_MAP 		codingMap;
	int 						codingCount;
	int							i,j,r,b,k;
	GenericAES					defaultAES;

	// Initialize IO coding map (networked fashion of mappings)
	generateCodingMap(codingMap, &codingCount, true);

	// Preparing all 4Bits internal encoding/decoding bijections
	CODING4X4_TABLE* pCoding04x04 = new CODING4X4_TABLE[codingCount+1];
	this->generate4X4Bijections(pCoding04x04, codingCount+1, useIO04x04Identity);

	// Generate mixing bijections
	MB08x08_TABLE eMB_L08x08 [MB_CNT_08x08_ROUNDS][MB_CNT_08x08_PER_ROUND];
	MB32x32_TABLE eMB_MB32x32[MB_CNT_32x32_ROUNDS][MB_CNT_32x32_PER_ROUND];
	this->generateMixingBijections(eMB_L08x08, MB_CNT_08x08_ROUNDS, eMB_MB32x32, MB_CNT_32x32_ROUNDS);

	// Encryption/Decryption dependent functions and tables
	int (&nextTbox)[N_BYTES]     = encrypt ? (this->shiftRowsLBijection) : (this->shiftRowsLBijectionInv);
	int (&shiftRowsOp)[N_BYTES]  = encrypt ? (this->shiftRows)    		 : (this->shiftRowsInv);
	W32XTB        (&genAES_edXTab)[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS]  = encrypt ? genAES.eXTab    : genAES.dXTab;
	AES_TB_TYPE2  (&genAES_edTab2)[N_ROUNDS][N_BYTES]			 		= encrypt ? genAES.eTab2    : genAES.dTab2;
	AES_TB_TYPE3  (&genAES_edTab3)[N_ROUNDS][N_BYTES]			 		= encrypt ? genAES.eTab3    : genAES.dTab3;
	W08x32Coding  (&codingMap_edT2)[N_ROUNDS][N_SECTIONS][4]	        = encrypt ? codingMap.eT2   : codingMap.dT2;
	W08x32Coding  (&codingMap_edT3)[N_ROUNDS][N_SECTIONS][4]            = encrypt ? codingMap.eT3   : codingMap.dT3;
	CODING        (&codingMap_edXOR1)[N_ROUNDS][N_SECTIONS][24]         = encrypt ? codingMap.eXOR1 : codingMap.dXOR1;
	CODING        (&codingMap_edXOR2)[N_ROUNDS][N_SECTIONS][24]         = encrypt ? codingMap.eXOR2 : codingMap.dXOR2;

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
		int rndPolynomial = useDualAESIdentity ? 0 : rand() % AES_IRRED_POLYNOMIALS;
		int rndGenerator  = useDualAESIdentity ? 0 : rand() % AES_GENERATORS;
		this->AESCipher[i].initFromIndex(rndPolynomial, rndGenerator);

		// convert BYTE[] to key
		BYTEArr_to_vec_GF2E(key, ksize, vecKey[i/N_SECTIONS][i%N_SECTIONS]);
		this->AESCipher[i].applyT(vecKey[i/N_SECTIONS][i%N_SECTIONS]);

		// Prepare key schedule from vector representation of encryption key
		this->AESCipher[i].expandKey(vecRoundKey[i/N_SECTIONS][i%N_SECTIONS], vecKey[i/N_SECTIONS][i%N_SECTIONS], ksize);
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
			// T table construction (Type2)
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

					// In first round we apply 128bit bijection composed of 8bit bijections (F_i^{-1})
					// Has to take into account ShiftRows() effect.
					if (r==0){
						bb = pCoding08x08[ shiftRowsOp[N_SECTIONS*i + j] ].invCoding[bb];
					} else {
						// Decode input with IO coding
						bb = iocoding_encode08x08(bb, codingMap_edT2[r][i][j].IC, true, pCoding04x04, pCoding08x08);
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
						tmpGF2E += vecRoundKey[r][i][16*r + idxTranspose(shiftRowsOp[ j*4 + i ])];
					} else if(r==0) {
						// Decryption & first round => add k_10 to state.
						// Same logic applies here
						// AddRoundKey(State, k_10)  | -> InvShiftRows(State)
						// InvShiftRows(State)       | -> AddRoundKey(State, InvShiftRows(k_10))
						tmpGF2E += vecRoundKey[r][i][16*N_ROUNDS + idxTranspose(shiftRowsOp[ j*4 + i ])];
					}

					// SBox transformation with dedicated AES for this round and section
					// Encryption: ByteSub
					// Decryption: ByteSubInv
					GF2E tmpE = encrypt ?
							  this->AESCipher[r*4 + i].ByteSub(tmpGF2E)
							: this->AESCipher[r*4 + i].ByteSubInv(tmpGF2E);

					// Decryption case:
					// T(x) = Sbox(x) + k
					if (!encrypt){
						tmpE += vecRoundKey[r][i][16*(N_ROUNDS-r-1) + idxTranspose(j*4 + i)];
					}

					// If we are in last round we also have to add k_10, not affected by ShiftRows()
					if (r==N_ROUNDS-1){
						// Adding last encryption key (k_10) by special way is performed only in encryption
						if (encrypt) {
							tmpE += vecRoundKey[r][i][16*(r+1) + idxTranspose(j*4 + i)];
						}

						// revert last dual AES transformation here
						this->AESCipher[4*r + i].applyTinv(tmpE);

						// Now we use output encoding G and quit, no MixColumn or Mixing bijections here.
						bb = getLong(tmpE);
						bb = pCoding08x08[ N_SECTIONS*i + j ].coding[bb];
						genAES_edTab2[r][i*4+j][b].B[0] = bb;
						genAES_edTab2[r][i*4+j][b].B[1] = bb;
						genAES_edTab2[r][i*4+j][b].B[2] = bb;
						genAES_edTab2[r][i*4+j][b].B[3] = bb;
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
						mcres = r<(N_ROUNDS-1) ? this->AESCipher[r*4 + i].mixColMat * zj : zj;
					} else {
						mcres = r<(N_ROUNDS-1) ? this->AESCipher[r*4 + i].mixColInvMat * zj : zj;
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
					// Build tables - for each byte
					for(b=0; b<256; b++){
						int	 bb = b;
						CODING & xorCoding = j > 2 ? codingMap_edXOR2[r][i][(j-3)*8+k] : codingMap_edXOR1[r][i][j*8+k];
						bb = iocoding_encode08x08(bb, xorCoding.IC, true, pCoding04x04, pCoding08x08);
//cout << "XTB["<<r<<"]["<<i<<"]["<<j<<"]["<<k<<"]["<<b<<"] = " << HI(bb) << " ^ " << LO(bb) << " = ";
						bb = HI(bb) ^ LO(bb);
						bb = iocoding_encode08x08(bb, xorCoding.OC, false, pCoding04x04, pCoding08x08);
//cout << bb << endl;
						//
						//             ________________________ ROUND
						//            |   _____________________ Section with same AES structure/MixCol stripe
						//            |  |   __________________ Master XOR table in section (6 in total, 3 up, 3 down)
						//            |  |  |   _______________ Slave XOR table in master table, 8 in total
						//            |  |  |  |   ____________ Mapping for input value
						//            |  |  |  |  |
						genAES_edXTab[r][i][j][k][b] = bb;

					}
				}
			}
		}
	}


	delete[] pCoding04x04;
}

int WBAESGenerator::generate4X4Bijections(CODING4X4_TABLE * tbl, size_t size, bool identity){
	unsigned int i,c;
	for(i=0; i<size; i++){
		c |= generate4X4Bijection(&tbl[i].coding, &tbl[i].invCoding, identity);
	}

	return c;
}

int WBAESGenerator::generate8X8Bijections(CODING8X8_TABLE * tbl, size_t size, bool identity){
	unsigned int i,c;
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


int WBAESGenerator::testWithVectors(bool coutOutput, WBAES &genAES){

	// generate table implementation for given key
	int i, err=0;

	CODING8X8_TABLE coding[16];
	W128b plain, cipher, state;
	generateIO128Coding(coding, true);

	if (coutOutput){
		cout << "Generating table implementation for testvector key: " << endl;
		dumpVectorT(GenericAES::testVect128_key, 16);
	}

	generateTables(GenericAES::testVect128_key, KEY_SIZE_16, genAES, coding, true);
	generateTables(GenericAES::testVect128_key, KEY_SIZE_16, genAES, coding, false);

	// see [http://csrc.nist.gov/publications/fips/fips197/fips-197.pdf]
	if (coutOutput){
		cout << "Testing Dual Whitebox AES generator implementation on test vectors..." << endl;
	}

	for(i=0; i<AES_TESTVECTORS; i++){
		arr_to_W128b(GenericAES::testVect128_plain[i], 0, plain);
		arr_to_W128b(GenericAES::testVect128_plain[i], 0, state);
		arr_to_W128b(GenericAES::testVect128_cipher[i], 0, cipher);

		// encryption
		genAES.encrypt(state);
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
		}

		genAES.decrypt(state);
		if (coutOutput){
			cout << "Dec(Enc(plaintext_test)): " << endl;
			dumpW128b(state); cout << endl;
		}

		if (compare_W128b(state, plain)){
			if (coutOutput) cout << "[  OK  ]    Dec(Enc(plaintext)) == plaintext_test" << endl;
		} else {
			err++;
			if (coutOutput) cout << "[ ERROR ]:  Dec(Enc(plaintext)) != plaintext_test" << endl;
		}

		if (coutOutput){
			cout << "==========================================================================================" << endl;
		}
	}

	return err;
}
