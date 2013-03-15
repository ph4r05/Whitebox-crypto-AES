/*
 * WBAES.cpp
 *
 *  Created on: Mar 10, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <ctime>


// NTL dependencies
#include "WBAES.h"
#include "NTLUtils.h"
#include <iomanip>
NTL_CLIENT


// Shift rows selector
int WBAES::shiftRows[N_BYTES] = {
		 0,  5, 10, 15,
		 4,  9, 14,  3,
		 8, 13,  2,  7,
		12,  1,  6, 11
};

int WBAES::shiftRowsInv[N_BYTES] = {
		0, 13, 10, 7,
		4,  1, 14, 11,
		8,  5,  2, 15,
		12, 9,  6,  3
};

WBAES::WBAES() {
	;
}

WBAES::~WBAES() {
	;
}

void WBAES::encrypt(W128b& state){
	encdec(state, true);
}

void WBAES::decrypt(W128b& state){
	encdec(state, false);
}

void WBAES::encdec(W128b& state, bool encrypt){
	int r=0, i=0;
	W32b ires[N_BYTES];				// intermediate result for T-boxes

	// encryption/decryption dependent operations and tables
	int (&shiftOp)[N_BYTES] 							 = encrypt ? (this->shiftRows) : (this->shiftRowsInv);
	W32XTB (&edXTab)[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS] = encrypt ? (this->eXTab) 	   : (this->dXTab);
	AES_TB_TYPE2 (&edTab2)[N_ROUNDS][N_BYTES]			 = encrypt ? (this->eTab2) 	   : (this->dTab2);
	AES_TB_TYPE3 (&edTab3)[N_ROUNDS][N_BYTES]			 = encrypt ? (this->eTab3) 	   : (this->dTab3);


	for(r=0; r<N_ROUNDS; r++){
		// Perform rest of the operations on 4 tuples.
		for(i=0; i<N_BYTES; i+=4){
			// Apply type 2 tables to all bytes, counting also shift rows selector.
			ires[i+0].l = edTab2[r][i+0][state.B[shiftOp[i+0]]].l;
			ires[i+1].l = edTab2[r][i+1][state.B[shiftOp[i+1]]].l;
			ires[i+2].l = edTab2[r][i+2][state.B[shiftOp[i+2]]].l;
			ires[i+3].l = edTab2[r][i+3][state.B[shiftOp[i+3]]].l;

			// In the last round, result is directly in T2 boxes
			if (r==(N_ROUNDS-1)){
				state.B[i+0] = ires[i+0].B[0];
				state.B[i+1] = ires[i+1].B[0];
				state.B[i+2] = ires[i+2].B[0];
				state.B[i+3] = ires[i+3].B[0];
				continue;
			}

			// XOR results of T2 boxes
			op8xor(ires[i+0], ires[i+1], edXTab[r][i/4][0], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], edXTab[r][i/4][1], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], edXTab[r][i/4][2], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage

			// Apply T3 boxes, valid XOR results are in ires[0], ires[4], ires[8], ires[12]
			// Start from the end, because in ires[i] is our XORing result.
			//
			//                    ________________________ ROUND
			//                   |    ____________________ T3 box for 1 section
			//                   |   |      ______________ (1 xor 2) xor (3 xor 4)
			//                   |   |     |         _____ 8bit parts of 32 bit result
			//                   |   |     |        |
			ires[i+3].l = edTab3[r][i+3][ires[i].B[i+3]].l;
			ires[i+2].l = edTab3[r][i+2][ires[i].B[i+2]].l;
			ires[i+1].l = edTab3[r][i+1][ires[i].B[i+1]].l;
			ires[i+0].l = edTab3[r][i+0][ires[i].B[i+0]].l;


			// Apply XORs again, now on T3 results
			// Copy results back to state
			op8xor(ires[i+0], ires[i+1], edXTab[r][i/4][3], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], edXTab[r][i/4][4], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], edXTab[r][i/4][5], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage

			// Copy results back to state
			// ires[i] now contains 32bit XOR result
			state.l[i/4] = ires[i].l;
		}
	}
}
