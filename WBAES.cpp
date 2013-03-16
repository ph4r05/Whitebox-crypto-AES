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
		 0,  1,  2,  3,
		 5,  6,  7,  4,
		10, 11,  8,  9,
		15, 12, 13, 14
};

int WBAES::shiftRowsInv[N_BYTES] = {
		 0,  1,  2,  3,
		 7,  4,  5,  6,
		10, 11,  8,  9,
		13, 14, 15, 12
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
//		cout << "inputState[" << r << "] dump: " << endl;
		dumpW128b(state);

		// Perform rest of the operations on 4 tuples.
		for(i=0; i<N_BYTES; i+=4){
			// Apply type 2 tables to all bytes, counting also shift rows selector.
			// One section ~ 1 column of state array, so select 1 column, first will
			// have indexes 0,4,8,12. Also take ShiftRows() into consideration.
			ires[i+0].l = edTab2[r][i+0][state.B[shiftOp[i/4+0*4]]].l;
			ires[i+1].l = edTab2[r][i+1][state.B[shiftOp[i/4+1*4]]].l;
			ires[i+2].l = edTab2[r][i+2][state.B[shiftOp[i/4+2*4]]].l;
			ires[i+3].l = edTab2[r][i+3][state.B[shiftOp[i/4+3*4]]].l;
//			cout << "Selecting: " << shiftOp[i/4+0*4] << ", "
//								  << shiftOp[i/4+1*4] << ", "
//								  << shiftOp[i/4+2*4] << ", "
//								  << shiftOp[i/4+3*4] << endl;
//
//			cout << "Selected bytes: " << CHEX(state.B[shiftOp[i/4+0*4]]) << ", "
//											  << CHEX(state.B[shiftOp[i/4+1*4]]) << ", "
//											  << CHEX(state.B[shiftOp[i/4+2*4]]) << ", "
//											  << CHEX(state.B[shiftOp[i/4+3*4]]) << endl;

			// In the last round, result is directly in T2 boxes
			if (r==(N_ROUNDS-1)){
				continue;
			}

//			cout << "T2[" << r << "][" << (i+0) << "] " << CHEX(ires[i].l) << endl;
//			cout << "T2[" << r << "][" << (i+1) << "] " << CHEX(ires[i+1].l) << endl;
//			cout << "T2[" << r << "][" << (i+2) << "] " << CHEX(ires[i+2].l) << endl;
//			cout << "T2[" << r << "][" << (i+3) << "] " << CHEX(ires[i+3].l) << endl;


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
			//                   |   |     |        ______ 8bit parts of 32 bit result
			//                   |   |     |       |
			ires[i+3].l = edTab3[r][i+3][ires[i].B[3]].l;
			ires[i+2].l = edTab3[r][i+2][ires[i].B[2]].l;
			ires[i+1].l = edTab3[r][i+1][ires[i].B[1]].l;
			ires[i+0].l = edTab3[r][i+0][ires[i].B[0]].l;

			// Apply XORs again, now on T3 results
			// Copy results back to state
			op8xor(ires[i+0], ires[i+1], edXTab[r][i/4][3], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], edXTab[r][i/4][4], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], edXTab[r][i/4][5], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage
		}

		//
		// Copy results back to state
		// ires[i] now contains 32bit XOR result
		// We have to copy result to column...
		for(i=0; i<N_BYTES; i+=4){
			state.B[i/4+ 0] = r<(N_ROUNDS-1) ? ires[i].B[0] : ires[i+0].B[0];
			state.B[i/4+ 4] = r<(N_ROUNDS-1) ? ires[i].B[1] : ires[i+1].B[0];
			state.B[i/4+ 8] = r<(N_ROUNDS-1) ? ires[i].B[2] : ires[i+2].B[0];
			state.B[i/4+12] = r<(N_ROUNDS-1) ? ires[i].B[3] : ires[i+3].B[0];
		}

		cout << "EndOfRound[" << r << "] dump: " << endl;
		dumpW128b(state);
	}
}
