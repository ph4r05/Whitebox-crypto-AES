/*
 * WBAES.cpp
 *
 *  Created on: Mar 10, 2013
 *      Author: ph4r05
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

int WBAES::shiftRowsLBijection[N_BYTES] = {
		0, 13, 10, 7,
		4,  1, 14, 11,
		8,  5,  2, 15,
		12, 9,  6,  3
};

WBAES::WBAES() {
	// TODO Auto-generated constructor stub

}

WBAES::~WBAES() {
	// TODO Auto-generated destructor stub
}

void WBAES::encrypt(W128b& state){
	int r=0, i=0;
	W32b ires[N_BYTES];				// intermediate result for T-boxes
	for(r=0; r<N_ROUNDS; r++){
		// Perform rest of the operations on 4 tuples.
		for(i=0; i<N_BYTES; i+=4){
			// Apply type 2 tables to all bytes, counting also shift rows selector.
			ires[i+0] = this->eTab2[r][i+0][state.B[shiftRows[i+0]]];
			ires[i+1] = this->eTab2[r][i+1][state.B[shiftRows[i+1]]];
			ires[i+2] = this->eTab2[r][i+2][state.B[shiftRows[i+2]]];
			ires[i+3] = this->eTab2[r][i+3][state.B[shiftRows[i+3]]];

			// XOR results of T2 boxes
			op8xor(ires[i+0], ires[i+1], this->eXTab[r][i/4][0], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], this->eXTab[r][i/4][1], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], this->eXTab[r][i/4][2], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage

			// Apply T3 boxes, valid XOR results are in ires[0], ires[4], ires[8], ires[12]
			// Start from the end, because in ires[i] is our XORing result.
			//
			//                       ________________________ ROUND
			//                      |    ____________________ T3 box for 1 section
			//                      |   |      ______________ (1 xor 2) xor (3 xor 4)
			//                      |   |     |         _____ 8bit parts of 32 bit result
			//                      |   |     |        |
			ires[i+3] = this->eTab3[r][i+3][ires[i].B[i+3]];
			ires[i+2] = this->eTab3[r][i+2][ires[i].B[i+2]];
			ires[i+1] = this->eTab3[r][i+1][ires[i].B[i+1]];
			ires[i+0] = this->eTab3[r][i+0][ires[i].B[i+0]];


			// Apply XORs again, now on T3 results
			// Copy results back to state
			op8xor(ires[i+0], ires[i+1], this->eXTab[r][i/4][3], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], this->eXTab[r][i/4][4], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], this->eXTab[r][i/4][5], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage

			// Copy results back to state
			// ires[i] now contains 32bit XOR result
			state.B[i+0] = ires[i].B[0];
			state.B[i+1] = ires[i].B[1];
			state.B[i+2] = ires[i].B[2];
			state.B[i+3] = ires[i].B[3];
		}
	}
}
