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

// NTL dependencies
#include "WBAES.h"
#include "NTLUtils.h"

//NTL_CLIENT

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
}

