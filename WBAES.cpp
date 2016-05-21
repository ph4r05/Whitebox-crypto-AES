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
const int WBAES::shiftRows[N_BYTES] = {
		 0,  1,  2,  3,
		 5,  6,  7,  4,
		10, 11,  8,  9,
		15, 12, 13, 14
};

const int WBAES::shiftRowsInv[N_BYTES] = {
		 0,  1,  2,  3,
		 7,  4,  5,  6,
		10, 11,  8,  9,
		13, 14, 15, 12
};

void arr_to_W128b(unsigned char * src, size_t offset, W128b& dst){
	int i=0;
	for(i=0; i<16; i++){
		dst.B[idxTranspose(i)] = src[offset+i];
	}
}

void arr_to_W128b(char * src, size_t offset, W128b& dst){
	int i=0;
	for(i=0; i<16; i++){
		dst.B[idxTranspose(i)] = src[offset+i];
	}
}

void W128b_to_arr(char * dst, size_t offset, W128b& src){
	int i=0;
	for(i=0; i<16; i++){
		dst[offset+i] = src.B[idxTranspose(i)];
	}
}

bool compare_W128b(const W128b& src, const W128b& dst){
	int i;
	for(i=0; i<16; i++){
		if (src.B[i]!=dst.B[i]) return false;
	}

	return true;
}

WBAES::WBAES() {
	dumpEachRound = false;
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
	W32b  ires[N_BYTES];				// intermediate result for T2,T3-boxes
	W128b ares[N_BYTES];				// intermediate result for T1-boxes

	// encryption/decryption dependent operations and tables
	const int (&shiftOp)[N_BYTES]                        = encrypt ? (this->shiftRows) : (this->shiftRowsInv);
	W32XTB (&edXTab)[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS] = encrypt ? (this->eXTab)     : (this->dXTab);
	W32XTB (&edXTabEx)[2][15][4]                         = encrypt ? (this->eXTabEx)   : (this->dXTabEx);
	AES_TB_TYPE1 (&edTab1)[2][N_BYTES]                   = encrypt ? (this->eTab1)     : (this->dTab1);
	AES_TB_TYPE2 (&edTab2)[N_ROUNDS][N_BYTES]            = encrypt ? (this->eTab2)     : (this->dTab2);
	AES_TB_TYPE3 (&edTab3)[N_ROUNDS][N_BYTES]            = encrypt ? (this->eTab3)     : (this->dTab3);
#ifdef AES_BGE_ATTACK
	GF256_func_t (&edOutputBijection)[N_ROUNDS][N_BYTES] = encrypt ? (this->eOutputBijection) : (this->dOutputBijection);
#endif

	// At first we have to put input to T1 boxes directly, no shift rows
	// compute result to ares[16]
	for(i=0; i<N_BYTES; i++){
		// Note: Tbox is indexed by cols, state by rows - transpose needed here
		W128CP(ares[i], edTab1[0][i][state.B[idxTranspose(i)]]);
	}

	// Now compute cascade of XOR tables
	// We have 8*32 XOR tables, they sum T1: [01] [23] [45] [67] [89] [1011] [1213] [1415]
	//                                        0     1   2     3   4     5      6      7
	// The task is to connect                  \   /     \   /     \   /        \    /
	//                                         [0123]    [4567]   [891011]    [12131415]
	//                                            8         9        10          11
	//                                             \       /          \          /
	//                                              \     /            \        /
	//                                             [01234567]       [89101112131415]
	//                                                 12                 13
	//                                                  \                 /
	//                                                   \               /
	//                                                [0123456789101112131415]
	//                                                           14
	//

	// 1st level of XORs
	for(i=0;i<N_BYTES;i+=2){
		op8xor_128(ares[i+0], ares[i+1], edXTabEx[0][i/2], ares[i+0]);  // 1 xor 2 --> 1
	}

	// Finish XOR cascade
	op8xor_128(ares[0],  ares[2],  edXTabEx[0][8],  ares[0]);  // 0  xor 2  --> 0
	op8xor_128(ares[4],  ares[6],  edXTabEx[0][9],  ares[4]);  // 4  xor 6  --> 4
	op8xor_128(ares[8],  ares[10], edXTabEx[0][10], ares[8]);  // 8  xor 10 --> 8
	op8xor_128(ares[12], ares[14], edXTabEx[0][11], ares[12]); // 12 xor 14 --> 12
	// 3. lvl
	op8xor_128(ares[0],  ares[4],  edXTabEx[0][12], ares[0]);  // 0 xor 4  --> 0
	op8xor_128(ares[8],  ares[12], edXTabEx[0][13], ares[8]);  // 8 xor 12 --> 8
	// 4. lvl - final stage. Result in ares[0]
	op8xor_128(ares[0],  ares[8],  edXTabEx[0][14], ares[0]);  // 0 xor 8 --> 0
	// Copy result from ares[0] to state
	W128CP(state, ares[0]);

	// Compute 9 rounds of T2 boxes
	for(r=0; r<(N_ROUNDS-1); r++){
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
			state.B[i/4+ 0] = ires[i].B[0];
			state.B[i/4+ 4] = ires[i].B[1];
			state.B[i/4+ 8] = ires[i].B[2];
			state.B[i/4+12] = ires[i].B[3];
		}

#ifdef AES_BGE_ATTACK
		// If we are performing attack, we modified output bijection for 1 byte from 2 concatenated 4x4 bijections to one 8x8
		for(i=0; i<N_BYTES; i+=4){
			state.B[i/4+ 0] = edOutputBijection[r][i/4+ 0][state.B[i/4+ 0]];
			state.B[i/4+ 4] = edOutputBijection[r][i/4+ 4][state.B[i/4+ 4]];
			state.B[i/4+ 8] = edOutputBijection[r][i/4+ 8][state.B[i/4+ 8]];
			state.B[i/4+12] = edOutputBijection[r][i/4+12][state.B[i/4+12]];
		}
#endif

		if (dumpEachRound){
			cout << "EndOfRound[" << r << "] dump: " << endl;
			dumpW128b(state);
		}
	}

	//
	// Final round is special -> T1 boxes
	//
	for(i=0; i<N_BYTES; i+=4){
		// Rules:
		//   1. i-th T1 table stores to ares[i]
		//   2. T1, T2 tables are indexed by column (0,1,2,3 = indexes for first column processing boxes)
		//   3. state is indexed by rows!
		W128CP(ares[i+0], edTab1[1][i+0][state.B[shiftOp[i/4+0*4]]]);
		W128CP(ares[i+1], edTab1[1][i+1][state.B[shiftOp[i/4+1*4]]]);
		W128CP(ares[i+2], edTab1[1][i+2][state.B[shiftOp[i/4+2*4]]]);
		W128CP(ares[i+3], edTab1[1][i+3][state.B[shiftOp[i/4+3*4]]]);
	}

	// and finally compute XOR cascade again, now for T1[1] - output T1
	// 1st level of XORs
	for(i=0;i<N_BYTES;i+=2){
		op8xor_128(ares[i+0], ares[i+1], edXTabEx[1][i/2], ares[i+0]);  // 1 xor 2 --> 1
	}

	// Finish XOR cascade by hand
	op8xor_128(ares[0],  ares[2],  edXTabEx[1][8],  ares[0]);  // 0  xor 2  --> 0
	op8xor_128(ares[4],  ares[6],  edXTabEx[1][9],  ares[4]);  // 4  xor 6  --> 4
	op8xor_128(ares[8],  ares[10], edXTabEx[1][10], ares[8]);  // 8  xor 10 --> 8
	op8xor_128(ares[12], ares[14], edXTabEx[1][11], ares[12]); // 12 xor 14 --> 12
	// 3. lvl
	op8xor_128(ares[0],  ares[4],  edXTabEx[1][12], ares[0]);  // 0 xor 4  --> 0
	op8xor_128(ares[8],  ares[12], edXTabEx[1][13], ares[8]);  // 8 xor 12 --> 8
	// 4. lvl - final stage. Result in ares[0]
	op8xor_128(ares[0],  ares[8],  edXTabEx[1][14], ares[0]);  // 0 xor 8 --> 0
	// Copy result from ares[0] to state, transpose, not W128CP(state, ares[0]);
	for(i=0; i<N_BYTES; i++){
		state.B[i] = ares[0].B[idxTranspose(i)];
	}

#ifdef AES_BGE_ATTACK
	// If we are performing attack, we modified output bijection for 1 byte from 2 concatenated 4x4 bijections to one 8x8
	for(i=0; i<N_BYTES; i+=4){
		state.B[i/4+ 0] = edOutputBijection[r][i/4+ 0][state.B[i/4+ 0]];
		state.B[i/4+ 4] = edOutputBijection[r][i/4+ 4][state.B[i/4+ 4]];
		state.B[i/4+ 8] = edOutputBijection[r][i/4+ 8][state.B[i/4+ 8]];
		state.B[i/4+12] = edOutputBijection[r][i/4+12][state.B[i/4+12]];
	}
#endif

	if (dumpEachRound){
		cout << "EndOfRound[" << r << "] dump: " << endl;
		dumpW128b(state);
	}
}

int WBAES::save(const char * filename){
#ifdef WBAES_BOOST_SERIALIZATION
	std::ofstream ofs(filename);
	int code = save(ofs);
	ofs.close();

	return code;
#else
	cerr << "WBAES::save: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

int WBAES::load(const char * filename){
#ifdef WBAES_BOOST_SERIALIZATION
	std::ifstream ifs(filename);
	int code = load(ifs);
	ifs.close();
	return code;
#else
	cerr << "WBAES::load: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

int WBAES::save(ostream& out){
#ifdef WBAES_BOOST_SERIALIZATION
	boost::archive::binary_oarchive oa(out);
	save(oa);

	return 0;
#else
	cerr << "WBAES::save: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

int WBAES::load(istream& ins){
#ifdef WBAES_BOOST_SERIALIZATION
	boost::archive::binary_iarchive ia(ins);
	load(ia);

	return 0;
#else
	cerr << "WBAES::load: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

std::string WBAES::save() {
#ifdef WBAES_BOOST_SERIALIZATION
	std::ostringstream out;
	save(out);
	return out.str();
#else
	cerr << "WBAES::save: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}

#ifdef WBAES_BOOST_SERIALIZATION
int WBAES::load(boost::archive::binary_iarchive& ins){
	ins >> *this;
	return 0;
}

int WBAES::save(boost::archive::binary_oarchive& out){
	out << *this;
	return 0;
}
#endif

int WBAES::loadString(std::string serialized) {
#ifdef WBAES_BOOST_SERIALIZATION
	std::istringstream inf(serialized);
	int code = load(inf);
	return code;
#else
	cerr << "WBAES::save: Boost is not enabled, use WBAES_BOOST_SERIALIZATION" << endl;
	return -1;
#endif
}
