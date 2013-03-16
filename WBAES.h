/*
 * WBAES.h
 *
 *  Created on: Mar 10, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#ifndef WBAES_H_
#define WBAES_H_
#include "NTLUtils.h"
#include <iomanip>

// 4bit operations support on 8bit storage
#define HI(x)   (((x) >> 4) & 0xF)                      // HI(xxxxyyyy) = 0000xxxx 
#define LO(x)   ((x) & 0xF)                             // LO(xxxxyyyy) = 0000yyyy
#define HILO(h,l)   ((((h) & 0xF) << 4) | (l & 0xF))    // HILO(qqqqwwww, rrrrtttt) = wwwwtttt

// operations - for XOR tables, 2x4bit argument
#define OP2HI(h1,h2) HILO(HI(h1), HI(h2))               // OP2HI(qqqqwwww, rrrrtttt) = qqqqrrrr
#define OP2LO(h1,h2) HILO(LO(h1), LO(h2))               // OP2LO(qqqqwwww, rrrrtttt) = wwwwtttt

// extended HI/LO operations
#define HIIDX(x,i) (((x) >> 4*(i)) & 0xF)               // HI(qqqqwwwwrrrrtttt, 2) = wwww
#define HILO8(a,b,c,d,e,f,g,h) ( \
            (((a) & 0xF) << 4*7) | \
            (((b) & 0xF) << 4*6) | \
            (((c) & 0xF) << 4*5) | \
            (((d) & 0xF) << 4*4) | \
            (((e) & 0xF) << 4*3) | \
            (((f) & 0xF) << 4*2) | \
            (((g) & 0xF) << 4*1) | \
            (  h  & 0xF)     ) 

// constructs 32bit unsigned long from 4x bytes.
#define HILO4(a,b,c,d) ( \
            (((a) & 0xFF) << 8*3) | \
            (((b) & 0xFF) << 8*2) | \
            (((c) & 0xFF) << 8*1) | \
            (  d  & 0xFF)     )

//
// Some AES constants
//
#define N_ROUNDS        10      // AES rounds
#define N_SECTIONS      4       // 4 independent groups in one round
#define N_XOR_GROUPS    6       // 3 XOR tables to combine MIX Col, 3 XOR tables to combine MB
#define N_BYTES         16      // 16 bytes in one round

//
// EXTENDED DATA TYPES
//
typedef unsigned char BYTE;
typedef unsigned long DWORD;       
typedef BYTE          BITS4;                // FORM OF BITS4 is 0000xxxx, ONLY LOWER 4 BITS ARE USED
typedef BYTE    MCSTRIP[4];                 // partitional strip obtained by multiplication with MC            

typedef union _W32B{
    BYTE B[4];
    unsigned int l;
} W32b;

typedef union _W128B{
    BYTE B[16];
    unsigned int l[4];
} W128b;

// XOR table is 8b->4b mapping. Thus simple array of size 2^8=256 with type BYTE.
// 4bit type would be enough, but smallest possible is char, thus 8.
typedef BITS4 XTB[256];

// TYPE 1 tables
// DEF: G * INP
// Input is 1 byte, output is 128bit wide
typedef W128b AES_TB_TYPE1[256];

// TYPE 2 tables (T, Ty, MB boxes) 
// DEF: MB * Tyi * T * L2 ^{-1} (x)
// Input is 1 byte (2x BITS4), output is 32bit wide (after MC)
typedef W32b AES_TB_TYPE2[256];

// TYPE 3 tables
// DEF: L * MB ^{-1} (x)
// Input is 1 byte (2x BITS4), output is 32bit wide
typedef W32b AES_TB_TYPE3[256];
 
// 8 XOR tables for XORing 2x32bit input to obtain 32bit output 
typedef XTB    W32XTB[8];

//
// WBACR AES TABLES DEFINITIONS    
//
//typedef MCSTRIP ROUND_TABLE[256];           // 8 x 32 bits table
//typedef BITS4   XOR_TABLE[256];             // 8 x 4 bits table
//typedef BYTE    FINALROUND_TABLE[256];      // 8 x 8 bits table
//typedef BYTE    INVFIRSTROUND_TABLE[256];   // 8 x 8 bits table

// 32bit wide XOR, o1,o2 are of type W32b, xtb is of type W32XTB. 
// Returns unsigned long int value directly.
#define OP8XORlong(o1,o2,xtb) ((    HILO(xtb[0][OP2HI(o1[0], o2[0])] ,          \
                                         xtb[1][OP2LO(o1[0], o2[0])] ))         \
                                 | (HILO(xtb[2][OP2HI(o1[1], o2[1])] ,          \
                                         xtb[3][OP2LO(o1[1], o2[1])] ) <<  8)   \
                                 | (HILO(xtb[4][OP2HI(o1[2], o2[2])] ,          \
                                         xtb[5][OP2LO(o1[2], o2[2])] ) << 16)   \
                                 | (HILO(xtb[6][OP2HI(o1[3], o2[3])] ,          \
                                         xtb[7][OP2LO(o1[3], o2[3])] ) << 24))  \

// Copies W32b data from destination (d) to source (s)
#define W32CP(s,d) { s.l = d.l; }
	
// Copies W128b data from destination (d) to source (s)
#define W128CP(s,d) { s.l[0]  = d.l[0];   s.l[1] = d.l[1];  s.l[2] = d.l[2];  s.l[3] = d.l[3]; }

// 
// Simple 32bit wide XOR operation.
// O1, O2 are W32b operands. 
// Result = O1 XOR O2
//
#define OP8XOR(o1, o2, xtb, res) {                                             \
    res.B[0] = HILO(xtb[0][OP2HI(o1.B[0], o2.B[0])] ,                          \
                    xtb[1][OP2LO(o1.B[0], o2.B[0])] );                         \
    res.B[1] = HILO(xtb[2][OP2HI(o1.B[1], o2.B[1])] ,                          \
                    xtb[3][OP2LO(o1.B[1], o2.B[1])] );                         \
    res.B[2] = HILO(xtb[4][OP2HI(o1.B[2], o2.B[2])] ,                          \
                    xtb[5][OP2LO(o1.B[2], o2.B[2])] );                         \
    res.B[3] = HILO(xtb[6][OP2HI(o1.B[3], o2.B[3])] ,                          \
                    xtb[7][OP2LO(o1.B[3], o2.B[3])] ); }

/**
 * Simple 32bit wide XOR operation.
 * O1, O2 are W32b operands. 
 * Result = O1 XOR O2
 * 
 * Typesafe wrapper for macro OP8XOR.
 */
inline void op8xor(const W32b& o1, const W32b& o2, const W32XTB& xtb, W32b& res){
	W32b tmpRes, a1, a2;
	a1.l = o1.l;
	a2.l = o2.l;
	tmpRes.l = o1.l ^ o2.l;

    OP8XOR(o1, o2, xtb, res);

    if (res.l != tmpRes.l){
    	cout << "XOR warning!!! expected: " << CHEX(tmpRes.l) << " but got: " << CHEX(res.l) << " = " << CHEX(a1.l) << " ^ " << CHEX(a2.l) << endl;
    }
}

inline void dumpW128b(W128b& a){
	int i;
	for(i=0; i<16; i++){
		std::cout << CHEX((int)a.B[i]) << " ";
		if ((i+1)%4 == 0) std::cout << std::endl;
	}
}

inline void dumpW32b(W32b& a){
	int i;
	for(i=0; i<4; i++){
		std::cout << CHEX((int)a.B[i]) << " ";
	}
	std::cout << std::endl;
}

class WBAES {
public:
	WBAES();
	virtual ~WBAES();
	
	// How shift rows affects state array - indexes.
	// Selector to state array in the beggining of encryption round
	// 
	// Effect of shiftRows operation:
	//
	// | 00 04 08 12 |                          | 00 04 08 12 |
	// | 01 05 09 13 | ---   Shift Rows   --->  | 05 09 13 01 |
	// | 02 06 10 14 |  (cyclic left shift)     | 10 14 02 06 |
	// | 03 07 11 15 |                          | 15 03 07 11 |
	//
	static int shiftRows[N_BYTES];

    // Inverse ShiftRows()
    // | 00 04 08 12 |                          | 00 04 08 12 |
	// | 01 05 09 13 | --- Shift Rows Inv --->  | 13 01 05 09 |
	// | 02 06 10 14 |  (cyclic left right)     | 10 14 02 06 |
	// | 03 07 11 15 |                          | 07 11 15 03 |
	//
    static int shiftRowsInv[N_BYTES];

		
	// XOR tables
	W32XTB eXTab[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS];
	
	// Type I - just first round
	// @deprecated - useless here
	AES_TB_TYPE1 eFirstRoundTab;
	
	// Type II tables
    AES_TB_TYPE2 eTab2[N_ROUNDS][N_BYTES];
    
    // Type III tables
    AES_TB_TYPE3 eTab3[N_ROUNDS][N_BYTES];
    
    // universal encryption/decryption method
    void encdec(W128b& state, bool encrypt);

    // pure table implementation of encryption of given state
    void encrypt(W128b& state);

    // pure table implementation of decryption of given state
    void decrypt(W128b& state);

    //
    // Decryption tables
    //
    // XOR tables
    W32XTB dXTab[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS];

    // Type II tables
    AES_TB_TYPE2 dTab2[N_ROUNDS][N_BYTES];

    // Type III tables
    AES_TB_TYPE3 dTab3[N_ROUNDS][N_BYTES];
};


#endif /* WBAES_H_ */
