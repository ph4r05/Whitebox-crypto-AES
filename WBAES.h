/*
 * WBAES.h
 *
 *  Created on: Mar 10, 2013
 *      Author: ph4r05
 */

#ifndef WBAES_H_
#define WBAES_H_

// 4bit operations support on 8bit storage
#define HI(x)   (((x) >> 4) & 0xF)                      // HI(xxxxyyyy) = 0000xxxx 
#define LO(x)   ((x) & 0xF)                             // LO(xxxxyyyy) = 0000yyyy
#define HILO(h,l)   ((((h) & 0xF) << 4) | (l & 0xF))    // HILO(qqqqwwww, rrrrtttt) = wwwwtttt

// operations - for XOR tables, 2x4bit argument
#define OP2HI(h1,h2) HILO(HI(h1), HI(h2))               // OP2HI(qqqqwwww, rrrrtttt) = qqqqrrrr
#define OP2LO(h1,h2) HILO(LO(h1), LO(h2))               // OP2LO(qqqqwwww, rrrrtttt) = wwwwtttt

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
    unsigned long int l;
} W32b;

typedef union _W128B{
    BYTE B[16];
    unsigned long int l[4];
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
                    xtb[5][OP2LO(o1.B[3], o2.B[2])] );                         \
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
    OP8XOR(o1, o2, xtb, res);
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
	
	// Effect of shift rows operation for L bijections in T3 tables.
	// DEF:  shiftRowsLBijection[i] = to which T2 table in next round
	//       will be i-th byte of state passed as input from this round.
	//
    // With this information we can construct L^r OUT bijection in T3 tables
    // to match L^{r+1, -1} IN bijection in T2 tables in next round.
    //
    // Every round operates on state byte in this way 
    // Upper row - which byte is selected, lower row - which byte is stored
    //
    // 00 05 10 15 | 04 09 14 03 | 08 12 02 07 | 12 01 06 11   |
    // -----------------------------------------------------   v
    // 00 01 02 03 | 04 05 06 07 | 08 09 10 11 | 12 13 14 15
    // 
    // In next round, shift rows operation will be used again, so 
    // This table gives prescript how input bytes will be 
    // mapped to T_boxes in next round (upper row), counting with 
    // shift operation in the beggining of the next round.
    //
    // Example: First 4 bytes taken from byte array 0,5,10,15 (shift rows effect) 
    //      will feed T2_0,T2_1,T2_2,T2_3 boxes.
    // 
    // 00 01 02 03 | 04 05 06 07 | 08 09 10 11 | 12 13 14 15   |
    // -----------------------------------------------------   v
    // 00 13 10 07 | 04 01 14 11 | 08 05 02 15 | 12 09 06 03
	static int shiftRowsLBijection[N_BYTES];
	
	// XOR tables
	W32XTB eXTab[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS];
	
	// Type I - just first round
	AES_TB_TYPE1 eFirstRoundTab;
	
	// Type II tables
    AES_TB_TYPE2 eTab2[N_ROUNDS][N_BYTES];
    
    // Type III tables
    AES_TB_TYPE3 eTab3[N_ROUNDS][N_BYTES];
    
    // pure table implementation of encryption of given state
    void encrypt(W128b& state);
};


#endif /* WBAES_H_ */
