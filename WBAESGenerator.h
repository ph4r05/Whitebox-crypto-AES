/*
 * WBAESGenerator.h
 *
 *  Created on: Mar 10, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#ifndef WBAESGENERATOR_H_
#define WBAESGENERATOR_H_

#include "base.h"
#include <assert.h>
#include <string.h>

#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_long.h>
#include <math.h>
#include <vector>
#include "NTLUtils.h"
#include "GenericAES.h"
#include "WBAES.h"
#include "MixingBijections.h"

#ifdef WBAES_BOOST_SERIALIZATION
#include <cstddef> // NULL
#include <iostream>
#include <fstream>
#include <string>
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/split_free.hpp>
#endif

#define NO_CODING           0x00000000  // IDENTITY CODING
#define UNASSIGNED_CODING   0xFFFFFFFF  // INVALID CODING
#define UNUSED_CODING       0xFFFFFFFE  // This coding is not in use (XOR tables use only lower 4 bits for result)
#define USE_IDENTITY_CODING(idx) ((idx) == NO_CODING || (idx) == UNASSIGNED_CODING || (idx) == UNUSED_CODING)
// VALID CODINGS ORDINARY NUMBER IS FROM 0x00000000 TO 0xFFFFFFFE (TOTAL COUNT == 2^32 - 1) 

// CODING SIZE TYPE
#define COD_BITS_UNASSIGNED 0x00
#define COD_BITS_4          0x01
#define COD_BITS_8          0x02
#define COD_BITS_8_EXT      0x03

// MIXING BIJECTION TYPE
#define MB_IDENTITY 0x00
#define MB_8x8      0x01
#define MB_32x32    0x02
#define MB_128x128  0x04

// MIXING BIJECTION COUNTS
#define MB_CNT_08x08_ROUNDS 9
#define MB_CNT_08x08_PER_ROUND 16
#define MB_CNT_32x32_ROUNDS 9
#define MB_CNT_32x32_PER_ROUND 4

// NUMBER OF XOR TABLES FOR ONE T1 TABLE
#define XTB_CNT_T1 480

//
//  HIGHLOW, DEFINE TWO 4-BITS CODING FOR 8-BITS ARGUMENT
//
typedef struct _HIGHLOW {
    BYTE    type;   // CODING SIZE TYPE. CURRENTLY DEFINED COD_BITS_4 & COD_BITS_8   
    
    DWORD   H;      // HIGH 4-BITS CODING (H == L for COD_BITS_8)
    DWORD   L;      // LOW 4-BITS CODING

    _HIGHLOW() {
        type = COD_BITS_4;            // DEFAULT IS COD_BITS_4
        H = UNASSIGNED_CODING;
        L = UNASSIGNED_CODING;
    }
} HIGHLOW;

//
//  CODING, DEFINE INPUT AND OUTPUT WBACR AES CODING FOR 8-BITS ARGUMENT
//
typedef struct _CODING {
    HIGHLOW IC;  
    HIGHLOW OC;
} CODING;

// BIJECTIONS DEFINITIONS
typedef BITS4   BIJECT4X4[16];   // 4 x 4 bits table => 16 row 
typedef BYTE    BIJECT8X8[256];  // 8 x 8 bits table => 256 row
typedef int     BIJECT8X8_EX[256];  // 8 x 8 bits table => 256 row,

//
//  4-BITS TO 4-BITS BIJECTION 
//
typedef struct _CODING4X4_TABLE {
    BIJECT4X4   coding;
    BIJECT4X4   invCoding;          // SPEED OPTIMALIZATION, CAN BE ALSO COMPUTED FROM coding MEMBER (DUE TO BIJECTION PROPERTY)
    
    _CODING4X4_TABLE(void) {
        memset(coding, 0, sizeof(BIJECT4X4));
        memset(invCoding, 0, sizeof(BIJECT4X4));
    }   
} CODING4X4_TABLE;

//
//  8-BITS TO 8-BITS BIJECTION 
//
typedef struct _CODING8X8_TABLE {
    BIJECT8X8   coding;
    BIJECT8X8   invCoding;          // SPEED OPTIMALIZATION, CAN BE ALSO COMPUTED FROM coding MEMBER (DUE TO BIJECTION PROPERTY)
    
    _CODING8X8_TABLE(void) {
        memset(coding, 0, sizeof(BIJECT8X8));
        memset(invCoding, 0, sizeof(BIJECT8X8));
    }   
} CODING8X8_TABLE;

//
//  8-BITS TO 8-BITS BIJECTION, LARGER DATA TYPE USED (int) INSTEAD OF BYTE 
//  SOME SPECIAL VALUES CAN BE ASSIGNED (NOT FROM RANGE <0, 255>) 
//
typedef struct _CODING8X8_TABLE_EX {
    BIJECT8X8_EX   coding;
    BIJECT8X8_EX   invCoding;      // SPEED OPTIMALIZATION, CAN BE ALSO COMPUTED FROM coding MEMBER (DUE TO BIJECTION PROPERTY)
    
    _CODING8X8_TABLE_EX(void) {
        clear();
    }   

    void clear() {
        memset(coding, 0xFF, sizeof(BIJECT8X8_EX));
        memset(invCoding, 0xFF, sizeof(BIJECT8X8_EX));
    }
} CODING8X8_TABLE_EX;

//
//  Mixing bijection (linear transformation represented as GF(2) matrix)
//
typedef struct _MB_TABLE {
	//int         type;
    NTL::mat_GF2     mb;
    NTL::mat_GF2     inv;          // SPEED OPTIMALIZATION, CAN BE ALSO COMPUTED FROM coding MEMBER (DUE TO BIJECTION PROPERTY)
    
    _MB_TABLE(void) {

    }   
} MB_TABLE;

typedef MB_TABLE MB08x08_TABLE;
typedef MB_TABLE MB32x32_TABLE;
typedef MB_TABLE MB128x128_TABLE;

//
// Coding for T2 and T3 boxes, 8bit -> 32bit
//
typedef struct _W08x32Coding {
    HIGHLOW IC;  
    HIGHLOW OC[4];
} W08x32Coding;

//
// Coding for T1 boxes, 8bit -> 128bit
//
typedef struct _W08x128Coding {
    HIGHLOW IC;
    HIGHLOW OC[16];
} W08x128Coding;

//
// Input/Output encoding. It is specification for T1 tables for apps using WBAES.
//
typedef struct _ExtEncoding {
	CODING4X4_TABLE lfC[2][2*N_BYTES];		// 0=>first input round bijection, 1=>last output round bijection
	MB128x128_TABLE IODM[2];                // 128 x 128 GF(2) matrix, input, output mixing bijection
	int flags = 0;							// flags the structure was created with
} ExtEncoding;

#define WBAESGEN_EXTGEN_fCID 1          // lfC[0]  in ExtEncoding will be identity
#define WBAESGEN_EXTGEN_lCID 2          // lfC[1]  in ExtEncoding will be identity
#define WBAESGEN_EXTGEN_IDMID 4         // IODM[0] in ExtEncoding will be identity
#define WBAESGEN_EXTGEN_ODMID 8         // IODM[1] in ExtEncoding will be identity

// whole ExtEncoding will be identity
#define WBAESGEN_EXTGEN_ID (WBAESGEN_EXTGEN_fCID | WBAESGEN_EXTGEN_lCID | WBAESGEN_EXTGEN_IDMID | WBAESGEN_EXTGEN_ODMID)

// NOTE: 
// Coding for XOR boxes can be done with CODING type easily
//
//  WBACR_AES_CODING_MAP DEFINES ASSIGNED CODING TO EACH VALUE OF WBACR_AES_TABLE
//
typedef CODING  CODING32W[4];
typedef struct _WBACR_AES_CODING_MAP {
    // ENCRYPT TABLES
	W08x128Coding       eT1[2][N_BYTES];                              // ENCRYPT ROUND CODING MAP
    W08x32Coding        eT2[N_ROUNDS][N_SECTIONS][4];                 // ENCRYPT ROUND CODING MAP
    W08x32Coding        eT3[N_ROUNDS][N_SECTIONS][4];                 // ENCRYPT ROUND CODING MAP   
    CODING              eXOR1[N_ROUNDS][N_SECTIONS][24];              // 24 == 8 4-BITS PARTS OF 32-BITS DWORD, used with T2 tables
    CODING              eXOR2[N_ROUNDS][N_SECTIONS][24];              // 24 == 8 4-BITS PARTS OF 32-BITS DWORD, used with T3 tables
    CODING              eXOR3[2][XTB_CNT_T1];					      // 15*4*8, used with T1 tables
    
    // DECRYPT TABLES
    W08x128Coding       dT1[2][N_BYTES];                              // ENCRYPT ROUND CODING MAP
    W08x32Coding        dT2[N_ROUNDS][N_SECTIONS][4];                 // ENCRYPT ROUND CODING MAP
	W08x32Coding        dT3[N_ROUNDS][N_SECTIONS][4];                 // ENCRYPT ROUND CODING MAP
	CODING              dXOR1[N_ROUNDS][N_SECTIONS][24];              // 24 == 8 4-BITS PARTS OF 32-BITS DWORD, used with T2 tables
	CODING              dXOR2[N_ROUNDS][N_SECTIONS][24];              // 24 == 8 4-BITS PARTS OF 32-BITS DWORD, used with T3 tables
	CODING              dXOR3[2][XTB_CNT_T1];					      // 15*4*8, used with T1 tables
} WBACR_AES_CODING_MAP, *PWBACR_AES_CODING_MAP;

//
// Allocates new 4X4 encodings for 08x32 tables (T2,T3) from given offset (can be used to allocate also T1)
// Allocation = generate unique bijection ID for particular IO box.
// Only OC (output coding) is generated = donor of the bijection. IC = acceptor and is set by CONNECT* macros
// From other tables OC fields.
//
#define ALLOCW08x32CodingEx(cod, ofs, idx) {                  \
    cod.OC[(ofs)+0].type = COD_BITS_4; cod.OC[(ofs)+1].type = COD_BITS_4; \
    cod.OC[(ofs)+2].type = COD_BITS_4; cod.OC[(ofs)+3].type = COD_BITS_4; \
    assert(cod.OC[(ofs)+0].H==UNASSIGNED_CODING); cod.OC[(ofs)+0].H = ++(idx);        \
    assert(cod.OC[(ofs)+0].L==UNASSIGNED_CODING); cod.OC[(ofs)+0].L = ++(idx);        \
    assert(cod.OC[(ofs)+1].H==UNASSIGNED_CODING); cod.OC[(ofs)+1].H = ++(idx);        \
    assert(cod.OC[(ofs)+1].L==UNASSIGNED_CODING); cod.OC[(ofs)+1].L = ++(idx);        \
    assert(cod.OC[(ofs)+2].H==UNASSIGNED_CODING); cod.OC[(ofs)+2].H = ++(idx);        \
    assert(cod.OC[(ofs)+2].L==UNASSIGNED_CODING); cod.OC[(ofs)+2].L = ++(idx);        \
    assert(cod.OC[(ofs)+3].H==UNASSIGNED_CODING); cod.OC[(ofs)+3].H = ++(idx);        \
    assert(cod.OC[(ofs)+3].L==UNASSIGNED_CODING); cod.OC[(ofs)+3].L = ++(idx);        };

#define ALLOCW08x32Coding(cod, idx) ALLOCW08x32CodingEx(cod, 0, idx)

//
// Allocate T1 tables - generate bijection IDs for output side of the table (128-bit wide)
//
#define ALLOCW08x128Coding(cod, idx) {                            \
	ALLOCW08x32CodingEx(cod, 0, idx);                             \
	ALLOCW08x32CodingEx(cod, 4, idx);                             \
	ALLOCW08x32CodingEx(cod, 8, idx);                             \
	ALLOCW08x32CodingEx(cod, 12, idx);                            };

//
// Allocates new output coding for 4-bit XOR boxes XTB[offset+0 - offset+7], altogether 32 bit XOR table
// Recall that output of XOR is stored in LOW part, thus upper is unused -> no allocation for upper part.
//
#define ALLOCXORCoding(xtb, offset, idx) {                                                                           \
    xtb[(offset)+0].OC.type = COD_BITS_4; xtb[(offset)+0].OC.H = UNUSED_CODING; assert(xtb[(offset)+0].OC.L==UNASSIGNED_CODING); xtb[(offset)+0].OC.L = ++(idx);      \
    xtb[(offset)+1].OC.type = COD_BITS_4; xtb[(offset)+1].OC.H = UNUSED_CODING; assert(xtb[(offset)+1].OC.L==UNASSIGNED_CODING); xtb[(offset)+1].OC.L = ++(idx);      \
    xtb[(offset)+2].OC.type = COD_BITS_4; xtb[(offset)+2].OC.H = UNUSED_CODING; assert(xtb[(offset)+2].OC.L==UNASSIGNED_CODING); xtb[(offset)+2].OC.L = ++(idx);      \
    xtb[(offset)+3].OC.type = COD_BITS_4; xtb[(offset)+3].OC.H = UNUSED_CODING; assert(xtb[(offset)+3].OC.L==UNASSIGNED_CODING); xtb[(offset)+3].OC.L = ++(idx);      \
    xtb[(offset)+4].OC.type = COD_BITS_4; xtb[(offset)+4].OC.H = UNUSED_CODING; assert(xtb[(offset)+4].OC.L==UNASSIGNED_CODING); xtb[(offset)+4].OC.L = ++(idx);      \
    xtb[(offset)+5].OC.type = COD_BITS_4; xtb[(offset)+5].OC.H = UNUSED_CODING; assert(xtb[(offset)+5].OC.L==UNASSIGNED_CODING); xtb[(offset)+5].OC.L = ++(idx);      \
    xtb[(offset)+6].OC.type = COD_BITS_4; xtb[(offset)+6].OC.H = UNUSED_CODING; assert(xtb[(offset)+6].OC.L==UNASSIGNED_CODING); xtb[(offset)+6].OC.L = ++(idx);      \
    xtb[(offset)+7].OC.type = COD_BITS_4; xtb[(offset)+7].OC.H = UNUSED_CODING; assert(xtb[(offset)+7].OC.L==UNASSIGNED_CODING); xtb[(offset)+7].OC.L = ++(idx);      };

//
// Allocates XOR table 128-bit wide
//
#define ALLOCXOR128Coding(xtb, offset, idx) {                                                                        \
		ALLOCXORCoding(xtb, (offset)+0,  idx);                                                                       \
		ALLOCXORCoding(xtb, (offset)+8,  idx);                                                                       \
		ALLOCXORCoding(xtb, (offset)+16, idx);                                                                       \
		ALLOCXORCoding(xtb, (offset)+24, idx);                                                                       };

//
// Connects OUTPUT coding of 32bit wide boxes (T2,T3) to INPUT coding of XOR boxes, 32bit wide. 
// Each XOR box accepts 2 arguments, first in HIGH part, second in LOW part, thus when associating
// mapping from one particular W32box we are using either HIGH or LOW parts. 
//
#define CONNECT_W08x32_TO_XOR_EX(cod, xtb, HL, offsetL, offsetR) { \
	assert(xtb[(offsetL)+0].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+0].H!=UNASSIGNED_CODING); xtb[(offsetL)+0].IC.HL = cod.OC[(offsetR)+0].H;                \
	assert(xtb[(offsetL)+1].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+0].L!=UNASSIGNED_CODING); xtb[(offsetL)+1].IC.HL = cod.OC[(offsetR)+0].L;                \
	assert(xtb[(offsetL)+2].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+1].H!=UNASSIGNED_CODING); xtb[(offsetL)+2].IC.HL = cod.OC[(offsetR)+1].H;                \
	assert(xtb[(offsetL)+3].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+1].L!=UNASSIGNED_CODING); xtb[(offsetL)+3].IC.HL = cod.OC[(offsetR)+1].L;                \
	assert(xtb[(offsetL)+4].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+2].H!=UNASSIGNED_CODING); xtb[(offsetL)+4].IC.HL = cod.OC[(offsetR)+2].H;                \
	assert(xtb[(offsetL)+5].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+2].L!=UNASSIGNED_CODING); xtb[(offsetL)+5].IC.HL = cod.OC[(offsetR)+2].L;                \
	assert(xtb[(offsetL)+6].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+3].H!=UNASSIGNED_CODING); xtb[(offsetL)+6].IC.HL = cod.OC[(offsetR)+3].H;                \
	assert(xtb[(offsetL)+7].IC.HL==UNASSIGNED_CODING && cod.OC[(offsetR)+3].L!=UNASSIGNED_CODING); xtb[(offsetL)+7].IC.HL = cod.OC[(offsetR)+3].L;                }

#define CONNECT_W08x32_TO_XOR_H_EX(cod, xtb, offsetL, offsetR)   CONNECT_W08x32_TO_XOR_EX(cod, xtb, H, offsetL, offsetR)
#define CONNECT_W08x32_TO_XOR_L_EX(cod, xtb, offsetL, offsetR)   CONNECT_W08x32_TO_XOR_EX(cod, xtb, L, offsetL, offsetR)

#define CONNECT_W08x32_TO_XOR(cod, xtb, HL, offset) CONNECT_W08x32_TO_XOR_EX(cod, xtb, HL, offset, 0)
#define CONNECT_W08x32_TO_XOR_H(cod, xtb, offset)   CONNECT_W08x32_TO_XOR_H_EX(cod, xtb, offset, 0)
#define CONNECT_W08x32_TO_XOR_L(cod, xtb, offset)   CONNECT_W08x32_TO_XOR_L_EX(cod, xtb, offset, 0)

//
// Connects OUTPUT coding for XOR tables to INPUT coding of XOR tables on lower layer.
// Has effect of combining result of 2XOR tables to input of 1 XOR table.
//
// Recall that XOR result is always stored in lower part of XOR, thus on the left side we
// are using OC.L;
//
// 1 XOR table accepts input from 2 sources. 
// In HIGH part is first argument, in LOW part is the second. Same functionality as
// in CONNECT_W08x32_TO_XOR macro
//
// This macro accepts XOR tables 32bit wide.
#define CONNECT_XOR_TO_XOR(xtb1, offset1, xtb3, offset3, HL) {                  \
    assert(xtb3[(offset3)+0].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+0].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+0].IC.HL = xtb1[(offset1)+0].OC.L;                          \
    assert(xtb3[(offset3)+1].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+1].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+1].IC.HL = xtb1[(offset1)+1].OC.L;                           \
    assert(xtb3[(offset3)+2].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+2].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+2].IC.HL = xtb1[(offset1)+2].OC.L;                           \
    assert(xtb3[(offset3)+3].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+3].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+3].IC.HL = xtb1[(offset1)+3].OC.L;                           \
    assert(xtb3[(offset3)+4].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+4].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+4].IC.HL = xtb1[(offset1)+4].OC.L;                           \
    assert(xtb3[(offset3)+5].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+5].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+5].IC.HL = xtb1[(offset1)+5].OC.L;                           \
    assert(xtb3[(offset3)+6].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+6].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+6].IC.HL = xtb1[(offset1)+6].OC.L;                           \
    assert(xtb3[(offset3)+7].IC.HL==UNASSIGNED_CODING && xtb1[(offset1)+7].OC.L!=UNASSIGNED_CODING); xtb3[(offset3)+7].IC.HL = xtb1[(offset1)+7].OC.L;                           }

#define CONNECT_XOR_TO_XOR_128(xtb1, offset1, xtb3, offset3, HL) {             \
        CONNECT_XOR_TO_XOR(xtb1, (offset1)+0,  xtb3, (offset3)+0,  HL);        \
        CONNECT_XOR_TO_XOR(xtb1, (offset1)+8,  xtb3, (offset3)+8,  HL);        \
        CONNECT_XOR_TO_XOR(xtb1, (offset1)+16, xtb3, (offset3)+16, HL);        \
        CONNECT_XOR_TO_XOR(xtb1, (offset1)+24, xtb3, (offset3)+24, HL);        }

#define CONNECT_XOR_TO_XOR_H(xtb1, offset1, xtb3, offset3) CONNECT_XOR_TO_XOR(xtb1, offset1, xtb3, offset3, H)
#define CONNECT_XOR_TO_XOR_L(xtb1, offset1, xtb3, offset3) CONNECT_XOR_TO_XOR(xtb1, offset1, xtb3, offset3, L)
#define CONNECT_XOR_TO_XOR_128_H(xtb1, offset1, xtb3, offset3) CONNECT_XOR_TO_XOR_128(xtb1, offset1, xtb3, offset3, H)
#define CONNECT_XOR_TO_XOR_128_L(xtb1, offset1, xtb3, offset3) CONNECT_XOR_TO_XOR_128(xtb1, offset1, xtb3, offset3, L)

//
// Connects 8bit output from 2 consecutive XOR tables to 8b input of W08x32 table
//
#define CONNECT_XOR_TO_W08x32(xtb, offset, cod) {                             \
    cod.IC.type = xtb[(offset)+0].OC.type;                                    \
    assert(cod.IC.H==UNASSIGNED_CODING && xtb[(offset)+0].OC.L!=UNASSIGNED_CODING); \
    assert(cod.IC.L==UNASSIGNED_CODING && xtb[(offset)+1].OC.L!=UNASSIGNED_CODING); \
    cod.IC.H = xtb[(offset)+0].OC.L;                                          \
    cod.IC.L = xtb[(offset)+1].OC.L;                                          }


// 
// Assembles 8bit number (BYTE / unsigned char) from bit representation in column vector. LSB first
//
#define ColBinaryVectorToByte(src,i,j) (                    \
                  ((src[(i)+0][(j)] == 1) ? 1<<0 : 0)       \
                | ((src[(i)+1][(j)] == 1) ? 1<<1 : 0)       \
                | ((src[(i)+2][(j)] == 1) ? 1<<2 : 0)       \
                | ((src[(i)+3][(j)] == 1) ? 1<<3 : 0)       \
                | ((src[(i)+4][(j)] == 1) ? 1<<4 : 0)       \
                | ((src[(i)+5][(j)] == 1) ? 1<<5 : 0)       \
                | ((src[(i)+6][(j)] == 1) ? 1<<6 : 0)       \
                | ((src[(i)+7][(j)] == 1) ? 1<<7 : 0))       

//
// Takes 8bit number (BYTE / unsigned char) and stores its bit representation to col vector
// starting at given coordinates to array (may be mat_GF2). LSB first
#define ByteToColBinaryVector(c,dst,i,j) {                                  \
                dst[(i)+0][(j)] = ((c) & 1<<0) ? 1:0; dst[(i)+1][(j)] = ((c) & 1<<1) ? 1:0; \
                dst[(i)+2][(j)] = ((c) & 1<<2) ? 1:0; dst[(i)+3][(j)] = ((c) & 1<<3) ? 1:0; \
                dst[(i)+4][(j)] = ((c) & 1<<4) ? 1:0; dst[(i)+5][(j)] = ((c) & 1<<5) ? 1:0; \
                dst[(i)+6][(j)] = ((c) & 1<<6) ? 1:0; dst[(i)+7][(j)] = ((c) & 1<<7) ? 1:0;}

// Positive modulo
#define POS_MOD(a,m) (((a) % (m)) < 0 ? ((a) % (m)) + (m) : (a) % (m))

class WBAESGenerator {
public:
	WBAESGenerator();
	virtual ~WBAESGenerator();

    // Effect of shift rows operation for L bijections in T3 tables.
    // DEF:  shiftRowsLBijection[i] = to which T2 table in next round
    //       will be i-th byte of state passed as input from this round.
    //
    // Recall that T2 boxes are indexed by columns, so in first column there
    // are boxes T2_0, T2_1, T2_2, T2_3. But state array is indexed by rows (note: not by design).
    //
    // With this information we can construct L^r OUT bijection in T3 tables
    // to match L^{r+1, -1} IN bijection in T2 tables in next round.
    //
    // Every round operates on state byte in this way 
    // Upper row - which byte is selected from state array to 1,2,3,4-th column
    //		(separated by "|")
    //
    // Lower row - which byte is stored
    //
    // 00 05 10 15 | 01 06 11 12 | 02 07 08 13 | 03 04 09 14   |
    // -----------------------------------------------------   v
    // 00 04 08 12 | 01 05 09 13 | 02 06 10 14 | 03 07 11 15
    //
    // Equals with:
    //                                   +------------------------- ShiftRows() in next round
    //               +-------------------|------------------------- L(T2(ShiftRows()))
    //               |                   |                   +----- Corresponding Tboxes
    //               |                   |                   |
    //  00 01 02 03  |  00' 01' 02' 03'  |  00' 01' 02' 03'  |  00 04 08 12
    //  04 05 06 07  |  05' 06' 07' 04'  |  06' 07' 04' 05'  |  01 05 09 13
    //  08 09 10 11  |  10' 11' 08' 09'  |  08' 09' 10' 11'  |  02 06 10 14
    //  12 13 14 15  |  15' 12' 13' 14'  |  14' 15' 12' 13'  |  03 07 11 15
    //                         |
    //                         +----------------------------------- Will feed next T2 box
    //                         |
    //                  00  04  08  12     |
    //                  13  01  05  09     |  = ShiftRowsInv(CorrespondingTboxes)
    //                  10  14  02  07     |
    //                  07  11  15  03     |
    //
    // For example 00',06',08',14' will be feed to T2_0,1,2,3 boxes in new round.
    //
    // Note: actual table is transposed since Tboxes are indexed by cols.
    //
    // In next round, shift rows operation will be used again, so 
    // this table gives prescript how input bytes will be
    // mapped to T_boxes in next round (upper row), counting with 
    // shift operation in the beginning of the next round.
    static int shiftRowsLBijection[N_BYTES];
    
    // Same principle as previous = ShiftRows(CorrespondingTboxes)
    static int shiftRowsLBijectionInv[N_BYTES];

    // How shift rows affects state array - indexes.
    // Selector to state array in the beginning of encryption round
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

    //
    // 40 Generic AES instances to generate resulting cipher.
    // 10 for each round times 4 in each "section" (meaning mix column stripe)
    //
    // If all initialized to default AES, you will get default WBAES, otherwise
    // you will get cipher using dual AES - should raise known attack to high complexities.
    GenericAES AESCipher[N_ROUNDS * N_SECTIONS];
    inline GenericAES& getAESCipher(int idx){ return this->AESCipher[idx]; };

    // use given protection or not?
    bool useDualAESARelationsIdentity;
    bool useDualAESIdentity;
    bool useDualAESSimpeAlternate;
    bool useIO04x04Identity;
    bool useIO08x08Identity;
    bool useMB08x08Identity;
    bool useMB32x32Identity;

    //
    // Mixing bijections
    // Round 2..10, 16x 08x08 MB (L)
    // Round 1..9,   4x 32x32 MB (MB for each MixColumn stripe)
    //
    MB08x08_TABLE MB_L08x08 [MB_CNT_08x08_ROUNDS][MB_CNT_08x08_PER_ROUND];
    MB32x32_TABLE MB_MB32x32[MB_CNT_32x32_ROUNDS][MB_CNT_32x32_PER_ROUND];

    //
    // Input output coding - for each byte of state array.
    // It is necessary to know this mappings to use ciphers, since cipher assumes plaintext/ciphertext is encoded
    // using this bijections.

    //
    // Coding map generated for lookup table network
    WBACR_AES_CODING_MAP *codingMap;
    CODING4X4_TABLE      *pCoding04x04;
    CODING8X8_TABLE      *pCoding08x08;

    //
    // Input/output bijection encoding
    // There are 4 layers of XOR tables.
    //    1 layer - produces XOR result of T2 table output
    //    2 layer - produces resulting XOR result from 4 state bytes
    //    3 layer - XOR result from T3 tables
    //    4 layer - XOR result from all 4 T3 tables
    //
    // Default input/output encoding size = 4bits
    // List of all encodings in one generic round, in one section/MC strip.
    //    1. XOR4 -> T2           8 x 4    (round boundary, from previous round; in first round - type 1 table instead of XOR4)
    //    2. T2   -> XOR1     4 x 8 x 4
    //    3. XOR1 -> XOR2     2 x 8 x 4
    //    4. XOR2 -> T3           8 x 4
    //    5. T3   -> XOR3     4 x 8 x 4
    //    6. XOR3 -> XOR4     2 x 8 x 4
    //    7. XOR4 -> T2           8 x 4   (round boundary, to next round)
    // -----------------------------------------------------------------------------
    //                       15 x 8 x 4 = 480 IO tables in 1 MC stripe
    //                                  = 1920 in one round, 19200 in whole AES
    //
    // Encryption and decryption has same set of tables - we can use same procedure.
    //
    void generateCodingMap(WBACR_AES_CODING_MAP* pCodingMap, int *codingCount, bool encrypt);

	//
	// Generate random mixing bijections and their inverses
	// Initializes:
	//      MB_L32x32 - 8x8 bit mixing bijection (invertible matrix), with 4x4 submatrices with full rank
	//      MB_MB08x08 - 32x32 bit mixing bijection (invertible matrix), with 4x4 submatrices with full rank
	int generateMixingBijections(
			MB08x08_TABLE L08x08[][MB_CNT_08x08_PER_ROUND], int L08x08rounds,
			MB32x32_TABLE MB32x32[][MB_CNT_32x32_PER_ROUND], int MB32x32rounds,
			bool MB08x08Identity=false, bool MB32x32Identity=false);
	int generateMixingBijections(bool identity=false);

	// generates new externalEncoding
	// takes flags WBAESGEN_EXTGEN_* determining which parts of ExtEncoding will be id.
	void generateExtEncoding(ExtEncoding * extc, int flags);

	//
	// Generate WB AES tables for encryption or decryption - MAIN method here
	void generateTables(BYTE *key, enum keySize ksize, WBAES * genAES, ExtEncoding * extc, bool encrypt);

	//
	// Helper method, generates only T1 tables from external encoding
	void generateT1Tables(WBAES * genAES, ExtEncoding * extc, bool encrypt);

	//
	// Helper method, generates XOR table
	void generateXorTable(CODING * xorCoding, XTB * xtb);

	//
	// Applies external encoding - after this, state can be passed to WB AES using this external encoding
	void applyExternalEnc(W128b& state, ExtEncoding * extc, bool input);
	void applyExternalEnc(BYTE* state, ExtEncoding * extc, bool input, size_t numBlocks = 1);

	//
	// Raw method for generating random bijections
	int generate4X4Bijections(CODING4X4_TABLE * tbl, size_t size, bool identity=false);
	int generate8X8Bijections(CODING8X8_TABLE * tbl, size_t size, bool identity=false);
	int generate4X4Bijection(BIJECT4X4 *biject, BIJECT4X4 *invBiject, bool identity=false);
	int generate8X8Bijection(BIJECT8X8 *biject, BIJECT8X8 *invBiject, bool identity=false);
 	
	// test whitebox implementation with test vectors
	int testWithVectors(bool coutOutput, WBAES * genAES, int extCodingFlags = 0);
	int testComputedVectors(bool coutOutput, WBAES * genAES, ExtEncoding * extc);

	inline void BYTEArr_to_vec_GF2E(const BYTE * arr, size_t len, NTL::vec_GF2E& dst){
		charArr_to_vec_GF2E(arr, len, dst);
 	}

 	// Converts column of 8 binary values to BYTE value
    inline BYTE matGF2_to_BYTE(NTL::mat_GF2& src, int row, int col){
    	return ColBinaryVectorToByte(src, row, col);
    }
    
    // Converts BYTE value to matGF
    inline void BYTE_to_matGF2(BYTE c, NTL::mat_GF2& ret, int row, int col){
    	ByteToColBinaryVector(c, ret, row, col);
    }
 	
 	// Converts column of 32 binary values to W32b value
 	inline void matGF2_to_W32b(NTL::mat_GF2& src, int row, int col, W32b& dst){
 		//assert((src.NumRows()) < (row*8));
 		//assert((src.NumCols()) < col);
 		dst.l = 0;
 		dst.B[0] = ColBinaryVectorToByte(src, row+8*0, col);
 		dst.B[1] = ColBinaryVectorToByte(src, row+8*1, col);
 		dst.B[2] = ColBinaryVectorToByte(src, row+8*2, col);
 		dst.B[3] = ColBinaryVectorToByte(src, row+8*3, col);
 	}
 	
 	// Converts W32b file to matGF2
 	inline void W32b_to_matGF2(W32b& src, NTL::mat_GF2& dst){
        ByteToColBinaryVector(src.B[0], dst, 8*0, 0);
        ByteToColBinaryVector(src.B[0], dst, 8*1, 0);
        ByteToColBinaryVector(src.B[0], dst, 8*2, 0);
        ByteToColBinaryVector(src.B[0], dst, 8*3, 0); 
    }
    
    inline BYTE iocoding_encode08x08(BYTE src, HIGHLOW& hl, bool inverse, CODING4X4_TABLE* tbl4, CODING8X8_TABLE* tbl8){
        if (hl.type == COD_BITS_4){
            return inverse ?
                  HILO(
                	USE_IDENTITY_CODING(hl.H) ? HI(src) : tbl4[hl.H].invCoding[HI(src)],
                	USE_IDENTITY_CODING(hl.L) ? LO(src) : tbl4[hl.L].invCoding[LO(src)])

                : HILO(
                	USE_IDENTITY_CODING(hl.H) ? HI(src) : tbl4[hl.H].coding[HI(src)],
                	USE_IDENTITY_CODING(hl.L) ? LO(src) : tbl4[hl.L].coding[LO(src)]);
        } else if (hl.type == COD_BITS_8){
        	assert(tbl8 != NULL);
            return inverse ?
                  (USE_IDENTITY_CODING(hl.L) ? src : tbl8[hl.L].invCoding[src])
                : (USE_IDENTITY_CODING(hl.L) ? src : tbl8[hl.L].coding[src]);
        }
        
        return src; 
    }
    
    inline BYTE iocoding_encode08x08(BYTE src, CODING& coding, bool encodeInput, CODING4X4_TABLE* tbl4, CODING8X8_TABLE* tbl8){
	    HIGHLOW * hl = encodeInput ? &(coding.IC) : &(coding.OC);
	    return iocoding_encode08x08(src, *hl, encodeInput, tbl4, tbl8);
    }
    
    inline void iocoding_encode32x32(W32b& dst, W32b& src, W08x32Coding& coding, bool encodeInput, CODING4X4_TABLE* tbl4, CODING8X8_TABLE* tbl8){
    	// encoding input - special case, input is just 8bit wide
		if (encodeInput){
			dst.B[0] = iocoding_encode08x08(src.B[0], coding.IC, encodeInput, tbl4, tbl8);
			dst.B[1] = iocoding_encode08x08(src.B[1], coding.IC, encodeInput, tbl4, tbl8);
			dst.B[2] = iocoding_encode08x08(src.B[2], coding.IC, encodeInput, tbl4, tbl8);
			dst.B[3] = iocoding_encode08x08(src.B[3], coding.IC, encodeInput, tbl4, tbl8);
		} else {
			dst.B[0] = iocoding_encode08x08(src.B[0], coding.OC[0], encodeInput, tbl4, tbl8);
			dst.B[1] = iocoding_encode08x08(src.B[1], coding.OC[1], encodeInput, tbl4, tbl8);
			dst.B[2] = iocoding_encode08x08(src.B[2], coding.OC[2], encodeInput, tbl4, tbl8);
			dst.B[3] = iocoding_encode08x08(src.B[3], coding.OC[3], encodeInput, tbl4, tbl8);
		}
    }

    inline void iocoding_encode128x128(W128b& dst, W128b& src, W08x128Coding& coding, bool encodeInput, CODING4X4_TABLE* tbl4, CODING8X8_TABLE* tbl8){
		// encoding input - special case, input is just 8bit wide
		if (encodeInput){
			for(int i=0; i<16; i++){
				dst.B[i] = iocoding_encode08x08(src.B[i], coding.IC, encodeInput, tbl4, tbl8);
			}
		} else {
			for(int i=0; i<16; i++){
				dst.B[i] = iocoding_encode08x08(src.B[i], coding.OC[i], encodeInput, tbl4, tbl8);
			}
		}
	}

	int save(const char * filename, WBAES * aes, ExtEncoding * extCoding);
	int load(const char * filename, WBAES * aes, ExtEncoding * extCoding);
	int save(ostream& out, WBAES * aes, ExtEncoding * extCoding);
	int load(istream& ins, WBAES * aes, ExtEncoding * extCoding);
};

#ifdef WBAES_BOOST_SERIALIZATION
// serialization functions
// CODING4X4_TABLE
namespace boost{ namespace serialization {
		template<class Archive> inline void serialize(Archive &ar, struct _CODING4X4_TABLE &i, const unsigned version){
			ar & i.coding;
			ar & i.invCoding;
		}}}

// MB_TABLE
namespace boost{ namespace serialization {
		template<class Archive> inline void serialize(Archive &ar, struct _MB_TABLE &i, const unsigned version){
			ar & i.mb;
			ar & i.inv;
		}}}

// ExtEncoding
namespace boost{ namespace serialization {
		template<class Archive> inline void serialize(Archive &ar, struct _ExtEncoding &i, const unsigned version){
			ar & i.IODM[0];
			ar & i.IODM[1];
			for(int k=0; k<2; k++){
				for(int l=0; l<2*N_BYTES; l++){
					ar & i.lfC[k][l];
				}
			}
		}}}

// NTL::GF2
namespace boost { namespace serialization {
		template<class Archive> inline void serialize(Archive & ar, NTL::GF2 & t, const unsigned int file_version){
			split_free(ar, t, file_version);
		}
		template<class Archive> void save(Archive & ar, const NTL::GF2 & t, unsigned int version)
		{
			char c = (char)rep(t);
			ar & c; // rep returns long, space optimization, store as 1B (GF2 is boolean)
		}
		template<class Archive> void load(Archive & ar, NTL::GF2 & t, unsigned int version)
		{
			char cur = 0;
			ar & cur;
			t = (long)cur;
		}
	}}

// NTL::mat_GF2
namespace boost { namespace serialization {
		template<class Archive> inline void serialize(Archive & ar, NTL::mat_GF2 & t, const unsigned int file_version){
			split_free(ar, t, file_version);
		}
		template<class Archive> void save(Archive & ar, const NTL::mat_GF2 & t, unsigned int version)
		{
			long i, j, n = t.NumRows(), m = t.NumCols();
			ar & n;
			ar & m;

			// Per-element serialization. Not very space-effective. Trivial implementation.
			for(i=0; i<n; i++){
				for(j=0; j<m; j++){
					ar & t.get(i, j);
				}
			}
		}
		template<class Archive> void load(Archive & ar, NTL::mat_GF2 & t, unsigned int version)
		{
			long i, j, n, m;
			NTL::GF2 cur;

			ar & n;
			ar & m;
			t.SetDims(n, m);
			for(i=0; i<n; i++){
				for(j=0; j<m; j++){
					ar & cur;
					t.put(i, j, cur);
				}
			}
		}
	}}

BOOST_CLASS_IMPLEMENTATION(struct _CODING4X4_TABLE, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(struct _MB_TABLE, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(struct _ExtEncoding, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(NTL::GF2, boost::serialization::object_serializable);
BOOST_CLASS_IMPLEMENTATION(NTL::mat_GF2, boost::serialization::object_serializable);
#endif

#endif /* WBAESGENERATOR_H_ */
