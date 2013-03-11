/*
 * WBAESGenerator.h
 *
 *  Created on: Mar 10, 2013
 *      Author: ph4r05
 */

#ifndef WBAESGENERATOR_H_
#define WBAESGENERATOR_H_

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

//  DEFINITIONS OF STRINGS USED AS INSERTION BEGIN OF GENERATED TABLES INTO HEADER FILES  
//  INSERTION BEHAVIOUR:
//  1. SEARCH FOR '#keyword#'
//  2. #keyword#user_data#  EXPECTED 
//  2. REPLACED BY 'keyworduser_data = {...TABLE DATA...}' STRING 
#define STR_ROUNDTABLES         "#roundTables#"
#define STR_XORTABLES           "#xorTables#"
#define STR_FINALROUNDTABLE     "#finalRoundTable#"
#define STR_INVROUNDTABLES      "#invRoundTables#"
#define STR_INVXORTABLES        "#invXorTables#"
#define STR_INVFIRSTROUNDTABLE  "#invFirstRoundTable#"

#define NO_CODING           0x00000000  // IDENTITY CODING
#define UNASSIGNED_CODING   0xFFFFFFFF  // INVALID CODING
#define UNUSED_CODING       0xFFFFFFFE  // This coding is not in use (XOR tables use only lower 4 bits for result)
// VALID CODINGS ORDINARY NUMBER IS FROM 0x00000000 TO 0xFFFFFFFE (TOTAL COUNT == 2^32 - 1) 

// CODING SIZE TYPE
#define COD_BITS_UNASSIGNED 0x00
#define COD_BITS_4          0x01
#define COD_BITS_8          0x02

// MIXING BIJECTION TYPE
#define MB_IDENTITY 0x00
#define MB_8x8      0x01
#define MB_32x32    0x02

// MIXING BIJECTION COUNTS
#define MB_CNT_08x08_ROUNDS 9
#define MB_CNT_08x08_PER_ROUND 16
#define MB_CNT_32x32_ROUNDS 9
#define MB_CNT_32x32_PER_ROUND 4


#define IOBLOCK_BASEID      10000   //FIRST CODING ID ASSIGNED TO 8-BITS IO BLOCK CONDING. INTRODUCED TO STOP MISMASH BETWEN 4-BITS AND 8-BITS CODING IDs

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
//  Mixing bijection
//
typedef struct _MB_TABLE {
	//int         type;
    NTL::mat_GF2     mb;
    NTL::mat_GF2     inv;          // SPEED OPTIMALIZATION, CAN BE ALSO COMPUTED FROM coding MEMBER (DUE TO BIJECTION PROPERTY)
    
    _MB_TABLE(void) {
        //type=MB_IDENTITY;
    }   
} MB_TABLE;

typedef MB_TABLE MB08x08_TABLE;
typedef MB_TABLE MB32x32_TABLE;

//
// Coding for T2 and T3 boxes, 8bit -> 32bit
//
typedef struct _W08x32Coding {
    HIGHLOW IC;  
    HIGHLOW OC[4];
} W08x32Coding;

// NOTE: 
// Coding for XOR boxes can be done with CODING type easily
//
//  WBACR_AES_CODING_MAP DEFINES ASSIGNED CODING TO EACH VALUE OF WBACR_AES_TABLE
//
typedef CODING  CODING32W[4];
typedef struct _WBACR_AES_CODING_MAP {
    // ENCRYPT TABLES
    W08x32Coding        eT2[N_ROUNDS][N_SECTIONS][4];                 // ENCRYPT ROUND CODING MAP
    W08x32Coding        eT3[N_ROUNDS][N_SECTIONS][4];                 // ENCRYPT ROUND CODING MAP   
    CODING              eXOR1[N_ROUNDS][N_SECTIONS][24];              // 24 == 8 4-BITS PARTS OF 32-BITS DWORD
    CODING              eXOR2[N_ROUNDS][N_SECTIONS][24];              // 24 == 8 4-BITS PARTS OF 32-BITS DWORD
    CODING              eFINALROUND[N_SECTIONS][4];
    
    // DECRYPT TABLES
    CODING32W           INVROUND[N_ROUNDS][N_SECTIONS][4]; // DECRYPT ROUND    
    CODING              INVXOR[N_ROUNDS][N_SECTIONS][24];           // 24 == 8 4-BITS PARTS OF 32-BITS DWORD
    CODING              INVFIRSTROUND[N_SECTIONS][4];
} WBACR_AES_CODING_MAP, *PWBACR_AES_CODING_MAP;




#define ALLOCW08x32Coding(cod, idx) {                       \
    cod.OC[0].type = COD_BITS_4; cod.OC[1].type = COD_BITS_4; \
    cod.OC[2].type = COD_BITS_4; cod.OC[3].type = COD_BITS_4; \
    cod.OC[0].H = ++(idx);                                  \
    cod.OC[0].L = ++(idx);                                  \
    cod.OC[1].H = ++(idx);                                  \
    cod.OC[1].L = ++(idx);                                  \
    cod.OC[2].H = ++(idx);                                  \
    cod.OC[2].L = ++(idx);                                  \
    cod.OC[3].H = ++(idx);                                  \
    cod.OC[3].L = ++(idx);                                  };

//
// Connects OUTPUT coding of 32bit wide boxes (T2,T3) to INPUT coding of XOR boxes, 32bit wide. 
// Each XOR box accepts 2 arguments, first in HIGH part, second in LOW part, thus when associating
// mapping from one particular W32box we are using either HIGH or LOW parts. 
//
#define CONNECT_W08x32_TO_XOR(cod, xtb, HL, offset) {             \
    xtb[(offset)+0].IC.HL = cod.OC[0].H;                          \
    xtb[(offset)+1].IC.HL = cod.OC[0].L;                          \
    xtb[(offset)+2].IC.HL = cod.OC[1].H;                          \
    xtb[(offset)+3].IC.HL = cod.OC[1].L;                          \
    xtb[(offset)+4].IC.HL = cod.OC[2].H;                          \
    xtb[(offset)+5].IC.HL = cod.OC[2].L;                          \
    xtb[(offset)+6].IC.HL = cod.OC[3].H;                          \
    xtb[(offset)+7].IC.HL = cod.OC[3].L;                          }

#define CONNECT_W08x32_TO_XOR_H(cod, xtb, offset) CONNECT_W08x32_TO_XOR(cod, xtb, H, offset)
#define CONNECT_W08x32_TO_XOR_L(cod, xtb, offset) CONNECT_W08x32_TO_XOR(cod, xtb, L, offset)

//
// Allocates new output coding for XOR boxes. 
// Recall that output of XOR is stored in LOW part, thus upper is unused -> no allocation for upper part.
//
#define ALLOCXORCoding(xtb, offset, idx) {                                                                           \
    xtb[(offset)+0].OC.type = COD_BITS_4; xtb[(offset)+0].OC.H = UNUSED_CODING; xtb[(offset)+0].OC.L = ++(idx);      \
    xtb[(offset)+1].OC.type = COD_BITS_4; xtb[(offset)+1].OC.H = UNUSED_CODING; xtb[(offset)+1].OC.L = ++(idx);      \
    xtb[(offset)+2].OC.type = COD_BITS_4; xtb[(offset)+2].OC.H = UNUSED_CODING; xtb[(offset)+2].OC.L = ++(idx);      \
    xtb[(offset)+3].OC.type = COD_BITS_4; xtb[(offset)+3].OC.H = UNUSED_CODING; xtb[(offset)+3].OC.L = ++(idx);      \
    xtb[(offset)+4].OC.type = COD_BITS_4; xtb[(offset)+4].OC.H = UNUSED_CODING; xtb[(offset)+4].OC.L = ++(idx);      \
    xtb[(offset)+5].OC.type = COD_BITS_4; xtb[(offset)+5].OC.H = UNUSED_CODING; xtb[(offset)+5].OC.L = ++(idx);      \
    xtb[(offset)+6].OC.type = COD_BITS_4; xtb[(offset)+6].OC.H = UNUSED_CODING; xtb[(offset)+6].OC.L = ++(idx);      \
    xtb[(offset)+7].OC.type = COD_BITS_4; xtb[(offset)+7].OC.H = UNUSED_CODING; xtb[(offset)+7].OC.L = ++(idx);      };


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
    xtb3[(offset3)+0].IC.HL = xtb1[(offset1)+0].OC.L;                           \
	xtb3[(offset3)+1].IC.HL = xtb1[(offset1)+1].OC.L;                           \
    xtb3[(offset3)+2].IC.HL = xtb1[(offset1)+2].OC.L;                           \
    xtb3[(offset3)+3].IC.HL = xtb1[(offset1)+3].OC.L;                           \
    xtb3[(offset3)+4].IC.HL = xtb1[(offset1)+4].OC.L;                           \
    xtb3[(offset3)+5].IC.HL = xtb1[(offset1)+5].OC.L;                           \
    xtb3[(offset3)+6].IC.HL = xtb1[(offset1)+6].OC.L;                           \
    xtb3[(offset3)+7].IC.HL = xtb1[(offset1)+7].OC.L;                           }

#define CONNECT_XOR_TO_XOR_H(xtb1, offset1, xtb3, offset3) CONNECT_XOR_TO_XOR(xtb1, offset1, xtb3, offset3, H)
#define CONNECT_XOR_TO_XOR_L(xtb1, offset1, xtb3, offset3) CONNECT_XOR_TO_XOR(xtb1, offset1, xtb3, offset3, L)

//
// Connects 8bit output from 2 consecutive XOR tables to 8b input of W08x32 table
//
#define CONNECT_XOR_TO_W08x32(xtb, offset, cod) {                             \
    cod.IC.type = xtb[(offset)+0].OC.type;                                    \
    cod.IC.H = xtb[(offset)+0].OC.L;                                          \
    cod.IC.L = xtb[(offset)+1].OC.L;                                          }


class WBAESGenerator {
public:
	WBAESGenerator();
	virtual ~WBAESGenerator();
	
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
    
	//
	// 40 Generic AES instances to generate resulting cipher.
	// 10 for each round times 4 in each "section" (meaning mix column stripe)
	//
	// If all initialized to default AES, you will get default WBAES, otherwise
	// you will get cipher using dual AES - should raise known attack to high complexities.
	GenericAES AESCipher[N_ROUNDS * N_SECTIONS];
	
	//
	// Mixing bijections
	// Round 2..10, 16x 08x08 MB (L)
	// Round 1..9,   4x 32x32 MB (MB for each MixColumn stripe)  
    //
 	MB08x08_TABLE MB_L08x08 [MB_CNT_08x08_ROUNDS][MB_CNT_08x08_PER_ROUND];
 	MB32x32_TABLE MB_MB32x32[MB_CNT_32x32_ROUNDS][MB_CNT_32x32_PER_ROUND];
 	
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
 	// TODO: define types & storage tables, use Petr's way of associating bijections to tabes (networking)
 	//  1 | T2 | T3 | T4 | X1 | X2 | X3 | B1 | B2 | B3 | B4 | X4 | X5 | X6
 	//  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 | 12 | 13
 	void encGenerateCodingMap(WBACR_AES_CODING_MAP& pCodingMap, int *codingCount);
 	
 	//
 	// Generate random mixing bijections and their inverses
 	// Initializes:
 	//      MB_L32x32 - 8x8 bit mixing bijection (invertible matrix), with 4x4 submatrices with full rank
 	//      MB_MB08x08 - 32x32 bit mixing bijection (invertible matrix), with 4x4 submatrices with full rank
 	void generateMixingBijections();
 	
 	//
 	// Generate random coding (bijections).
 	void generateIOCoding();
};

#endif /* WBAESGENERATOR_H_ */
