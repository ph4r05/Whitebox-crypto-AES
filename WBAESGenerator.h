/*
 * WBAESGenerator.h
 *
 *  Created on: Mar 10, 2013
 *      Author: ph4r05
 */

#ifndef WBAESGENERATOR_H_
#define WBAESGENERATOR_H_

#include "WBAES.h"

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
// VALID CODINGS ORDINARY NUMBER IS FROM 0x00000000 TO 0xFFFFFFFE (TOTAL COUNT == 2^32 - 1) 

// CODING SIZE TYPE
#define COD_BITS_UNASSIGNED 0x00
#define COD_BITS_4          0x01
#define COD_BITS_8          0x02

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
	// Round 1..9,   4x 32x32 MB (MB for each MixCOlumn stripe)  
    // 	
 	// TODO: mixing bijections definition & storage tables, no mapping
 	// is needed, easy mappings.
 	
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
 	//    1. XOR4 -> T2           8 x 4    (round boundary, from previous round)
 	//    2. T2   -> XOR1     4 x 8 x 4
 	//    3. XOR1 -> XOR2     2 x 8 x 4  
 	//    4. XOR2 -> T3           8 x 4
 	//    5. T3   -> XOR3     4 x 8 x 4
 	//    6. XOR3 -> XOR4     2 x 8 x 4
 	//    7. XOR4 -> T2           8 x 4   (round boundary, to next round)
 	// 
 	// TODO: define types & storage tables, use Petr's way of associating bijections to tabes (networking)
};

#endif /* WBAESGENERATOR_H_ */
