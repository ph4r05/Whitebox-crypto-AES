/*
 * WBAESGenerator.h
 *
 *  Created on: Mar 10, 2013
 *      Author: ph4r05
 */

#ifndef WBAESGENERATOR_H_
#define WBAESGENERATOR_H_

#include "WBAES.h"

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
	
};

#endif /* WBAESGENERATOR_H_ */
