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

// NTL dependencies
#include "NTLUtils.h"
#include <iomanip>
NTL_CLIENT

int WBAESGenerator::shiftRowsLBijection[N_BYTES] = {
		0, 13, 10, 7,
		4,  1, 14, 11,
		8,  5,  2, 15,
		12, 9,  6,  3
};

WBAESGenerator::WBAESGenerator() {
	// TODO Auto-generated constructor stub

}

WBAESGenerator::~WBAESGenerator() {
	// TODO Auto-generated destructor stub
}

