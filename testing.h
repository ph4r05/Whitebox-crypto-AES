/*
 * MGR_NTL.h
 *
 *  Created on: Mar 29, 2013
 *      Author: ph4r05
 */

#ifndef MGR_NTL_H_
#define MGR_NTL_H_

//#include <stdio.h>
//#include <stdlib.h>
#include "base.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <set>
#include <map>

// NTL dependencies
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_long.h>
#include <cmath>
#include <vector>

#include "GenericAES.h"
#include "NTLUtils.h"
#include "MixingBijections.h"
#include "WBAES.h"
#include "WBAESGenerator.h"
#include "md5.h"
#include "LinearAffineEq.h"
#include "BGEAttack.h"

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>

// Structure for A1 A2 relations generator
typedef struct _a1a2rec {
    unsigned long int id;
    int count;
    std::vector<unsigned long int> * vec;
} a1a2rec;
typedef boost::unordered_map<std::string, a1a2rec> a1a2map;

// Structure for AES affine rel. finding
typedef struct _AESAffineElement {
	int square, multi, type;
} AESAffineElement;
typedef boost::unordered_multimap<std::string, AESAffineElement> AESAffineMap;
typedef std::pair<std::string, AESAffineElement> AESAffineMapElem;

std::string hashLookupTable(NTL::vec_GF2E s);
int MBgen(void);
int dualAESTest(void);
int A1A2relationsGenerator(void);
int AESAffineRelationsVerify(bool inverseSbox = true);
int AffCallbackCorrespondence(  wbacr::laeqv::affineEquiv_t * el, 
                                wbacr::laeqv::affineEquivalencesList * lish, 
                                boost::unordered_set<std::string> * hashes, 
                                wbacr::laeqv::LinearAffineEq * eqCheck, 
                                void * usrData);

#endif /* MGR_NTL_H_ */
