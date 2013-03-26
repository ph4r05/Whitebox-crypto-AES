/*
 * LinearAffineEq.h
 *
 *  Created on: Mar 26, 2013
 *      Author: ph4r05
 */

#ifndef LINEARAFFINEEQ_H_
#define LINEARAFFINEEQ_H_

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
#include "MixingBijections.h"

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
#include <set>
#include <map>

namespace wbacr {
namespace laeqv {

//
// Basic types
//
typedef unsigned int bsetElem;
typedef boost::unordered_set<bsetElem> bset;
typedef std::map<bsetElem, bsetElem> smap;
typedef std::pair<bsetElem, bsetElem> smapElem;

typedef struct _linearEquiv_t {
	mat_GF2 Ta;
	mat_GF2 Tb;
} linearEquiv_t;

typedef struct _affineEquiv_t {
	mat_GF2 Ta;
	mat_GF2 Tb;
	GF2X a;
	GF2X b;
} affineEquiv_t;


typedef struct _linEqGuess_t {
	bsetElem guessKey;
	bsetElem guessVal;
	unsigned int idx;
	unsigned int guessCn;

	bset Ca, Cb, Ua, Ub;	// values backup for reverting bad guess
	smap mapA, mapB;		// mapping backup for reverting bad guess
} linEqGuess_t;

typedef vector<linEqGuess_t> recStack_t;

class LinearAffineEq {
public:
	LinearAffineEq();
	virtual ~LinearAffineEq();

	// returns set C = A \ B
	static bset setDiff(const bset &A, const bset &B);
	static void dumpMap(const smap& mp);
	static void dumpSet(const bset s);

	// extract linearly independent vectors from input map - keys
	static bset extractLinearlyIndependent(const smap& mp);
	// converts set of vectors to matrix, vector is in column repr.
	static mat_GF2 vectorSet2GF2matrix(const bset & s, int dim);
	// converts output values of mapping smap given by values bset to GF2 matrix
	static mat_GF2 values2GF2matrix(const bset & s, const smap & m, int dim);

	int checkInvertibleLinear(const bset & Ua,   const bset & Ub,
							  smap & mapA,      smap & mapB,
							  bsetElem * S1, 	bsetElem * S1inv,
							  bsetElem * S2, 	bsetElem * S2inv,
							  mat_GF2 & Ta,		 mat_GF2 & Tb);

	int findLinearEquivalences(bsetElem * S1, 	bsetElem * S1inv,
							bsetElem * S2, 		bsetElem * S2inv);

	inline void setDimension(unsigned int dim){
		this->dim  = dim;
		this->size = (int) pow(2.0, (double)dim);
	}

	inline unsigned int getDim(){ return this->dim; }
	inline unsigned int getSize(){ return this->size; }
	unsigned int verbosity;


	//
	// Result attributes
	//
	int relationsCount;

protected:
	unsigned int dim;
	unsigned int size; // 2^dim
};



} /* namespace laeqv */
} /* namespace wbacr */
#endif /* LINEARAFFINEEQ_H_ */
