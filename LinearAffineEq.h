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
#include "base.h"
#include "NTLUtils.h"
#include "MixingBijections.h"

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
#include <set>
#include <map>
#include <deque>
#include <iostream>

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
	mat_GF2 Tainv;
	mat_GF2 Tbinv;
	smap TaV;
	smap TbV;
	smap TbinvV;
} linearEquiv_t;

typedef deque<linearEquiv_t> linearEquivalencesList;

typedef struct _affineEquiv_t {
	linearEquiv_t linPart;  // results from linear equivalence algorithm
	bsetElem a;				// affine constant a
	bsetElem b;				// affine constant b
	smap L1;				// actual affine mapping L1
	smap L2;				// actual affine mapping L2
	bool checkPassed;		// was final check successfull?
	std::string totalHash;	// md5(L1||L2)
} affineEquiv_t;

typedef deque<affineEquiv_t> affineEquivalencesList;


// Recursive stack structure - state stored in one virtual
// recursive call, corresponds to guessing A(x) value part of code.
typedef struct _linEqGuess_t {
	bsetElem guessKey;		// x, for which was value guessed
	bsetElem guessVal;		// A(x), guessed mapping value for x
	unsigned int idx;		// index in randomly permuted guess array - to exhaust all possibilities
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
	static void dumpMapS(ostream& out, const smap& mp, bool newline);
	static void dumpSetS(ostream& out, const bset s);
	static void dumpMap(const smap& mp);
	static void dumpSet(const bset s);
	static std::string dumpMapT(const smap& mp, bool newline);
	static std::string dumpSetT(const bset s);
	static inline std::string hashSmap(const smap & map){
		std::string inputBuffer = LinearAffineEq::dumpMapT(map, false);
		return hashString(inputBuffer);
	}

	// extract linearly independent vectors from input map - keys
	bset extractLinearlyIndependent(const smap& mp);
	// converts set of vectors to matrix, vector is in column repr.
	static mat_GF2 vectorSet2GF2matrix(const bset & s, int dim);
	// converts output values of mapping smap given by values bset to GF2 matrix
	static mat_GF2 values2GF2matrix(const bset & s, const smap & m, int dim);

	/**
	 * Linear Equivalence algorithm
	 */
	int findLinearEquivalences(bsetElem * S1,   bsetElem * S1inv,
							   bsetElem * S2,   bsetElem * S2inv,
							   linearEquivalencesList * list);

	/**
	 * Affine equivalence algorithm - naive one, just using all possible affine constants with linear equivalences
	 */
	int findAffineEquivalences(bsetElem * S1,   bsetElem * S1inv,
			   	   	   	   	   bsetElem * S2,   bsetElem * S2inv,
			   	   	   	   	   affineEquivalencesList * list, bool inverseAffineConsts=false,
			   	   	   	   	   int (*callback) (affineEquiv_t *, affineEquivalencesList *, boost::unordered_set<std::string> *, LinearAffineEq *, void *) = NULL,
			   	   	   	   	   void * usrData = NULL);

	/**
	 * Builds lookup table for transformation given as matrix multiplication and compares with
	 * already obtained mapping values in mapA.
	 *
	 * Computes affine transformation: Ta * x + cst
	 *
	 * If everything is OK, 0 is returned, if at least one value mismatches, negative number is returned.
	 */
	int buildLookupTableAndCheck(mat_GF2 & Ta, bsetElem cst, smap & mapA);

	inline void setDimension(unsigned int dim){
		this->dim  = dim;
		this->size = (int) pow(2.0, (double)dim);
	}

	inline unsigned int getDim(){ return this->dim; }
	inline unsigned int getSize(){ return this->size; }
	unsigned int verbosity;
	unsigned int verbosityAffine;
	bool randomizeXGuess;


	//
	// Result attributes
	//
	int relationsCount;

protected:
	unsigned int dim;
	unsigned int size; // 2^dim

	/**
	 * Helper function for Linear Equivalence algorithm.
	 * Corresponds to part :
	 * 	If B is invertible linear mapping then
	 * 	   Derive A and check A,B at al points, that are still left in Ua and Ub
	 */
	int checkInvertibleLinear(const bset & Ua,   const bset & Ub,
							  smap & mapA,       smap & mapB,
							  bsetElem * S1,     bsetElem * S1inv,
							  bsetElem * S2,     bsetElem * S2inv,
							  mat_GF2 & Ta,      mat_GF2 & Tb,
							  mat_GF2 & Tbinv,   smap & mapBinv,
							  bool AisA);
};



} /* namespace laeqv */
} /* namespace wbacr */
#endif /* LINEARAFFINEEQ_H_ */
