/*
 * BGEAttack.h
 *
 *  Created on: Apr 7, 2013
 *      Author: ph4r05
 */

#ifndef BGEATTACK_H_
#define BGEATTACK_H_
#include "base.h"
#include "NTLUtils.h"
#include "MixingBijections.h"
#include "LinearAffineEq.h"
#include "GenericAES.h"
#include "WBAES.h"
#include "WBAESGenerator.h"
NTL_CLIENT
using wbacr::laeqv::affineEquiv_t;

namespace wbacr {
namespace attack {

// Function map hash -> function idx
typedef boost::unordered_map<std::string, int>fctionMap_t;
typedef std::pair<std::string, int> fctionMap_elem_t;

// Psi mapping, function idx -> [\Beta]
typedef BYTE psiFction_t[256];

// R set in algorithm for PSI recovery
typedef boost::unordered_map<BYTE, BYTE> Rmap_t;
typedef std::pair<BYTE, BYTE> Rmap_elem_t;

// Hashes GF256 function with MD5
inline std::string hashFunction(GF256_func_t f){
	std::ostringstream out;
	for(int i=0; i<GF256; i++){
		out << "["<<CHEX(i)<<"]="<<CHEX(f[i])<<",";
	}

	return hashString(out.str());
}

// composes two functions and return hash of resulting function. Computes (f \ocirc g)(x)
std::string composeFunction(GF256_func_t f, GF256_func_t g);

// Element of Sset - set of functions;
// S = {Q \op \oplus_\beta \op Q^{-1}}_{\beta \in GF(2^8)}
//
// or just ordinary function of Rbox, x2=0, x3=0
typedef struct fction_t_ {
	// Rbox function
	// f(x, c1, 0, 0)

	BYTE c1;
	GF256_func_t f;
	GF256_func_t finv;
	std::string hash;

	fction_t_(void) {
		c1=0;

		for(int i=0; i<GF256; i++){
			f[i]    = i;
			finv[i] = i;
		}

		initHash();
	}

	void initHash(void){
		hash=hashFunction(f);
	}

} fction_t;

// S set of functions.
// Note: all functions in this set should form vector space, thus
// by composing any two of them one should obtain another from set
typedef struct Sset_t_ {
	// Each function in array below is constructed as
	//    f_c1(f^{-1}_00(x))
	//  = f(f^{-1}(x,0,0,0), c1, 0, 0)
	fction_t_ fctions[GF256];
	fctionMap_t fmap;
	fction_t f_00;

	// PSI mapping (isomorphism S -> GF(2)^8)
	psiFction_t psi;
	int a;
} Sset_t;

// set S per each round (each for one state array element)
typedef struct SsetPerRound_t_{
	//
	//       +------- Column (entering to MixCol operation)
	//       |  +---- Row (y0, y1, y2, y3)
	//       |  |
	Sset_t S[4][4];
} SsetPerRound_t;

// Q~
typedef struct Qaffine_t_{
	//         +--------------- 9 Rounds with usable output encoding
	//         |      +-------- For each state array element (y0, y1, y2, y3, for each column)
	//         |      |
	fction_t Q[9][AES_BYTES];
}Qaffine_t;

// set for beta coefficients from 3.3 section
typedef boost::unordered_set<BYTE> Bset;

// characteristic polynomial multimap for resolving to beta
typedef boost::unordered_multimap<int, BYTE> charPolynomialMultimap_t;
typedef std::pair<int, BYTE> charPolynomialMultimap_elem_t;

// square matrix 8x8 (maximum) of GF2X elements
typedef struct mat_GF2X_t_ {
	NTL::GF2X x[8][8];
	int n;
}mat_GF2X_t;

// Affine equivalence
typedef std::map<BYTE, BYTE> lmap_t;
typedef std::pair<BYTE, BYTE> lmap_elem_t;
typedef struct affineEquiv_t_ {
	mat_GF2 L;
	mat_GF2 Linv;
	lmap_t  Lm;
	lmap_t  Lminv;
	BYTE   c;
} affineEquiv_t;

// vector of base vectors - used for solving singular system of homogeouns equations
typedef std::vector<NTL::vec_GF2> baseVectors_t;

// Structure for 1 affine relation in proposition 3
typedef struct prop3affineRelation_t_{
	BYTE delta;             // delta_i according to the paper, defining affine mapping
	BYTE c;                 // c_i     according to the paper, defining affine mapping

	BYTE affineConst;       // affine constant of final affine mapping determined
	mat_GF2 L;				// linear part of the mapping
	mat_GF2 Linv;			// inversion of linear part of the mapping
	GF256_func_t affineMap;	// resulting affine mapping P~_i
} prop3affineRelation_t;

// Return structure for proposition 3
typedef struct prop3struct_t_{
	prop3affineRelation_t P[4];
} prop3struct_t;

class BGEAttack {
public:
	BGEAttack();
	virtual ~BGEAttack();

	// AES to attack on
	WBAES wbaes;

	int run(void);
	void Rbox(W128b& state, bool encrypt=true, int r=1, bool noShift=false);
	void recoverPsi(Sset_t & set);
	int deriveBset(Bset & bset, GenericAES & aes, bool permissive=true);

	// Computes characteristic polynomial of matrix with coefficients from GF2.
	// Uses recursive coMatrices algorithm.
	GF2X characteristicPolynomial(mat_GF2 m);

	// Proposition 1 solver - finding affine equivalences between (yi, yj)
	// yi is computed with Rbox() function using current wbaes.
	//
	// @param r       round to compute on
	// @param col     column in given round
	// @param syi     1..4, which row to use for yi
	// @param syj     1..4, which row to use for yj
	// @param sx      which x 1..4 should be used to compute yi, yj
	int proposition1(affineEquiv_t & ret, int r, int col, int syi, int syj, int sx);

	// Proposition 2 solver - Finding A~0 for given matrix L and Beta.
	// L = A0 * beta * A0^{-1}
	// A~0 = A0 * gamma 		for some unique gamma \in GF(2^8)
	//
	// Function returns result in out variable.
	// If a system of equations obtained is non-singular, out will contain only
	// one vector - solution of the system, but can be trivial (null).
	//   If the only one found solution is trivial, return value = -1
	//   If found solution is invalid, return value = -2  (something wrong, maybe bug)
	//   Otherwise return value = 0
	//
	// If the system is singular, the out contains the set of orthogonal base vectors
	// that span subspace that solves the system. Return value = 1.
	int proposition2(mat_GF2 & inp, baseVectors_t & out, mat_GF2 beta);

	// Proposition 3 according to the paper.
	// Finding affine relations with embedded round key
	//
	// @param aes     AES to use for Sbox inverse
	// @param r       computing prop3 for given round
	// @param col     column of state array to compute on
	// @param row     y0...y3 to compute on
	// @param Atild   result from proposition 2, A~ relation
	// Returns:
	//                -1 in case Atild is empty
	//                -2 in case A~Inv cannot be constructed from Atild
	//                Flags determining whether P0,P1,P2,P3 were computed successfully (bit 1 if yes)
	int proposition3(prop3struct_t * out, GenericAES & aes, int r, int col, int row, baseVectors_t & Atild);

	// just identity on 16 elements - used when shift rows operation ignored
	static int shiftIdentity[16];
	static int shiftT2[16];

protected:
	// used internally for comptuting characteristic polynomial
	GF2X characteristicPolynomial(mat_GF2X_t m);
};

} /* namespace attack */
} /* namespace wbacr */
#endif /* BGEATTACK_H_ */
