/*
 * MixingBijections.h
 *
 *  Created on: Mar 7, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#ifndef MIXINGBIJECTIONS_H_
#define MIXINGBIJECTIONS_H_

#include "base.h"
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

// field size
#define FSIZE 2
// block matrix default size
#define QSIZE 4

/**
 * Generates random matrix of dimension pxp that is invertible in GF(2)
 */
long generateInvertiblePM(mat_GF2& M, int p);
void AddToCol(mat_GF2& x, long j, const vec_GF2& a);

/**
 * Extended Inversion version - should return also invertible P,Q matrices in
 * matrix A decomposition PAQ = R where R is in canonical form.
 *
 * Returns rank of matrix
 */
long invP(ref_GF2 d, mat_GF2& X, mat_GF2& Q, const mat_GF2& A);

/**
 * Generates n x n matrix M in canonical form for given rank.
 */
void canonical(mat_GF2& M, int rank, int n);

/**
 * Generates matrix A according to paper [http://eprint.iacr.org/2002/096.pdf]
 * From lemma 1.
 *
 * T = canonical(rank,m) + A is invertible, according to this paper.
 */
void generateARankMatrix(mat_GF2& A, int rank, int n);

/**
 * Generates mixing bijection matrix according to paper [http://eprint.iacr.org/2002/096.pdf].
 * p | t. Will compute matrix A s.t. dimension = t x t and is composed from block of size p x p
 * submatrices.
 */
int generateMixingBijection(mat_GF2& RES, int t, int p);

/**
 * Generates random bijection in set of given size.
 * Generated bijections are in form of lookup tables. size = 2^dim
 *
 * Random bijection is generated with "Knuth Shuffles" algorithm providing uniform
 * distribution of random permutations.
 *
 * See also:
 * 	Efficient Sampling of Random Permutations, http://hal.inria.fr/docs/00/11/00/56/PDF/perm.pdf
 * 	http://www.cs.tufts.edu/comp/250P/classpages/randperm.html
 * 	http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 */
int generateRandomBijection(vec_GF2X& bijection, vec_GF2X& inverse, int size, int dim);
int generateRandomBijection(unsigned char *bijection, unsigned char *inverse, int size, int init);

template<typename T> int randomPermutationT(T * bijection, int size, int init){
	int i;
	if (init>=0) {
		for(i=0;i<size;i++) bijection[i] = i+init;
	}

	// yes, we start from second element on purpose, to produce uniform distribution
	for(i=1; i<size; i++){
		// rnd is index from interval [0, i]
		int rnd = phrand() % (i+1);

		// swap values
		T tmp = bijection[rnd];
		bijection[rnd] = bijection[i];
		bijection[i] = tmp;
	}

	return 0;
}

template<typename T> int generateRandomBijectionT(T * bijection, T * inverse, int size, int init){
	int i;
	if (init) {
		for(i=0;i<size;i++){
			bijection[i] = i;
			inverse[i] = i;
		}
	}

	// yes, we start from second element on purpose, to produce uniform distribution
	for(i=1; i<size; i++){
		// rnd is index from interval [0, i]
		int rnd = phrand() % (i+1);

		// swap indexes
		T idx = inverse[bijection[rnd]];
		inverse[bijection[rnd]] = inverse[bijection[i]];
		inverse[bijection[i]]   = idx;

		// swap values
		T tmp = bijection[rnd];
		bijection[rnd] = bijection[i];
		bijection[i] = tmp;
	}

	return 0;
}

#endif /* MIXINGBIJECTIONS_H_ */
