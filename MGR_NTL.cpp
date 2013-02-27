//============================================================================
// Name        : MGR_NTL.cpp
// Author      : Dusan Klinec (ph4r05)
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

// NTL dependencies
#include <NTL/mat_GF2.h>
#include <NTL/vec_long.h>
#include <NTL/new.h>
#include <math.h>
NTL_CLIENT

// field size
#define FSIZE 2
// block matrix default size
#define QSIZE 4

// prototypes / forward declarations
long gaussP(mat_GF2& M, mat_GF2& P, long w);
long gaussP(mat_GF2& M, mat_GF2& P);
long generateInvertiblePM(mat_GF2& M, int p);

using namespace std;
using namespace NTL;
int main(void) {
	long i;
	puts("Hello World!!!");
	cout << "Wazzup?" <<endl;

	// very poor PRNG seeding, but just for now
	srand((unsigned)time(0));

	// sample matrix stuff
	mat_GF2 A;
	A.SetDims(QSIZE, QSIZE);
	cout << "I will show you a nice matrix: " << endl << A << endl << endl;

	// now generate invertible matrix
	i = generateInvertiblePM(A, QSIZE);
	if (i>=0){
		cout << "found invertible matrix in [" << i << "] iterations: " << endl << A << endl << endl;
	} else {
		cout << "Invertible matrix was not found" << endl;
	}

	// Now try bigger matrix - mixing bijection 8x8
	i = generateInvertiblePM(A, 8);
	if (i>=0){
		cout << "found invertible matrix in [" << i << "] iterations: " << endl << A << endl << endl;
	} else {
		cout << "Invertible matrix was not found" << endl;
	}



	return EXIT_SUCCESS;
}

/**
 * Generates random matrix of dimension pxp that is invertible in GF(2)
 */
long generateInvertiblePM(mat_GF2& M, int p){
	int rounds=0;
	long i, j;
	GF2 det;

	// Initialize M as square matrix pxp
	M.SetDims(p,p);

	// Iterate until we have some invertible matrix, or to some boundary.
	// Reaching this boundary is highly improbable for small p.
	for(rounds=0; rounds < 100; rounds++){
		// Fill matrix with random values and then compute determinant.
		for(i=0; i<p; i++){
			for(j=0; j<p; j++){
				M.put(i,j,rand()%2);
			}
		}

		// test for determinant. If determinant != 0 then matrix is non-singular, invertible
		determinant(det, M);
		if (det!=0){
			return rounds;
		}
	}

	return -1;
}

/**
 * Extended Gauss version - should return also P matrix in
 * matrix A decomposition PAQ = R where R is in canonical form.
 */
long gaussP(mat_GF2& M, mat_GF2& P, long w)
{
   long k, l;
   long i, j;
   long pos;

   long n = M.NumRows();
   long m = M.NumCols();

   if (w < 0 || w > m)
      Error("gauss: bad args");

   long wm = (m + NTL_BITS_PER_LONG - 1)/NTL_BITS_PER_LONG;

   l = 0;
   for (k = 0; k < w && l < n; k++) {
      long wk = k/NTL_BITS_PER_LONG;
      long bk = k - wk*NTL_BITS_PER_LONG;
      _ntl_ulong k_mask = 1UL << bk;


      pos = -1;
      for (i = l; i < n; i++) {
         if (M[i].rep.elts()[wk] & k_mask) {
            pos = i;
            break;
         }
      }

      if (pos != -1) {
         if (l != pos)
        	 swap(M[pos], M[l]);

         _ntl_ulong *y = M[l].rep.elts();

         for (i = l+1; i < n; i++) {
            // M[i] = M[i] + M[l]*M[i,k]

            if (M[i].rep.elts()[wk] & k_mask) {
               _ntl_ulong *x = M[i].rep.elts();

               for (j = wk; j < wm; j++)
                  x[j] ^= y[j];
            }
         }

         l++;
      }
   }

   return l;
}

long gaussP(mat_GF2& M, mat_GF2& P)
{
   return gaussP(M, P, M.NumCols());
}

