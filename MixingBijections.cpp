/*
 * MixingBijections.cpp
 *
 *  Created on: Mar 7, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#include "MixingBijections.h"
#include <iostream>
NTL_CLIENT

using namespace std;
using namespace NTL;
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
				M.put(i,j,phrand()%2);
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
 * Extended Inversion version - should return also invertible P,Q matrices in
 * matrix A decomposition PAQ = R where R is in canonical form.
 *
 * Returns rank of matrix
 */
long invP(ref_GF2 d, mat_GF2& X, mat_GF2& Q, const mat_GF2& A)
{
   long n = A.NumRows();
   if (A.NumCols() != n)
      Error("solve: nonsquare matrix");

   if (n == 0) {
      X.SetDims(0, 0);
      set(d);
   }

   long i, j, k, pos;
   long rank=n;

   //
   // Gauss Jordan Elimination. Matrix M is extended version
   // with 2*N columns. Copy of A is in the left half, unity
   // matrix is in the right half.
   //
   mat_GF2 M;
   M.SetDims(n, 2*n);

   vec_GF2 aa;
   aa.SetLength(2*n);

   // Initializing Q matrix as unit matrix, will correspond to
   // column operations performed to obtain canonical form.
   // Since matrix is represented as array of vectors (rows),
   // we will work with transpose version of matrix.
   mat_GF2 I;
   ident(I, n);
   transpose(Q, I);
   for (i = 0; i < n; i++) {
      aa = A[i];
      aa.SetLength(2*n);
      aa.put(n+i, 1);
      M[i] = aa;
   }

   long wn = ((2*n) + NTL_BITS_PER_LONG - 1)/NTL_BITS_PER_LONG;
   for (k = 0; k < n; k++) {
      long wk = k/NTL_BITS_PER_LONG;
      long bk = k - wk*NTL_BITS_PER_LONG;
      _ntl_ulong k_mask = 1UL << bk;

#ifdef DEBUGOUT
      cout << "Intermediate result in step=" << k <<  "; Matrix" << endl << M << endl;
#endif
      // Find leading one in rows on k-th position in row.
      // Search in interval [k,n] (thus from current row down).
      pos = -1;
      for (i = k; i < n; i++) {
         if (M[i].rep.elts()[wk] & k_mask) {
            pos = i;
            break;
         }
      }
#ifdef DEBUGOUT
      cout << "Line pos: [" << pos << "] has leading 1 on [" << k << "]. position" << endl;
#endif
      if (pos == -1) {
		  // If here it means there is no row in matrix that has leading
    	  // 1 on k-th position.
    	  //
    	  // Thus now look in rows [k,n] and find some row that has
    	  // 1 element on position > k. Then we will perform column swap
    	  // to obtain 1 element on desired position = k. This change has to be
    	  // reflected to Q matrix.
    	  //
    	  // Finding unit element on position k+1 in all rows [k,n].
    	  // If fails, look for unit element on position k+2 in all rows [k,n]...
    	  int kk, ii, colpos=-1;
			for (kk = k+1; kk < n; kk++) {
				long kwk = kk / NTL_BITS_PER_LONG;
				long kbk = kk - kwk * NTL_BITS_PER_LONG;
				_ntl_ulong kk_mask = 1UL << kbk;
				colpos=kk;

				// Find leading one in rows on kk-th position in row.
				// Search in interval [k,n] (thus from current row down).
#ifdef DEBUGOUT
				cout << "Looking for leading 1 element in column: " << kk << "; mask: " << kk_mask << endl;
#endif
				pos = -1;
				for (ii = k; ii < n; ii++) {
					if (M[ii].rep.elts()[kwk] & kk_mask) {
						pos = ii;
						break;
					}
				}
				if (pos!=-1) break;
			}


			if (pos==-1){
				// No such column exists, thus just simply null rest of columns in Q matrix
				// to obtain canonical form of product PAQ.
				rank = k;
#ifdef DEBUGOUT
				cout << "No such column exists, we are already in canonical form;"\
						"nulling all columns from: " << k << "; Rank: " << rank << endl;
#endif
				for(kk=k; kk<n; kk++){
					for(ii=0; ii<n; ii++){
						Q.put(kk, ii, 0);
					}
				}

				break;
			}

#ifdef DEBUGOUT
			cout << "Swaping column [" << k << "] with column [" << colpos << "]. Matrix: " <<endl;
#endif
			// Do column swap to obtain 1 on desired k-th position.
			for(ii=0; ii<n; ii++){
				GF2 tmp = M.get(ii, k);
				M.put(ii, k, M.get(ii, colpos));
				M.put(ii, colpos, tmp);
			}

			// reflect this swap to Q matrix, swap rows (transpose form)
			swap(Q[colpos], Q[k]);
#ifdef DEBUGOUT
			cout << M << endl << "Qmatrix: " << endl << Q << endl << endl;
#endif
		}

      if (pos != -1) {
    	  // row number <pos> has leading one on k-th position
         if (k != pos) {
#ifdef DEBUGOUT
        	 cout << "Swap line " << pos << " with line " << k << endl;
#endif
            swap(M[pos], M[k]);
         }

         // obtain bit representation of vector in i-th row
         _ntl_ulong *y = M[k].rep.elts();

         for (i = k+1; i < n; i++) {
            // M[i] = M[i] + M[k]*M[i,k]

        	 // By another words, we are re-seting other 1s
        	 // in rows > k (down).
            if (M[i].rep.elts()[wk] & k_mask) {
               _ntl_ulong *x = M[i].rep.elts();

               // simple element-by-element addition
               for (j = wk; j < wn; j++)
                  x[j] ^= y[j];
            }
         }
      }
   }

   vec_GF2 XX;
   XX.SetLength(2*n);

   X.SetDims(n, n);
   clear(X);

   for (j = 0; j < n; j++) {
      XX.SetLength(n+j+1);
      clear(XX);
      XX.put(n+j, to_GF2(1));

      for (i = n-1; i >= 0; i--) {
         XX.put(i, XX*M[i]);
      }

      XX.SetLength(n);
      AddToCol(X, j, XX);
   }

   // transpose Q matrix finally
   Q = transpose(Q);

   // determinant=0 <=> rank == n
   if (rank==n) set(d);
   else clear(d);

   return rank;
}

void AddToCol(mat_GF2& x, long j, const vec_GF2& a)
// add a to column j of x
// ALIAS RESTRICTION: a should not alias any row of x
{
   long n = x.NumRows();
   long m = x.NumCols();

   if (a.length() != n || j < 0 || j >= m)
      Error("AddToCol: bad args");

   long wj = j/NTL_BITS_PER_LONG;
   long bj = j - wj*NTL_BITS_PER_LONG;
   _ntl_ulong j_mask = 1UL << bj;

   const _ntl_ulong *ap = a.rep.elts();
   _ntl_ulong a_mask = 1;

   long i;
   for (i = 0; i < n; i++) {
      if (*ap & a_mask)
         x[i].rep.elts()[wj] ^= j_mask;

      a_mask <<= 1;
      if (!a_mask) {
         a_mask = 1;
         ap++;
      }
   }
}

/**
 * Generates n x n matrix M in canonical form for given rank.
 */
void canonical(mat_GF2& M, int rank, int n){
	long i=0;
	ident(M, n);
	for(i=rank+1; i<n; i++){
		M.put(i,i,0);
	}
}


/**
 * Generates matrix A according to paper [http://eprint.iacr.org/2002/096.pdf]
 * From lemma 1.
 *
 * T = canonical(rank,m) + A is invertible, according to this paper.
 */
void generateARankMatrix(mat_GF2& A, int rank, int n){
	long i=0, offset=0;

	A.SetDims(n,n);
	clear(A);
	if (rank==1){
		// On rank = 1 matrix has special form [1 1; 1 0] and then I
		A.put(0,0,1);
		A.put(0,1,1);
		A.put(1,0,1);
		for(i=2; i<n; i++){
			A.put(i,i,1);
		}
		return;
	}

	if ((rank % 2) == 1){
		// First block of matrix is 3x3 in special form [1 1 1; 1 1 0; 1 0 0]
		A.put(0,0,1);
		A.put(0,1,1);
		A.put(0,2,1);
		A.put(1,0,1);
		A.put(1,1,1);
		A.put(2,0,1);
		offset=3;
	}

	//
	// Merged case - r is odd or even and >= 3
	//

	// For even rank it is easy to construct
	// On diagonals is <rank> copies of matrix [0 1; 1 1]
	// filled with I_2 on rest of blocks
	for(i=0; i<rank/2; i++){
		A.put(2*i   + offset, 2*i+1 + offset, 1);
		A.put(2*i+1 + offset, 2*i   + offset, 1);
		A.put(2*i+1 + offset, 2*i+1 + offset, 1);
	}
	// the rest fill with 1 on diagonals (I_{n-r} matrix)
	for(i=rank+offset-1; i<n; i++){
		if (i<0) continue;
		A.put(i,i,1);
	}
	return;
}

/**
 * Generates mixing bijection matrix according to paper [http://eprint.iacr.org/2002/096.pdf].
 * p | t. Will compute matrix A s.t. dimension = t x t and is composed from block of size p x p
 * submatrices.
 */
int generateMixingBijection(mat_GF2& RES, int t, int p){
	// validate parameters
	if (t<p || (t%p) != 0){
		return -1;
	}
	RES.SetDims(t,t);

	// 0. generate M matrix pxp that is invertible
	mat_GF2 M;
	long res = generateInvertiblePM(M, p);
	if (res < 0) {
		// matrix was not found in 100 steps, weeeeird. HIGHLY UNPROBABLE.
		return -1;
	}
#ifdef DEBUGOUT
	cout << "generated M0 invertible matrix: " << endl << M << endl << endl;
#endif
	// some matrices that we will need, naming according to the paper
	mat_GF2 X; 	mat_GF2 Y;
	mat_GF2 P;	mat_GF2 Pinv;
	mat_GF2 Q;	mat_GF2 Qinv;
	mat_GF2 A;
	mat_GF2 TMP;
	mat_GF2 Minv;
	mat_GF2 N;
	GF2 d;
	ref_GF2 dd(d);

	int i,j,k;
	int curT = p;			// current size of matrix M
	int tmp;				// current column/row
	for(; curT < t; curT+=p){
		int pBlocksInM=curT/p;	// number of pxp sub-matrices in M

		// 1. X matrix - p x t matrix, generated from M matrix using some row
		X.SetDims(p, curT);
		tmp = phrand() % pBlocksInM;		// current row
		for(i=p*tmp,k=0; k<p; i++, k++){
			for(j=0; j<curT; j++){
				X.put(k,j, M.get(i,j));
			}
		}

		// 2. Y matrix - t x p matrix, generated from M matrix using some column
		Y.SetDims(curT, p);
		tmp = phrand() % pBlocksInM;
		for(i=0; i<curT; i++){
			for(j=p*tmp,k=0; k<p; j++, k++){
				Y.put(i,k, M.get(i,j));
			}
		}

		// 3. computing invertible P,Q matrices
		inv(Minv, M);
		TMP = X * Minv * Y;
#ifdef DEBUGOUT
		cout << "X matrix:" << endl << X << endl << endl;
		cout << "Y matrix:" << endl << Y << endl << endl;
		cout << "Generated M inverse: " << endl << Minv << endl << endl;
		cout << "TMP: " << endl << TMP << endl << endl;
#endif
		int rank = invP(dd, P, Q, TMP);
#ifdef DEBUGOUT
		cout << "Rank of TMP: " << rank;
		cout << "; P=" << endl << P << endl << endl;
		cout << "; Q=" << endl << Q << endl << endl;
#endif
		// 4. A matrix
		generateARankMatrix(A, rank, p);
#ifdef DEBUGOUT
		cout << "; A=" << endl << A << endl << endl;
#endif
		// 5. resulting matrix
		mat_GF2 TMP2;
		N.SetDims(curT + p, curT + p);
		inv(Pinv, P);
		inv(Qinv, Q);
		TMP2 = TMP + Pinv*A*Qinv;

		// copy M matrix, M is curT x curT matrix
		for(i=0;i<curT;i++){
			for(j=0;j<curT;j++){
				N.put(i,j,M.get(i,j));
			}
		}
		// copy X matrix, p x curT
		for(i=0;i<p;i++){
			for(j=0;j<curT;j++){
				N.put(curT+i,j,X.get(i,j));
			}
		}
		// copy Y matrix, curT x p
		for(i=0;i<curT;i++){
			for(j=0;j<p;j++){
				N.put(i,curT+j,Y.get(i,j));
			}
		}
		// copy TMP2 matrix, p x p
		for(i=0;i<p;i++){
			for(j=0;j<p;j++){
				N.put(curT+i,curT+j,TMP2.get(i,j));
			}
		}
#ifdef DEBUGOUT
		cout << "Intermediate result for curT=" << curT << "; Matrix = " << endl << N << endl << endl;
#endif
		M = N;
	}

	RES = M;
	return 0;
}

int generateRandomBijection(vec_GF2X& bijection, vec_GF2X& inverse, int size, int dim){
	int i;
	bijection.SetLength(size);
	inverse.SetLength(size);
	for(i=0; i<size; i++){
		bijection.put(i, GF2XFromLong(i, dim));
		inverse.put(i, GF2XFromLong(i, dim));
	}

	// yes, we start from second element on purpose, to produce uniform distribution
	for(i=1; i<size; i++){
		// rnd is index from interval [0, i]
		int rnd = phrand() % (i+1);
		swap(inverse[getLong(bijection[i])], inverse[getLong(bijection[rnd])]);
		swap(bijection[i], bijection[rnd]);
	}

	return 0;
}

int generateRandomBijection(unsigned char *bijection, unsigned char *inverse, int size, int init){
	return generateRandomBijectionT(bijection, inverse, size, init);
}
