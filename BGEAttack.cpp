/*
 * BGEAttack.cpp
 *
 *  Created on: Apr 7, 2013
 *      Author: ph4r05
 */

#include "BGEAttack.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

namespace wbacr {
namespace attack {

BGEAttack::BGEAttack() {
	this->coding = NULL;
	this->wbaes = NULL;
}

BGEAttack::~BGEAttack() {
	;
}

using namespace std;
using namespace NTL;
using namespace boost;
using namespace wbacr::attack;
using namespace wbacr::laeqv;

const int BGEAttack::shiftIdentity[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
const int BGEAttack::shiftT2[16] =
		{0,  4,  8, 12,
	    13,  1,  5,  9,
	    10, 14,  2,  6,
	     7, 11, 15,  3
		};


std::string composeFunction(GF256_func_t f, GF256_func_t g){
	unsigned int x;	// warning, if is type BYTE, then it will overflow and loop forever
	fction_t nf;
	for(x=0; x<=(GF256-1); x++){
		//cout << "; X=" << x << endl;
		//cout << "; g["<<CHEX(x)<<"] = " CHEX(g[x]) << ";" << endl;
		//cout << "; f[g["<<CHEX(x)<<"]] = " CHEX(f[g[x]]) << ";" << endl;

		nf.f[x] = f[g[x]];
		nf.finv[f[g[x]]] = x;
	}

	return hashFunction(nf.f);
}


void BGEAttack::Rbox(W128b& state, bool encrypt, int r, bool noShift, int colMask2compute){
	int i=0;
	W32b  ires[N_BYTES];				// intermediate result for T2,T3-boxes
	W128b ares[N_BYTES];				// intermediate result for T1-boxes

	// encryption/decryption dependent operations and tables
	const int (&shiftOp)[N_BYTES] = noShift ? this->shiftIdentity : (encrypt ? (this->wbaes->shiftRows)   : (this->wbaes->shiftRowsInv));
	W32XTB (&edXTab)[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS] = encrypt ? (this->wbaes->eXTab) 	         : (this->wbaes->dXTab);
	W32XTB (&edXTabEx)[2][15][4]                         = encrypt ? (this->wbaes->eXTabEx)          : (this->wbaes->dXTabEx);
	AES_TB_TYPE1 (&edTab1)[2][N_BYTES]                   = encrypt ? (this->wbaes->eTab1)            : (this->wbaes->dTab1);
	AES_TB_TYPE2 (&edTab2)[N_ROUNDS][N_BYTES]			 = encrypt ? (this->wbaes->eTab2) 	         : (this->wbaes->dTab2);
	AES_TB_TYPE3 (&edTab3)[N_ROUNDS][N_BYTES]			 = encrypt ? (this->wbaes->eTab3) 	         : (this->wbaes->dTab3);
#ifdef AES_BGE_ATTACK
	GF256_func_t (&edOutputBijection)[N_ROUNDS][N_BYTES] = encrypt ? (this->wbaes->eOutputBijection) : (this->wbaes->dOutputBijection);
#endif
	// Last round = special case. Just T1 boxes && 128-bot XOR cascade
	if (r==(N_ROUNDS-1)){
		for(i=0; i<N_BYTES; i+=4){
			W128CP(ares[i+0], edTab1[1][i+0][state.B[shiftOp[i/4+0*4]]]);
			W128CP(ares[i+1], edTab1[1][i+1][state.B[shiftOp[i/4+1*4]]]);
			W128CP(ares[i+2], edTab1[1][i+2][state.B[shiftOp[i/4+2*4]]]);
			W128CP(ares[i+3], edTab1[1][i+3][state.B[shiftOp[i/4+3*4]]]);
		}

		// and finally compute XOR cascade again, now for T1[1] - output T1
		// 1st level of XORs
		for(i=0;i<N_BYTES;i+=2){
			op8xor_128(ares[i+0], ares[i+1], edXTabEx[1][i/2], ares[i+0]);  // 1 xor 2 --> 1
		}

		// Finish XOR cascade by hand
		op8xor_128(ares[0],  ares[2],  edXTabEx[1][8],  ares[0]);  // 0  xor 2  --> 0
		op8xor_128(ares[4],  ares[6],  edXTabEx[1][9],  ares[4]);  // 4  xor 6  --> 4
		op8xor_128(ares[8],  ares[10], edXTabEx[1][10], ares[8]);  // 8  xor 10 --> 8
		op8xor_128(ares[12], ares[14], edXTabEx[1][11], ares[12]); // 12 xor 14 --> 12
		op8xor_128(ares[0],  ares[4],  edXTabEx[1][12], ares[0]);  // 0 xor 4  --> 0
		op8xor_128(ares[8],  ares[12], edXTabEx[1][13], ares[8]);  // 8 xor 12 --> 8
		op8xor_128(ares[0],  ares[8],  edXTabEx[1][14], ares[0]);  // 0 xor 8 --> 0
		for(i=0; i<N_BYTES; i++){
			state.B[i] = ares[0].B[idxTranspose(i)];
		}
	} else {
		// Firt round -> need to apply T1 boxes & 128-bit XOR cascade before
		if (r==0){
			// At first we have to put input to T1 boxes directly, no shift rows
			// compute result to ares[16]
			for(i=0; i<N_BYTES; i++){
				// Note: Tbox is indexed by cols, state by rows - transpose needed here
				W128CP(ares[i], edTab1[0][i][state.B[idxTranspose(i)]]);
			}

			// 1st level of XORs
			for(i=0;i<N_BYTES;i+=2){
				op8xor_128(ares[i+0], ares[i+1], edXTabEx[0][i/2], ares[i+0]);  // 1 xor 2 --> 1
			}

			// Finish XOR cascade
			op8xor_128(ares[0],  ares[2],  edXTabEx[0][8],  ares[0]);  // 0  xor 2  --> 0
			op8xor_128(ares[4],  ares[6],  edXTabEx[0][9],  ares[4]);  // 4  xor 6  --> 4
			op8xor_128(ares[8],  ares[10], edXTabEx[0][10], ares[8]);  // 8  xor 10 --> 8
			op8xor_128(ares[12], ares[14], edXTabEx[0][11], ares[12]); // 12 xor 14 --> 12
			op8xor_128(ares[0],  ares[4],  edXTabEx[0][12], ares[0]);  // 0 xor 4  --> 0
			op8xor_128(ares[8],  ares[12], edXTabEx[0][13], ares[8]);  // 8 xor 12 --> 8
			op8xor_128(ares[0],  ares[8],  edXTabEx[0][14], ares[0]);  // 0 xor 8 --> 0
			W128CP(state, ares[0]);
		}


		// Perform rest of the operations on 4 tuples.
		for(i=0; i<N_BYTES; i+=4){
			// Apply type 2 tables to all bytes, counting also shift rows selector.
			// One section ~ 1 column of state array, so select 1 column, first will
			// have indexes 0,4,8,12. Also take ShiftRows() into consideration.
			ires[i+0].l = edTab2[r][i+0][state.B[shiftOp[i/4+0*4]]].l;
			ires[i+1].l = edTab2[r][i+1][state.B[shiftOp[i/4+1*4]]].l;
			ires[i+2].l = edTab2[r][i+2][state.B[shiftOp[i/4+2*4]]].l;
			ires[i+3].l = edTab2[r][i+3][state.B[shiftOp[i/4+3*4]]].l;

			// XOR results of T2 boxes
			op8xor(ires[i+0], ires[i+1], edXTab[r][i/4][0], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], edXTab[r][i/4][1], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], edXTab[r][i/4][2], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage

			// Apply T3 boxes, valid XOR results are in ires[0], ires[4], ires[8], ires[12]
			// Start from the end, because in ires[i] is our XORing result.
			ires[i+3].l = edTab3[r][i+3][ires[i].B[3]].l;
			ires[i+2].l = edTab3[r][i+2][ires[i].B[2]].l;
			ires[i+1].l = edTab3[r][i+1][ires[i].B[1]].l;
			ires[i+0].l = edTab3[r][i+0][ires[i].B[0]].l;

			// Apply XORs again, now on T3 results
			// Copy results back to state
			op8xor(ires[i+0], ires[i+1], edXTab[r][i/4][3], ires[i+0]);  // 1 xor 2
			op8xor(ires[i+2], ires[i+3], edXTab[r][i/4][4], ires[i+2]);  // 3 xor 4
			op8xor(ires[i+0], ires[i+2], edXTab[r][i/4][5], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage
		}
		//
		// Copy results back to state
		// ires[i] now contains 32bit XOR result
		// We have to copy result to column...
		for(i=0; i<N_BYTES; i+=4){
			state.B[i/4+ 0] = ires[i].B[0];
			state.B[i/4+ 4] = ires[i].B[1];
			state.B[i/4+ 8] = ires[i].B[2];
			state.B[i/4+12] = ires[i].B[3];
		}
	}

#ifdef AES_BGE_ATTACK
	// If we are performing attack, we modified output bijection for 1 byte from 2 concatenated 4x4 bijections to one 8x8
	for(i=0; i<N_BYTES; i+=4){
		state.B[i/4+ 0] = edOutputBijection[r][i/4+ 0][state.B[i/4+ 0]];
		state.B[i/4+ 4] = edOutputBijection[r][i/4+ 4][state.B[i/4+ 4]];
		state.B[i/4+ 8] = edOutputBijection[r][i/4+ 8][state.B[i/4+ 8]];
		state.B[i/4+12] = edOutputBijection[r][i/4+12][state.B[i/4+12]];
	}
#endif
}

void BGEAttack::recoverPsi(Sset_t & S){
	Rmap_t R;
	int e = 1;

	// copy function map - writable copy (removable); fction map is hash -> idx
	fctionMap_t Stmp = S.fmap;

	// R = {'id', ['00']}; psi[id] = 0
	// R is fIdx -> beta
	R.insert(Rmap_elem_t(0,0));
	S.psi[0] = 0;

	// remove identity function from set
	Stmp.erase(S.fctions[0].hash);

	// start finding vector base & psi
	while(R.size() < GF256 && Stmp.size() > 0){
		// S <- S \ {f} // pick f from S and remove from S
		fctionMap_t::iterator it1 = Stmp.begin();
		BYTE fIdx = it1->second;
		Stmp.erase(it1);

		// if is already generated by some base, skip it
		if (R.count(fIdx)>0) continue;

		S.psi[fIdx] = e;
		//cout << "Taking fIdx="<<CHEX(fIdx)<<" as new base element e_i="<<CHEX(e)<<endl;

		Rmap_t Rcopy = R;
		for(Rmap_t::const_iterator it=Rcopy.begin(); it != Rcopy.end(); ++it){
			// Now compute elements (f \ocirc g, [e] ^ [n])
			// Functions should form vector space thus composition two of them should
			// result in another function from vector space. Thus compute hash of composition.

			// DEBUG
			//cout << " fIdx hash = " << S.fctions[fIdx].hash << endl;
			//cout << " firsthash = " << S.fctions[it->first].hash << endl;

			std::string nhash = composeFunction(S.fctions[fIdx].f, S.fctions[it->first].f);

			// hash should be contained in Stmp
			if (S.fmap.count(nhash)==0){
				cerr << "nhash["<<nhash<<"] was not found in VectorSpace; composition f \\ocirc g = " << CHEX(fIdx) << " o " << CHEX(it->first)
					 << "; betarepr: e="<<CHEX(e)<<"; ni="<<CHEX(it->second)<< endl;
				assert(S.fmap.count(nhash)>0);
			}

			BYTE nIdx = S.fmap[nhash];
			R.insert(Rmap_elem_t(nIdx, e ^ it->second));
			S.psi[nIdx] = e ^ it->second;
		}

		e *= 0x2;			// base move
	}

	if (R.size() < GF256 && Stmp.size() == 0){
		cerr << "Something bad happened, S is empty, R does not span whole vector space. size=" << R.size() << endl;
	}
}

int BGEAttack::deriveBset(Bset & bset, GenericAES & aes, bool permissive){
	struct coef_t {
		GF2E gfCoef;
		BYTE intCoef;
		GF2E gfCoefInv;
		BYTE intCoefInv;
		short int count;					  // number of occurences in one column in mixcol matrix
	};

	int dIdxs=0;							  // number of distincs indexes in one column in mixcol
	struct coef_t coefs[4];
	boost::unordered_map<BYTE, int> coefsMap; // simple map: index -> array of indexes

	// Determine coefficients used in MixCol in current AES used, compute inverses
	for(int i=0; i<4; i++){
		GF2E gfCoef = aes.mixColMat[i][0];
		BYTE curCoef = (BYTE) getLong(gfCoef);
		if(coefsMap.count(curCoef)==0){
			coefs[dIdxs].gfCoef = gfCoef;
			coefs[dIdxs].intCoef = curCoef;
			inv(coefs[dIdxs].gfCoefInv, gfCoef);
			coefs[dIdxs].intCoefInv = (BYTE) getLong(coefs[dIdxs].gfCoefInv);
			coefs[dIdxs].count=1;
			coefsMap.insert(std::pair<BYTE,int>(curCoef, dIdxs++));
		} else {
			coefs[coefsMap[curCoef]].count++;
		}
	}

	// just use MixCol matrix rows - more strict variant
	if (permissive==false){
		for(int i=0; i<4; i++){
			int a = coefsMap[(BYTE)getLong(aes.mixColMat[i][0])];
			int c = coefsMap[(BYTE)getLong(aes.mixColMat[i][1])];
			for(int j=0; j<4; j++){
				if (i==j) continue;
				int d = coefsMap[(BYTE)getLong(aes.mixColMat[j][0])];
				int b = coefsMap[(BYTE)getLong(aes.mixColMat[j][1])];

				GF2E beta = (coefs[a].gfCoef * coefs[b].gfCoef) * (coefs[c].gfCoefInv * coefs[d].gfCoefInv);
				BYTE intBeta = (BYTE) getLong(beta);
				if (doCout) cout << "i="<<i<<"; Beta=["<<CHEX(intBeta)
						<<"]; a="<<CHEX(coefs[a].intCoef)
						<<"; b="<<CHEX(coefs[b].intCoef)
						<<"; c="<<CHEX(coefs[c].intCoef) << ", inv="<<CHEX(coefs[c].intCoefInv)
						<<"; d="<<CHEX(coefs[d].intCoef) << ", inv="<<CHEX(coefs[d].intCoefInv)<<endl;

				if (bset.count(intBeta)==0) bset.insert(intBeta);
			}
		}
	} else {
		// Now generate combinations, 4*4*4*4 possibilities, I want to avoid 4 nested for loops, so with just one
		for(int i=0; i<GF256; i++){
			int a =  i       & 0x3; // index
			int b = (i >> 2) & 0x3;
			int c = (i >> 4) & 0x3;	// c is from same col as a, c is inverse
			int d = (i >> 6) & 0x3; // d is from same col as b, d is inverse
			if (a>=dIdxs || b>=dIdxs || c>=dIdxs || d>=dIdxs) continue; // in case there are some coefs with count>1
			if (a==c && coefs[a].count==1) continue;	// they are from same col, can have same value only if count is bigger than 1
			if (b==d && coefs[b].count==1) continue;    // they are from same col, can have same value only if count is bigger than 1
			if (a==d && coefs[a].count==1) continue;    // they are from same row, can have same value only if count is bigger than 1
			if (b==c && coefs[b].count==1) continue;    // they are from same row, can have same value only if count is bigger than 1
			if (a==d && b==c) continue; // not 2 rows same
			GF2E beta = (coefs[a].gfCoef * coefs[b].gfCoef) * (coefs[c].gfCoefInv * coefs[d].gfCoefInv);
			BYTE intBeta = (BYTE) getLong(beta);
			if (doCout) cout << "i="<<i<<"; Beta=["<<CHEX(intBeta)
					<<"]; a="<<CHEX(coefs[a].intCoef)
					<<"; b="<<CHEX(coefs[b].intCoef)
					<<"; c="<<CHEX(coefs[c].intCoef) << ", inv="<<CHEX(coefs[c].intCoefInv)
					<<"; d="<<CHEX(coefs[d].intCoef) << ", inv="<<CHEX(coefs[d].intCoefInv)<<endl;

			if (bset.count(intBeta)==0) bset.insert(intBeta);
		}
	}

	return 0;
}

NTL::GF2X BGEAttack::characteristicPolynomial(mat_GF2X_t & m){
	// This is recursive function, so if size is 2, return result directly
	if (m.n==1){
		return m.x[0][0];	// should never reach this point, but lets be defensive
	} else if (m.n==2){
		return (m.x[0][0] * m.x[1][1]) + (m.x[0][1] * m.x[1][0]);
	}

	NTL::GF2X poly = GF2XFromLong(0, 9);
	// Recursive expansion by rows
	mat_GF2X_t nm;
	nm.n = m.n-1;
	for(int k=0; k<m.n; k++){
		// optimization, multiplication by zero - save recursion steps
		if (m.x[0][k]==GF2X::zero()) continue;

		// create submatrix, we expand every time by first row
		for(int i=1; i<m.n; i++){
			for(int j=0,col=0; j<m.n; j++){
				if (j==k) continue;
				nm.x[i-1][col++] = m.x[i][j];
			}
		}

		poly = poly + m.x[0][k] * characteristicPolynomial(nm);
	}

	return poly;
}

NTL::GF2X BGEAttack::characteristicPolynomial(mat_GF2 & m){
	mat_GF2X_t subMatrix;
	int n = m.NumRows();
	subMatrix.n=n;

	// Transform to GF2X matrix with variable on main diagonal
	NTL::GF2X var = GF2XFromLong(2, 9);
	assert(m.NumCols()==n);	 		  // we can do this only for square matrix
	if (n==1) return (m[0][0] + var); // such a degenerate case

	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			NTL::GF2X cur = GF2XFromLong(m[i][j] == NTL::GF2::zero() ? 0:1, 9);
			subMatrix.x[i][j] = i==j ? (cur+var) : cur;
		}
	}

	return characteristicPolynomial(subMatrix);

}

int BGEAttack::proposition1(affineEquiv_t & ret, int r, int col, int syi, int syj, int sx){
	GF256_func_t yi, yj, yj_inv;        // pre-computed functions

	W128b state;
	for(int x=0; x<=0xff; x++){
		memset(&state, 0, sizeof(state));		// put 0 everywhere
		state.B[0+4*sx]=x;   state.B[1+4*sx]=x; // init with x values for y_0 in each column
		state.B[2+4*sx]=x;   state.B[3+4*sx]=x; // recall that state array is indexed by rows.

		this->Rbox(state, true, r, true);		// perform R box computation on input & output values
		yi[x]         = state.B[col + 4*syi];   //yi_inv[yi[x]] = x;
		yj[x]         = state.B[col + 4*syj];
		yj_inv[yj[x]] = x;
	}

	// We are looking for relation in a form:
	// yi(x, 00, 00, 00) ^ c = L(yj(x, 00, 00, 00))
	//
	//
	// Thus we are iterating by affine constant
	int c;
	for(c=0; c<=0xff; c++){
		// We are looking for linear relation between yi+c and yj. Thus just get values for L(e1), ..., L(e8)
		affineEquiv_t L;
		L.c = c;
		L.L.SetDims(8,8);         // mat_GF2 must be initialized to certain size before use
		L.Linv.SetDims(8,8);      // mat_GF2 must be initialized to certain size before use
		L.Lm[0] = 0;			  // must always hold for linear mapping

		// Problem is to find L:
		// yi(x, 00, 00, 00) ^ c = L(yj(x, 00, 00, 00))
		// It is easy to see that L is already defined if we have values for yi, yj, we just
		// have to:
		// 1. determine transformation for base vectors e1,...,e8
		// 2. test whether these transformation holds linearity property for whole space spanned by these base vectors
		//
		// Observe x is same for yi, yj, thus: a.) determine inverse for yj[x] = e_i; b.) L[ei] = yi[x] + c
		for(int i=0; i<8; i++){
			int lei = yi[yj_inv[1<<i]] ^ c;
			L.Lm[1<<i]   = lei;
			L.Lminv[lei] = 1<<i;

			// build transformation matrix
			colVectorT(lei, L.L, i);
		}

		// If transformation L is linear, then it has to have proper matrix inverse and determinant!=0
		GF2 determinant;
		NTL::inv(determinant, L.Linv, L.L);
		if (determinant==GF2::zero()){
			continue;	// matrix is not invertible -> not linear
		}

		// Now check whether already determined mapping is linear for whole space, just check.
		bool works = true;
		for(int x=0; x<=0xff; x++){
			BYTE ix = yj[x];
			BYTE iy = yi[x] ^ c;
			mat_GF2 mx = colVector(ix);    // mx for L * mx
			mat_GF2 my = colVector(iy);    // my for my =? L * mx
			if ((L.L * mx) != my){
				works=false;
				break;
			}

			L.Lm[ix] = iy;
			L.Lminv[iy] = ix;
		}

		// mapping is ok, return it
		if (works){
			ret = L;
			return 1;
		}
	}

	return 0; // not found
}

int BGEAttack::proposition2(mat_GF2 & L, baseVectors_t & out, mat_GF2 beta){
	// We have matrix L and matrix Beta, following equation holds from proposition 2:
	//
	// L * A~ = A~ * beta = R
	//
	// A~ is 8x8 matrix, has 64 elements. We take every such element as variable and
	// we are going to construct 64 equations with 64 unknowns and to solve this system
	// to obtain nontrivial solution, that will be returned as A~.
	//
	// Matrix multiplication gives us those equations. I will show this on first element of R matrix
	//
	// r_{00} = l_{00}*a_{00} + l_{01}*a_{10} + ... + l_{07}*a_{70} = a_{00}*b_{00} + a_{01}*b_{10} + ... + a_{07}*b_{70}
	// r_{ij} = \Sum_{k} l_{ik} * a_{kj} = \Sum_{k} a_{ik} * b_{kj}
	mat_GF2 eqSystem(INIT_SIZE, 64, 64);

	// Iterate over all 64 elements of R matrix and generate new equation in each iteration.
	int row=0;
	for(int i=0; i<8; i++){
		for(int j=0; j<8; j++, row++){
			for(int k=0; k<64; k++) eqSystem[row][k] = GF2::zero();

			// 1. do first summation (L*A~) and simultaneously second summation (A~ * beta)
			for(int k=0; k<8; k++){
				// Implemented by adding. Under some conditions both following lines will update
				// the same a variable a. If right side is 0, nothing happens (matrix is
				// initialized by default to 0.
				eqSystem[row][8*k + j] += L[i][k];
				eqSystem[row][8*i + k] += beta[k][j];
			}
		}
	}

	// Transform to row echelon form - we will need this if system is singular
	int rank = gauss(eqSystem);

	// Here we have system of 64 equations with 64 unknowns, we will solve it now.
	vec_GF2 b(INIT_SIZE, 64);
	vec_GF2 x(INIT_SIZE, 64);
	GF2 determinant;

	solve(determinant, x, eqSystem, b);
	if (determinant == GF2::zero()){
		int freeVariables = 64 - rank;

		// Return set of orthogonal vectors that solve this homogenous equation in case of singular system of
		// equations. Subspace generated by aforementioned vectors then solves equation.
		//
		// We have <freeVariables> free variables. Thus dimension of subspace that solves given system
		// has dimension = freeVaraibles.
		//
		// Return set of orthogonal vectors that span space that solves system. Each vector in this set
		// has only one free variable set to 1, others are set to zero. One vector has 64 elements.
		// If freeVariables=8 then we will return 8 vectors with 64 elements.
		//

		// Find orthogonal subspace that solves system
		for(int k=0; k<freeVariables; k++){
			vec_GF2 curV(INIT_SIZE, 64);  	// current vector

			// Exactly one free variable will be turned to 1 in base
			for(int i=0; i<freeVariables; i++){
				curV.put(63-i, i==k ? 1 : 0);
			}

			// Express other elements having assigned values by free variables
			for(int i=63-freeVariables; i>=0; i--){
				GF2 curVe = GF2::zero();
				// Express current element in current vector by lower vectors
				for(int j=i+1; j<64; j++){
					curVe += eqSystem[i][j] * curV[j];
				}

				curV.put(i, curVe);
			}

			// do self-test, this vector should solve system
			mat_GF2 outM(INIT_SIZE,8,8);
			for(int i=0; i<64; i++){
				outM.put(i/8, i%8, curV[i]);
			}

			// Final verification of given solution - is equation on the top correct now?
			mat_GF2 lhs = L * outM;
			mat_GF2 rhs = outM * beta;
			if(lhs != rhs){
				if (doCout) cout << "Something wrong with the result, equation does not hold in self-test..." << endl;
				if (doCout) cout << "Dimension of solution="<<freeVariables<<"; Faulty vector k="<<k<<"; Vector: " << endl;
				dumpMatrix(outM);
				return -4;
			}

			out.push_back(curV);
		}

		return 1;
	}

	// Vector X should contain resulting values for variables a_{00} .. a_{77}
	// We have only 1 solution here (can be trivial - homogenous system)
	bool trivialSolution=true;
	mat_GF2 outM(INIT_SIZE,8,8);
	for(int i=0; i<64; i++){
		outM.put(i/8, i%8, x[i]);
		if (x[i] != GF2::zero()) trivialSolution = false;
	}

	if (trivialSolution){
		return -1;
	}

	// Final verification of given solution - is equation on the top correct now?
	mat_GF2 lhs = L * outM;
	mat_GF2 rhs = outM * beta;
	if(lhs != rhs){
		if (doCout) cout << "Something wrong with the result, equation does not hold..." << endl;
		return -2;
	}

	out.push_back(x);
	return 0;
}

int BGEAttack::proposition3(prop3struct_t * out, GenericAES & aes, int r, int col, int row, baseVectors_t & Atild, long int vectorIdx){
	// Solving of this problem is based on iterating over 2 variables in GF(2^8) and checking
	// if resulting mapping is affine.
	//
	// Now prepare A~^{-1} from proposition 2 result. Just take first vector in Atild,
	// make matrix, take inverse, build lookup table
	if (Atild.size()==0) return -1;
	NTL::GF2 determinant;
	NTL::mat_GF2 AtildMat(INIT_SIZE, 8, 8);
	NTL::mat_GF2 AtildMatInv(INIT_SIZE, 8, 8);

	// Very defensive approach here, if Atild is composed of set of orthogonal vectors, try each from subspace
	int vectSpaceDim  = Atild.size();
	long int vIdxMax  = (long int) pow(2.0,vectSpaceDim);
	for(long int vIdx = 1; vIdx < vIdxMax ; vIdx++){

		// Reset matrix
		for(int x=0; x<64; x++){
			AtildMat.put(x/8, x%8, 0);
		}

		// Construct vector from subspace of solutions
		for(int d=0; d<vectSpaceDim; d++){
			if ((vIdx & (1<<d)) == 0) continue;		    // vIdx does not permit to use Atild[d] vector
			for(int x=0; x<64; x++){
				AtildMat.put(x/8, x%8, Atild[d].get(x));
			}
		}

		// Check if this matrix is OK
		NTL::inv(determinant, AtildMatInv, AtildMat);  // construct inverse - thats what we are looking for
		if (determinant==GF2::zero()){
			if (doCout) cout << "One Atild vector is not invertible but it should be: " << vIdx << endl;
			if (vectorIdx==vIdx) return -5;
		}

		if (determinant!=GF2::zero() && vectorIdx==vIdx){
			break;
		}
	}

	// No solution form invertible matrix
	if (determinant==GF2::zero()){
		if (doCout) cout << "There was problem with finding Atild inversion from given set of base vectors, dim=" << vectSpaceDim << endl;
		return -2;
	}

	// We determined Atild which to use, so continue with original proposition3
	return proposition3(out, aes, r, col, row, AtildMatInv);
}

int BGEAttack::proposition3(prop3struct_t * out, GenericAES & aes, int r, int col, int row, mat_GF2 AtildMatInv){
	// Build lookup table for Atild
	GF256_func_t AtildInvLoT;
	for(int x=0; x<=0xff; x++){
		mat_GF2 xMat   = colVector(x);
		mat_GF2 res    = AtildMatInv * xMat;
		GF2X    resGF  = colVector_GF2X(res, 0);
		AtildInvLoT[x] = getLong(resGF);
	}

	//
	// Construct those mappings y_row. This can be similar to prop1, but
	// here we are iterating over 2 variables and CHECKING if resulting mapping is affine.
	//
	int PiOK=0;
	for(int p=0; p<4; p++){
		GF256_func_t y_row;                         // pre-computed functions
		W128b state;
		for(int x=0; x<=0xff; x++){
			memset(&state, 0, sizeof(state));		// put 0 everywhere
			state.B[col + 4*p] = x;                 // In each iteration P_i, we set X to different position, namely x0,x1,x2,x3

			this->Rbox(state, true, r, true);		// perform R box computation on input & output values
			y_row[x] = state.B[col + 4*row];
		}

		// It is important to set dimensions to mat_GF2 before starting working with it, otherwise --> segfault
		out->P[p].L.SetDims(8,8);
		out->P[p].Linv.SetDims(8,8);
		out->P[p].valid=false;

		// Now iterate over 2 variables in GF(256), construct mapping and check affinity
		bool foundAffine=false;
		for(int delta=1; delta<=0xff && foundAffine==false; delta++){
			// multiplication matrix
			mat_GF2 deltMultM = aes.makeMultAMatrix(delta);
			// multiplication lookup table
			GF256_func_t deltMultLoT;
			for(int x=0; x<=0xff; x++){
				mat_GF2 xMat   = colVector(x);
				mat_GF2 res    = deltMultM * xMat;
				GF2X    resGF  = colVector_GF2X(res, 0);
				deltMultLoT[x] = getLong(resGF);
			}

			// Iterate over c constant \in GF(2^8)
			for(int c=0; c<=0xff && foundAffine==false; c++){
				out->P[p].delta = delta;
				out->P[p].c     = c;

				// Constructing Pi mapping
				// P0(x) = (SboxInv * deltMult_0 * AtildInv)(y_row( x, 00, 00, 00) + c0)
				// P1(x) = (SboxInv * deltMult_1 * AtildInv)(y_row(00,  x, 00, 00) + c1)
				// P2(x) = (SboxInv * deltMult_2 * AtildInv)(y_row(00, 00,  x, 00) + c2)
				// P3(x) = (SboxInv * deltMult_3 * AtildInv)(y_row(00, 00, 00,  x) + c3)
				GF256_func_t linPart; // Linear part of above equations - to test for linearity later
				for(int x=0; x<=0xff; x++){
					out->P[p].affineMap[x] = aes.sboxAffineInv[deltMultLoT[AtildInvLoT[y_row[x] ^ c]]];
					linPart[x]             = x==0 ? 0 : out->P[p].affineMap[x] ^ out->P[p].affineMap[0];
				}

				// Now check, if constructed mapping is affine. Store affine constant for later = Af(0)
				out->P[p].affineConst = out->P[p].affineMap[0];

				// Check, if given linear transformation works
				// in the same way on each element in space, using transformation results on canon. base
				bool isTestOK=true;
				for(int x=0; x<=0xff; x++){
					int shouldBe = linPart[x];
					int realResult = 0;
					for(int d=0; d<8; d++){
						realResult ^= ((x & (1<<d)) > 0 ? linPart[1<<d] : 0);
					}

					if (shouldBe != realResult){
						isTestOK=false;
						break;
					}
				}

				if (isTestOK==false){
					continue;
				}

				// Verify linearity of linPart. Easy procedure: 1. construct matrix representation by
				// applying transformation on canon. base,      2. test invertibility
				for(int x=0; x<64; x++){
					out->P[p].L.put(x/8, x%8, 0);
					out->P[p].Linv.put(x/8, x%8, 0);
				}

				for(int x=0; x<8; x++){
					colVectorT(linPart[1<<x], out->P[p].L, x);
				}

				NTL::GF2 determinant;
				NTL::inv(determinant, out->P[p].Linv, out->P[p].L);
				if (determinant == GF2::zero()){
					if (doCout) cout << "Strange situation, determinant=0, but transformation worked on space... p= "<<p<<";Delta="<<delta<<"; c="<<c<<endl;
					continue;
				}

				out->P[p].valid = true;
				foundAffine = true;
				break;
			}
		}

		PiOK |= foundAffine ? 1<<p : 0;
	}

	// Finally, compute c4 = y_row(00,00,00,00) for q_row = c0 ^ c1 ^ c2 ^ c3 ^ c4
	W128b state;
	memset(&state, 0, sizeof(state));  // put 0 everywhere
	this->Rbox(state, true, r, true);  // perform R box computation on input & output values
	out->c4 = state.B[4*row + col];    // In order to make sense, compute y_row(00,00,00,00)

	// Additional stuff - derive gamma and mix column coefficients
	// delta^{-1}_i = gamma * alfa_{0,i}
	// There are some values twice -> 01 (because 01 is also twice in each row in MC).
	if(PiOK==0xf){
		int idxTwice=-1;              // index that has 01 as MC coefficient
		for(int d1=0;d1<4;d1++){
			for(int d2=0;d2<4;d2++){
				if (d1!=d2 && out->P[d1].delta == out->P[d2].delta){ idxTwice = d1; break; }
			}
		}

		if (idxTwice==-1){
			if (doCout) cout << "Weird thing happened, there should be 2 deltas with same value, is MC as it should be?" << endl;
			return PiOK;
		}

		// Now we are able to derive gamma, delta^{-1}_{idxTwice} = gamma * 01 = gamma
		// alfa = delta^{-1} * gamma^{-1}
		GF2E gammaInvGF = GF2EFromLong(out->P[idxTwice].delta, 8);
		GF2E gammaGF    = NTL::inv(gammaInvGF);
		out->gamma      = (BYTE) getLong(gammaGF);
		for(int d=0; d<4; d++){
			if (out->P[d].valid==false) continue;
			GF2E deltaCurGF    = GF2EFromLong(out->P[d].delta, 8);
			GF2E deltaCurInvGF = NTL::inv(deltaCurGF);
			GF2E alfaGF        = gammaInvGF * deltaCurInvGF;
			out->P[d].alfa_0   = (BYTE) getLong(alfaGF);
		}
	}

	return PiOK;
}

int BGEAttack::recoverQj(GenericAES & aes, int r, int col, int row, const mat_GF2 A0, BYTE alfa00, BYTE alfa_row_0, mat_GF2 & Aj, BYTE * qj){
	assert(r>=0 && r<N_ROUNDS);
	assert(col>=0 && col<=3);
	assert(row>=1 && row<=3);

	// From the Proposition 1, if we know already linear part A0 of Q0, we can this way determine
	// Aj linear part of Qj for other mappings
	//
	// Proposition 1 solves:
	// y_i(x0,00,00,00) = L(y_j(x0,00,00,00)) + c    while L and c is returned
	// 	 L   = A_i * (\alfa_{i,0} * \alfa_{j,0}^{-1}) * A_j^{-1}  while we know L, A_i, alfa...
	//   A_j = L^{-1} * A_i * (\alfa_{i,0} * \alfa_{j,0}^{-1})
	affineEquiv_t L;
	int prop1result = proposition1(L, r, col, 0, row, 0);
	if (prop1result==0) {
		if (doCout) cout << "recoverAj fail: Proposition1 failed for r="<<r<<"; col="<<col<<"; (y0, y"<<row<<") " << endl;
		return -1;
	}

	// Here we need inverse of alfa_row_0
	GF2E alfa00GF      = GF2EFromLong(alfa00, 8);
	GF2E alfaRowGF     = GF2EFromLong(alfa_row_0, 8);
	GF2E alfaRowInvGF  = NTL::inv(alfaRowGF);
	GF2E alfaGF        = alfa00GF * alfaRowInvGF;
	BYTE alfa          = (BYTE) getLong(alfaGF);
	mat_GF2 alfaMatrix = aes.makeMultAMatrix(alfa);

	// compute linear part from L
	Aj = L.Linv * A0 * alfaMatrix;

	// Part 2 - constant part extraction
	// If we put known Aj linear part to proposition 3 we don't have to bother with A~_j and gamma.
	// Then delta_i^{-1} = afa_j,i
	mat_GF2 AjInv;
	GF2 determinant;
	NTL::inv(determinant, AjInv, Aj);
	if (determinant == GF2::zero()){
		if (doCout) cout << "Cannot compute inverse - should never happen!!" << endl;
		return -2;
	}

	prop3struct_t * prop3 = new prop3struct_t;
	int prop3result = proposition3(prop3, aes, r, col, row, AjInv);
	if (prop3result != 0xf) {
		delete prop3;
		if (doCout) cout << "some problem with proposition 3 in recover Qj occurred. result="<<prop3result<<endl;
		return -3;
	}

	// We don't really care about delta, important is, that Aj should be correct, so we
	// can directly recover qj - see end of section 3.3
	(*qj) =   prop3->c4
			^ prop3->P[0].c
			^ prop3->P[1].c
			^ prop3->P[2].c
			^ prop3->P[3].c;

	// Here we can do another test to verify that we are doing well - just visually, not checking against MC
	// delta_i^{-1} = afa_j,i
	for(int x=0; x<4; x++){
		GF2E deltaGF    = GF2EFromLong(prop3->P[x].delta, 8);
		GF2E deltaInvGF = NTL::inv(deltaGF);
		BYTE deltaInv   = getLong(deltaInvGF);
		if (doCout) cout
			 << "   recoverQj self-test; r="<<r<<"; col="<<col<<"; (y0, y"<<row<<"); P["<<x
			 <<"].deltaInv="<<CHEX(deltaInv)
			 <<"; alfa_{"<<row<<","<<x<<"}="<< CHEX(prop3->P[x].alfa_0) << endl;
	}

	if (doCout) cout << "   recoverQj; q = " << CHEX(*qj) << "; gamma=" << CHEX(prop3->gamma) << "; " << endl;

	delete prop3;
	return 0;
}

int BGEAttack::recoverCipherKey(GenericAES & aes, BYTE roundKeys[2][16], vec_GF2E& encKey){
	encKey.SetLength(16);

	// At first determine correct round and Rcon constant
	// Verify that given round key is OK
	// w9  = w8  + w5
	// w10 = w9  + w6
	// w11 = w10 + w7
	for(int col=1; col<4; col++){
		for(int row=1; row<4; row++){
			if (roundKeys[1][4*row+col] != (roundKeys[1][4*row+col-1] ^ roundKeys[0][4*row+col])){
				if (doCout) cout << "Error in round key verification, you passed invalid round keys; col="<<col<<"; row="<<row << endl;
				return -1;
			}
		}
	}

	//
	// Keys are probably OK, now recover RCON constant and round to which key belongs
	//
	vec_GF2E gw7(INIT_SIZE, 4);
	vec_GF2E  w7(INIT_SIZE, 4);
	for(int i=0; i<4; i++)  w7[i] = GF2EFromLong(roundKeys[0][4*i+3], 8);
	for(int i=0; i<4; i++) gw7[i] = GF2EFromLong(roundKeys[1][4*i] ^ roundKeys[0][4*i], 8); // w8 = g(w7) + w4

	int rconIdx = -1;
	for(int rc=0; rc<16; rc++){
		vec_GF2E gw7_copy(INIT_SIZE, 4);
		for(int i=0; i<4; i++) gw7_copy[i] = gw7[i];

		// Revert RC[rc] application on first element
		gw7_copy[0] = GF2EFromLong(getLong(gw7[0] + aes.RC[rc]),8);
		// Inverse Sbox application
		for(int i=0; i<4; i++) gw7_copy[i] = GF2EFromLong(aes.sboxAffineInv[getLong(gw7_copy[i])], 8);
		// Inverse left shift -> shift to the right, respectively test it with taking shift into account
		bool correctRcon=true;
		for(int i=0; i<4; i++){
			if (w7[(i+1) % 4] != gw7_copy[i]) {
				correctRcon=false;
				break;
			}
		}

		if(correctRcon){
			rconIdx = rc;
			if (doCout) cout << "We have correct Rcon! rconIdx="<<rc<<endl;
			break;
		}
	}

	if (rconIdx==-1) return -3;

	// now perform reverse key schedule to get cipher key
	mat_GF2E w(INIT_SIZE, 4, 4);
	for(int i=0; i<16; i++) w[i/4][i%4] = GF2EFromLong(roundKeys[0][i], 8);
	for(int rc=rconIdx-1; rc>=0; rc--){
		mat_GF2E wp(INIT_SIZE, 4, 4); // previous

		// we can simply derive w1, w2, w3
		for(int col=1; col<4; col++){
			for(int row=0; row<4; row++){
				wp[row][col] = GF2EFromLong(getLong(w[row][col] + w[row][col-1]), 8);
			}
		}

		// w0 = w4 + g(w3), g(w3) we already have, so compute g(w3)
		// 1. copy w3 to w0
		// 2. apply g to w0
		// 3. add w4 to w0
		for(int i=0; i<4; i++){
			wp[i][0] = GF2EFromLong(aes.sboxAffine[getLong(wp[(i+1) % 4][3])], 8); // take shift into account + apply sbox
		}

		// Apply RCON
		wp[0][0] = GF2EFromLong(getLong(wp[0][0] + aes.RC[rc]), 8);

		// Add w4 and we have complete w0
		for(int i=0; i<4; i++) wp[i][0] = GF2EFromLong(getLong(wp[i][0] + w[i][0]), 8);

		// copy back
		for(int i=0; i<16; i++) w[i/4][i%4] = GF2EFromLong(getLong(wp[i/4][i%4]), 8);

		if (doCout) cout << "RC=" << rc << "; previousKey: "<<endl;
		dumpMatrix(w);
	}

	for(int i=0; i<16; i++) encKey[i] = GF2EFromLong(getLong(w[i%4][i/4]), 8);
	return 0;
}

int BGEAttack::run(BYTE * key, keySize keyLen) {
	GenericAES defAES;
	defAES.init(0x11B, 0x03);

#ifndef AES_BGE_ATTACK
	cerr << "Cannot proceed with attack if \"AES_BGE_ATTACK\" is not defined, we are missing required additions"<<endl;
	exit(1);
#endif

	WBAESGenerator generator;
	if (doCout) cout << "Generating AES..." << endl;

	generator.useDualAESARelationsIdentity=false;
	generator.useDualAESIdentity=false;
	generator.useDualAESSimpeAlternate=false;
	generator.useIO04x04Identity=false;
	generator.useIO08x08Identity=false;
	generator.useMB08x08Identity=false;
	generator.useMB32x32Identity=false;


	this->wbaes  = new WBAES;
	this->coding = new ExtEncoding;
	generator.generateExtEncoding(this->coding, WBAESGEN_EXTGEN_ID);
	generator.generateTables(key ? key : GenericAES::testVect128_key, keyLen,
							 this->wbaes, this->coding, true);  if (doCout) cout << "AES ENC generated" << endl;
	generator.generateTables(key ? key : GenericAES::testVect128_key, keyLen,
							 this->wbaes, this->coding, false); if (doCout) cout << "AES DEC generated" << endl;

	// WBAES changed to state with affine matching bijections at round boundaries.
	if (doCout) cout << "Going to test WBAES before modifying tables" << endl;
	generator.testComputedVectors(true, this->wbaes, this->coding);

	if (doCout) cout << "Starting the BGE attack" << endl;
	int toReturn = this->attack();

	delete this->wbaes;
	delete this->coding;

	this->wbaes  = NULL;
	this->coding = NULL;

	return toReturn;
}

int BGEAttack::attack(void) {
	// The aim of this attack is to extract round keys, so generate them here to be able to compare results

	GenericAES defAES;
	defAES.init(0x11B, 0x03);

	bool encrypt = true;
	WBAESGenerator generator;
	W128b state;

	assert(this->wbaes!=NULL);

	vec_GF2E defaultKey;						// key for default AES in GF2E representation
	vec_GF2E expandedKey;						// expanded key for default AES
	const int * nextTbox = shiftT2;             // attack is not yet implemented for decryption
	generator.BYTEArr_to_vec_GF2E(GenericAES::testVect128_key, KEY_SIZE_16, defaultKey);			// convert BYTE key to GF2E key
	defAES.expandKey(expandedKey, defaultKey, KEY_SIZE_16);	// key schedule for default AES
	if (doCout) cout << "Expanded key: " << endl;
	if (doCout) dumpVector(expandedKey);


	//
	//
	// Attack below...
	//
	//

	// Recover affine parts of output mapping Q for round 0..9
	int r = 0;
	int row=0;
	int x = 0;
	int i = 0;
	int j = 0;
	int c1;

	// At first compute base function f_{00}, we will need it for computing all next functions,
	// to be exact, its inverse, f^{-1}_{00}
	if (doCout) cout << "Allocating memory for the attack" << endl;
	SsetPerRound_t * Sr = new SsetPerRound_t[10];
	Qaffine_t * Qaffine = new Qaffine_t;
	if (doCout) cout << "Memory allocated; Sr=" << dec << (sizeof(SsetPerRound_t)*10) << "; Qaffine=" << dec << sizeof(Qaffine_t) << endl;
	if (doCout) cout << "Memory allocated totally: " << (((sizeof(SsetPerRound_t)*10) + sizeof(Qaffine_t)) / 1024.0 / 1024.0) << " MB" << endl;

	if (doCout) cout << "Starting attack phase 1 ..." << endl;
	for(r=0; r<9; r++){
		// Init f_00 function in Sr
		for(i=0; i<AES_BYTES; i++){
			Sr[r].S[i%4][i/4].f_00.c1 = 0;
		}

		//
		// Compute f(x,0,0,0) function for each Q_{i,j}
		//
		// x x x x        y_{0,0} y_{1,0} ..
		// 0 0 0 0   R    y_{0,1} y_{1,1} ..
		// 0 0 0 0  --->  y_{0,2} y_{1,2} ..
		// 0 0 0 0        y_{0,3} y_{1,3} ..
		//
		if (doCout) cout << "Generating f_00 for round r="<<r<<endl;
		for(x=0; x<=0xff; x++){
			memset(&state, 0, sizeof(state));		// put 0 everywhere
			state.B[0]=x;   state.B[1]=x; 			// init with x values for y_0 in each column
			state.B[2]=x;   state.B[3]=x;           // recall that state array is indexed by rows.

			this->Rbox(state, true, r, true);		// perform R box computation on input & output values
			for(i=0; i<AES_BYTES; i++){
				fction_t & f00 = Sr[r].S[i%4][i/4].f_00;
				f00.f[x] = state.B[i];
				f00.finv[state.B[i]] = x;
			}
		}


		// f(x,0,0,0) finalization - compute hash of f00 function
		for(i=0; i<AES_BYTES; i++){
			Sr[r].S[i%4][i/4].f_00.initHash();
		}

		// now generate f(x,0,0,0) .. f(x,0xff,0,0) functions, generate then corresponding sets and whole Sr for round r
		if (doCout) cout << "Generating set S..." << endl;
		for(c1=0; c1<=0xff; c1++){
			// init f_c1 functions
			for(i=0; i<AES_BYTES; i++){
				Sr[r].S[i%4][i/4].fctions[c1].c1 = c1;
			}

			// Now generating function f(x,c1,0,0)
			//
			//  x  x  x  x        y_{0,0} y_{1,0} ..
			// c1 c1 c1 c1   R    y_{0,1} y_{1,1} ..
			//  0  0  0  0  --->  y_{0,2} y_{1,2} ..
			//  0  0  0  0        y_{0,3} y_{1,3} ..
			//
			//
			// Generating S[col][row] fctions as f(f^{-1}(x,0,0,0), c1, 0, 0)
			// To be exact it is f( f^{-1}(x,0,0,0)[y0], c1, 0, 0)[y0]   for y0
			//
			// So y0 must match in both functions (inverse f and f), thus we have to iterate y_i, i \in [0,3]
			// One calculation for y_0 (4 columns simultaneously), another iteration for y_1 etc...
			// This loop iterates over ROWS of functions f.
			for(row=0; row<4; row++){
				for(x=0; x<=0xff; x++){
					// Init input arguments to (x,c1,0,0) for each column
					memset(&state, 0, sizeof(state));

					state.B[0]=Sr[r].S[0][row].f_00.finv[x]; // f^{-1}(x,0,0,0)[y0], 1.st column
					state.B[1]=Sr[r].S[1][row].f_00.finv[x]; // f^{-1}(x,0,0,0)[y1], 2.nd column
					state.B[2]=Sr[r].S[2][row].f_00.finv[x]; // f^{-1}(x,0,0,0)[y2], 3.rd column
					state.B[3]=Sr[r].S[3][row].f_00.finv[x]; // f^{-1}(x,0,0,0)[y3], 4.th column

					// init second argument as c1
					state.B[4]=c1;  state.B[5]=c1;
					state.B[6]=c1;  state.B[7]=c1;

					this->Rbox(state, true, r, true);		// perform R box computation on input & output values
					for(i=0; i<4; i++){  					// if we are in row=3, take result from 3rd row of state array: 8,9,10,11. Iterating over columns
						Sr[r].S[i][row].fctions[c1].f[x] = state.B[4*row + i];
						Sr[r].S[i][row].fctions[c1].finv[state.B[4*row + i]] = x;
					}
				}
			}
		}

		// Insert hash of functions to corresponding S set
		for(i=0; i<AES_BYTES; i++){
			Sr[r].S[i%4][i/4].fmap.clear();
			for(j=0; j<GF256; j++){
				Sr[r].S[i%4][i/4].fctions[j].initHash();
				std::string & hash = Sr[r].S[i%4][i/4].fctions[j].hash;

				// By iterating c1 over GF256 we constructed vector space, with 256 elements, thus each function
				// should be different. We are checking, whether this assumption is OK and vector space is correct.
				if (Sr[r].S[i%4][i/4].fmap.count(hash)!=0){
					cerr << "Sr["<<r<<"].["<<(i%4)<<"]["<<(i/4)<<"] set already contains hash[" << hash << "]; for j=" << j << std::endl;
					assert(Sr[r].S[i%4][i/4].fmap.count(hash)!=0);
				}

				Sr[r].S[i%4][i/4].fmap.insert(fctionMap_elem_t(hash, j));
			}
		}


		//
		// Now we have Sr[r] generated (set S for whole round r)
		// According to the paper, we now have to construct isomorphism
		// \phi S --> GF(2)^8
		//      Q \op \oplus_{\beta} \op Q^{-1}  --> [\beta]
		//
		// But we don't know this isomorphism. Instead we will construct isomorphism \psi
		// \psi = L^{-1} \op \phi, where L is base change matrix (linear transformation) [e_i] --> [\Beta_i]
		//
		// We can find base set in S and express every element in S by means of base. This way we will find \psi.
		//
		for(i=0; i<AES_BYTES; i++){
			if (doCout) cout << "Recovering psi for Sidx="<<i<<endl;
			Sset_t & curS = Sr[r].S[i%4][i/4];
			recoverPsi(curS);

			// derive Q~
			//                   +-------------------------------- Q^{-1} \ocirc L^{-1}; 1,2 part of commutative graph
			//                   |                          +----- 3.rd part of commutative graph (on purpose, to be able to construct PSI)
			//                   |                          |
			// 1. Q~('00') = L^{-1}(Q^{-1}('00')) \oplus (Q \ocirc L)^{-1}('00') = ['00']
			// 2. f = Q~ \ocirc \oplus_{PSI(f)} \ocirc Q~{-1}
			//
			// (1., 2.) ==>f('00') = Q~(PSI(f)) ==>
			// Q~: x is PSI(f) for some f (we have 256 of them)
			// Q~: y is f('00')  --- || ---
			boost::unordered_set<BYTE> keySet;
			boost::unordered_set<BYTE> valSet;
			for(j=0; j<GF256; j++){
				    x = curS.psi[j];
				int y = curS.fctions[j].f[0];

				// correctness, verification that it spans completely
				if (keySet.count(x)>0){
					cerr << "Key x=["<<x<<"] already defined in Q["<<r<<"]["<<i<<"].f["<<x<<"]"<< endl;
				}
				if (valSet.count(y)>0){
					cerr << "Val y=["<<y<<"] already defined in Q["<<r<<"]["<<i<<"].f["<<y<<"]"<< endl;
				}

				Qaffine->Q[r][i].f[x]    = y;
				Qaffine->Q[r][i].finv[y] = x;
				keySet.insert(x);
				valSet.insert(y);
				//cout << "Q["<<r<<"]["<<i<<"].f["<<x<<"] = " << y << endl;
			}

			Qaffine->Q[r][i].initHash();
			if (doCout) cout << "Q~ recovered; hash=" << Qaffine->Q[r][i].hash << endl;
			if (doCout) cout << "PSI function: " << hashFunction(curS.psi) << endl;
		}

		if (doCout) cout << "PSI recovered for all sets in given round" << endl;
	}

	// Now we can reduce P, Q are non-linear & matching state to P,Q are affine and matching
	// P^{r+1}_{i,j}  = Q^r_{i,j} (they are matching)
	// P~^{r+1}_{i,j} = Q~^r_{i,j} (they are matching, also holds for new affine mappings)
	//
	// Reduction:
	// Q~^{-1}(Q(x)) = L^{-1}(x) \oplus L^{-1}(Q^{-1}('00'))  .... what is affine (last term is constant)
	//
	// Thus it is enough to apply Q~^{-1} to the end of R-box,
	// To preserve matching relation, we also apply Q~ before T2 boxes,
	// Before it was: MixCol -> Q           -> |new round|       -> Q^{-1} -> T2box
	// Now it will be MixCol -> Q -> Q~{-1} -> |new round| Q~    -> Q^{-1} -> T2box
	// It is easy to see that it matches...
	//
	// Thus P^{r+1} conversion to affine matching: P^{r+1}(Q~(x)) =
	// = Q^{-1}(Q~(x)) = Q^{-1}( Q(L( x   \oplus L^{-1}(Q^{-1}('00')))))
	//                 =           L( x   \oplus L^{-1}(Q^{-1}('00')))
	//                 =           L( x ) \oplus        Q^{-1}('00')      ==> affine

	// Applying above description and reduce to affine matching transformations
	// This would also work up to N_ROUNDS-1, but we would have to change input encoding of T1
	// in next roud - code overhead. Not needed...

	// Q will be on transfered to AES output bijection layer - additional one
	// P will be on this->wbaes.eTab2[r][?][0..256]
	if (doCout) cout << "Going to transform AES input/output matching bijection to affine" << endl;
	for(r=0; r<N_ROUNDS-2; r++){
		for(i=0; i<AES_BYTES; i++){
			// Choice of indexes for eTab2 explained below:
			// Since we are interested in output encoding, we don't really care about input to Rbox (shiftrows is ignored).
			// We have Q[r][i \ in [1..16]] which corresponds to state array indexed by rows.
			// In WBAES implementation, shift rows will follow, thus first 4 T2boxes will be feeded by indexes of state array:
			// state[0,7,10,13]; next 4 T2 boxes will be feeded by state[4,1,14,11]
			// The point is that nextTbox \ocirc shiftRows (during evaluation) = identity
			AES_TB_TYPE2 & curT2     = this->wbaes->eTab2[(r+1)][nextTbox[i]];    // T2 table in next round connected to particular HILO(xtb1, xtb2)
			XTB & curXtb1            = this->wbaes->eXTab[r][i%4][5][2*(i/4)+0];  // 3.rd index - XOR table number 5, last in round
			XTB & curXtb2            = this->wbaes->eXTab[r][i%4][5][2*(i/4)+1];
#ifdef AES_BGE_ATTACK
			GF256_func_t & outBiject = this->wbaes->eOutputBijection[r][i];
#endif

			// Copy original values of tables to be consistent in changes, not to change value that will be
			// needed in future in original form (before application of transformation)
			AES_TB_TYPE2 tmpT2;
			memcpy(tmpT2, curT2, sizeof(AES_TB_TYPE2));

			// Debug bijection for testing networking
			//BIJECT8X8 tmpF;
			//BIJECT8X8 tmpFinv;
			//generator.generate8X8Bijection(&tmpF, &tmpFinv, false);

			// which bijection to use? For debug we can use random one, or computed Q
			GF256_func_t & bijectionF    = Qaffine->Q[r][i].f;    //tmpFinv;
			GF256_func_t & bijectionFinv = Qaffine->Q[r][i].finv; //tmpF;

			// Set of really paranoid asserts, before any change in WBAES
			for(x=0; x<GF256; x++){
				assert(tmpT2[x].l==curT2[x].l);
				assert(bijectionFinv[bijectionF[x]] == x);
				assert(bijectionF[bijectionFinv[x]] == x);
			}

			//
			// Remember already seen values from XOR Boxes.
			//
			// Changing XOR tables is a bit tricky in implementation, because they are 8->4 bit functions.
			// Thus if we want to apply 8bit transformation on output on 2 XOR boxes, we have to concatenate
			// result of 2 XOR boxes together and apply transformation, but...
			//   1. For sure we have x!=y s.t.: XTB[X] == XTB[Y], thus XTB is not injective (cannot be) - probem with inverse in T2.
			//   2. App. transf. on (XTB1[x] || XTB2[x]) is wrong, since iterating by x does not have to exhaust whole 2^8 values
			// Here may be cases where x!=y : (XTB1[x] || XTB2[y]) \notin {(XTB1[x] || XTB2[x]); x  \in 0..255}
			// this will cause problem after applying bijection, not all values are covered, thus if we in T2 box apply inverse
			// it is wrong. For this reason we need following set to remember already mapped values and iterate in 2 nested for loops.
			//
			// Main reason for this is: HILO(xtb1, xtb2) is of type 16->8
			//
			// We also have to realize that if we want to add 8x8 bijection at the output of AES round, which is constructed
			// as concatenation of 4x4 bijections, we cannot do it with just modifying XOR tables since they cant express encoding
			// 1 byte with 8x8 bijection - look on WBAES.h for more details & further explanation.
			//
			int y=0;
			for(y=0; y<GF256; y++){
#ifdef AES_BGE_ATTACK
				for(x=0; x<GF256; x++){
					BYTE newXtbX = HILO(curXtb1[x], curXtb2[y]);
					outBiject[newXtbX] = bijectionFinv[newXtbX];	// this will be repeated 256 times, but with same value, so OK
				}
#endif

				// For updating T2 boxes we only need 1 for loop iterating over field since function is 8->8
				curT2[y].l = tmpT2[bijectionF[y]].l;
			}
		}
	}

	// WBAES changed to state with affine matching bijections at round boundaries.
	if (doCout) cout << "Going to test WBAES after modifying tables" << endl;
	if (this->coding != NULL) {
		generator.testComputedVectors(true, this->wbaes, this->coding);
	} else {
		ExtEncoding tmpCoding;
		generator.generateExtEncoding(&tmpCoding, WBAESGEN_EXTGEN_ID);
		generator.testComputedVectors(true, this->wbaes, &tmpCoding);
	}

	//
	// Attack, proceeding to phase 2 - removing affine parts for output encoding
	//

	// Derive B set
	// According to the paper, section 3.3
	// B = { \Alfa_{0,0}\Alfa_{1,1}\Alfa{0,1}^{-1}\Alfa_{1,0}^{-1} | \Alfa_{i,j} \ in {01, 02, 03}}
	//
	// We will need this later, to determine which element from B set forms particular linear transformation
	Bset bset2;
	deriveBset(bset2, defAES, true);

	if (doCout) {
		cout << "Bset: { ";
		for (Bset::const_iterator it = bset2.begin(); it != bset2.end(); ++it) {
			cout << CHEX(*it) << ", ";
		}
		cout << "} " << endl;
	}

	// Now we can generate multiplication matrices from given constants \in B, then compute
	// characteristic polynomials for them.
	//
	// L = A0 * beta * A0^{-1} ~ beta ==> poly(L) = poly(beta)
	// From this we can determine possible beta values for given L as matrix.
	//
	charPolynomialMultimap_t charPolyMap;
	for(Bset::const_iterator it = bset2.begin(); it!=bset2.end(); ++it){
		mat_GF2 m = defAES.makeMultAMatrix(*it);
		if (doCout) cout << "Computing characteristic polynomial for: " << CHEX(*it) << "; ....";

		GF2X poly = characteristicPolynomial(m);
		if (doCout) cout << "Characteristic polynomial=" << poly << "; hex=" << GF2EHEX(poly) << endl;

		// Insert characteristic polynomial to multi-map. Later it will be used to determine B
		charPolyMap.insert(charPolynomialMultimap_elem_t(getLong(poly), *it));
	}

	//
	// Out attack now recovers round keys in rounds 2,3,4,5.
	//
	const int roundStart  = 2;
	const int rounds2hack = 4;
	const int vectorIdx   = 1;		 // which vector to use if proposition 2 returned subspace
	BYTE roundKeys[rounds2hack][16];

	// Proposition 3 data structure for 3 rounds and 4 columns in each round
	std::vector<prop3struct_t*> prop3structVector;
	prop3structVector.reserve(rounds2hack*4);
	for(int x=0; x<rounds2hack*4; x++){
		prop3structVector[x] = new prop3struct_t;
	}

	for(r=roundStart; r<(roundStart+rounds2hack); r++){
		// Our attack now determines Q0 for given round and col totaly, with round key
		// extraction for k0,k1,k2,k3. We want to determine round key for whole round
		// so iterate over columns
		for(int col=0; col<4; col++){
			if (doCout) cout << endl << "Extracting round keys for r="<<r<<"; col="<<col<<endl;
			// We are using section 3.2, proposition 1 to compute L0, L1 linear parts of affine
			// transformation between (yi, yj).
			// Prop1: yi(x0, 00, 00, 00) = L(yj(x0, 00, 00, 00)) + c
			//
			// Then we compute L0 using prop1 on (y_0, y_1) when x0 is varying
			// In the same way L1 using prop1 on (y_0, y_1) when x1 is varying
			// Then we get L = L0 * L1^{-1} = A0 * beta * A0^{-1}     where A0 is linear part of Q0 output mapping.
			//
			affineEquiv_t L0, L1;
			int resL0 = proposition1(L0, r, col, 0, 1, 0);
			int resL1 = proposition1(L1, r, col, 0, 1, 1);
			if (resL0==0 || resL1==0){
				if (doCout) cout << "One of relations is not affine! r="<<r<<"; col="<<col<<"; l0="<<resL0 << "; l1="<<resL1 << endl;
				return -1;
			}

			// Compute L = L0 * L1^{-1} = A0 * beta * A0^{-1}, section 3.2 and its
			// characteristic polynomial, what will be used later to determine its beta part.
			mat_GF2 L = L0.L * L1.Linv;
			GF2X Lpoly = characteristicPolynomial(L);
			int LpolyInt = getLong(Lpoly);

			// Structure needed in proposition3
			prop3struct_t * prop3struct = prop3structVector[(r-roundStart)*4+col];
			baseVectors_t bVectors;
			BYTE beta;

			// Determine beta values for L using characteristic polynomial - explained above
			if (charPolyMap.count(LpolyInt)==0){
				if (doCout) cout << "Error, cannot find beta for characteristic polynomial: " << Lpoly << endl;
				if (doCout) cout << "Attack cannot continue, please check implementation for possible bugs..." << endl;
				return -5;
			}

			bool betaValid=false;
			auto its = charPolyMap.equal_range(LpolyInt);
			for (auto it = its.first; it != its.second; ++it) {
				int prop3Result=-1;
				beta = it->second;
				if (doCout) cout << "Possible beta for L is " << CHEX(beta) << "; poly(L)=" << Lpoly << endl;

				// Try this value of beta with Proposition 2 to obtain A~_0 = A_0 * gamma for L=A0 * beta * A0^{-1}
				// From proposition 2 we will get A~0 = A0 * gamma, where A0 is linear part of Q0
				// Assuming we know L = A0 * beta * A0^{-1}  and beta.
				int AtildResult = proposition2(L, bVectors, defAES.makeMultAMatrix(beta));
				if (AtildResult < 0 || bVectors.size()==0){
					if (doCout) cout << "Something bad happen, attack does not work here - Atild computation problem" << endl;
					continue;
				}

				// Just debug for user - we may get some error if beta was not correctly guessed, or
				// subspace of solutions.
				betaValid = true;
				if (doCout) cout << "We have number of solutions=" << bVectors.size() << endl;

				// Continue with proposition 3 to extract info about affine part Q0.
				// Recovers deltaMult_i and c_i from following equations:
				// P~_0(x) = (SboxInv * deltMult_0 * AtildInv)(y_row( x, 00, 00, 00) + c0)
				// P~_1(x) = (SboxInv * deltMult_1 * AtildInv)(y_row(00,  x, 00, 00) + c1)
				// P~_2(x) = (SboxInv * deltMult_2 * AtildInv)(y_row(00, 00,  x, 00) + c2)
				// P~_3(x) = (SboxInv * deltMult_3 * AtildInv)(y_row(00, 00, 00,  x) + c3)
				//
				// Also holds the following:
				// P~_i = P_i(x) + k_i  , i \in [0,3]
				//     k_i is corresponding round key byte
				//     P_i(x) is affine input transformation (we changed it to affine in 1. phase)
				// Generally:
				// P~_i = L_i(x) + Pc + k_i
				//
				prop3Result = proposition3(prop3struct, defAES, r, col, 0, bVectors, vectorIdx);
				if (doCout) cout << "Prop3 done(v="<<vectorIdx<<") with result=" << CHEX(prop3Result) << "; gamma=" << CHEX(prop3struct->gamma) << endl;
				for(int x=0; doCout && x<4; x++){
					cout << "  Prop3["<<x<<"] valid="<<(((prop3Result & (1<<x))>0) ? 1:0)
						 << "; .valid=" << prop3struct->P[x].valid
						 << "; delta="<<CHEX(prop3struct->P[x].delta)
						 << "; c="<<CHEX(prop3struct->P[x].c)
						 << "; affConst="<<CHEX(prop3struct->P[x].affineConst)
						 << "; alfa_{0,"<<x<<"}="<<CHEX(prop3struct->P[x].alfa_0) << endl;
				}

				// If proposition 3 failed, reason may be invalid beta, so try another one.
				// From paper we know that our beta could be from set {b, b^2} so try next one.
				// But definitely one of betas should be correct.
				if (prop3Result != 0xf){
					if (doCout) cout << "This value of Beta was probably not valid, since proposition 3 failed for some relations" << endl;
					betaValid=false;
					continue;
				}
			}

			// Things didn't go as expected, beta was not correctly guessed from characteristics polynomial.
			// This should never happen if implementation is correct.
			if (!betaValid){
				if (doCout) cout << "Problem with beta, cannot continue with attack." << endl;
				return -2;
			}

			// Construct also q0 ( constant part of Q0 affine mapping )
			// q0 = c0 ^ c1 ^ c2 ^ c3 ^ c4
			// c4 = y0(00,00,00,00)
			//
			// Computing q0 is essential for round key extraction.
			// Recall that P,Q are matching and AFFINE transformations.
			// P~_i = P_i(x) + k_i ==> P~^{r+1}_i * Q^{r}_i = P_i(Q(x)) + k_i
			//
			// P~^{r+1}_i * Q^{r}_i (0) =		// P^{r+1}_i * Q^{r}_i = identity  (affine MATCHING)
			//   =  P_i(Ai * 0 + qi) + k_i
			//   =  P_i(qi) + k_i
			//   =  P_iA * qi + P_iC + k_i
			//   =  k_i
			Qaffine->qj[r][col] =
					prop3struct->c4
					^ prop3struct->P[0].c
					^ prop3struct->P[1].c
					^ prop3struct->P[2].c
					^ prop3struct->P[3].c;

			// delta_i^{-1} = gamma * alfa_{0,i}
			// Now we know MC coefficients alfa, 01 occurs twice => two values in delta_i will be same
			// so they correspond to 01 MC coefficients. Others can be derived by finding row from
			// MC that has two 01 on particular positions in multiplication with corresponding x_j
			//
			// Small redefinition was here from previous sections
			// A~_0 = A_0 * 1/gamma    (previously it was gamma)
			//
			// We know that last MC coef is 01 (described above) => delta_3^{-1} = gamma
			//
			// If we already know gamma, we can compute A0, linear part of affine output transformation Q0
			// A_0 = A~_0 * gamma
			//
			// Now we can derive A_0, linear part of Q0
			mat_GF2 gammaMatrix = defAES.makeMultAMatrix(prop3struct->gamma);

			// Thing mentioned in original paper here!
			// Proposition 2 can return subspace that solves given system of equations
			// but real A_0, real linear part of Q0 is only one, since Q0 is only one. Beta is determined correctly
			// from comparison of characteristic polynomial.
			// But it seems that:
			// A_0 = \text{\~{A}}^{(1)}_0 \cdot \gamma^{(1)} = \text{\~{A}}^{(2)}_0 \cdot \gamma^{(2)} = \dots = \text{\~{A}}^{(2^{\text{dim}})}_0 \cdot \gamma^{(2^{\text{dim}})}
			// A_0 = A~0 * gamma_0 = A'~0 * gamma'_0 = ... etc. for whole vector space.
			mat_GF2 A0Matrix(INIT_SIZE, 8, 8);
			mat_GF2 A0InvMatrix(INIT_SIZE, 8, 8);
			assert(vectorIdx==1); // here we just assume that vector from sulution subspace choosen is #1
			for(int x=0; x<64; x++){
				A0Matrix.put(x/8, x%8, bVectors[0].get(x));
			}

			A0Matrix = A0Matrix * gammaMatrix;
			Qaffine->Aj[r][col] = A0Matrix;

			// test A0 for invertibility, should be invertible
			GF2 determinant;
			NTL::inv(determinant, A0InvMatrix, A0Matrix);
			if(determinant==GF2::zero()){
				if (doCout) cout << "Problem with the attack, A0 (linear part of Q0) is not invertible" << endl;
				return -6;
			}

			// We can check if everything is OK by comparing now L = A0 * beta * A0^{-1},
			// Just one of the possible assertions that should hold.
			mat_GF2 betaMatrix = defAES.makeMultAMatrix(beta);
			mat_GF2 Lcons = A0Matrix * betaMatrix * A0InvMatrix;
			if (doCout && Lcons != L){
				cout << "L matrix determined in proposition 2 is not same as should be" << endl;
				cout << "L matrix: " << endl; dumpMatrix(L);
				cout << "A0 * beta * A0^{-1} matrix: " << endl; dumpMatrix(Lcons);
				cout << endl;
			}

			// From the Proposition 1, if we know already linear part A0 of Q0, we can this way determine
			// Aj linear part of Qj for other mappings
			//
			// Proposition 1 solves:
			// y_i(x0,00,00,00) = L(y_j(x0,00,00,00)) + c    while L and c is returned
			// 	 L   = A_i * (\alfa_{i,0} * \alfa_{j,0}^{-1}) * A_j^{-1}  while we know L, A_i, alfa...
			//   A_j = L^{-1} * A_i * (\alfa_{i,0} * \alfa_{j,0}^{-1})
			//
			// Then with the Proposition 3 we obtain also constant parts - plugging known matrices to prop3 we get c_j
			recoverQj(defAES, r, col, 1, A0Matrix, prop3struct->P[0].alfa_0, 0x1, Qaffine->Aj[r][4  + col], &(Qaffine->qj[r][4  + col]));
			recoverQj(defAES, r, col, 2, A0Matrix, prop3struct->P[0].alfa_0, 0x1, Qaffine->Aj[r][8  + col], &(Qaffine->qj[r][8  + col]));
			recoverQj(defAES, r, col, 3, A0Matrix, prop3struct->P[0].alfa_0, 0x3, Qaffine->Aj[r][12 + col], &(Qaffine->qj[r][12 + col]));

			// We could have removed affine parts Q0 and P0 from WBAES implementation
			// to verify our approach is correct.
			// k^{r}_{i,j} = P~^{3}_{i,j} * Q|^{2}_{i,j}  where Q|^{2}_{i,j} is affine part of Q
		}
	}

	// From proposition3 we obtained P~ affine equivalences
	// P~_i = P_i(x) + k_i   where k_i is corresponding byte of round key and P_i(x) is affine transformation
	// P^{r+1}_{i,j} * Q^{r}_{i,j} = identity  ==> derive round key by composition of those two.
	if (doCout) cout << endl << "Going to reconstruct encryption key from extracted round keys..." << endl;

	// Recover keys as described. Round keys are indexed by rows, but in AES implementation round keys are
	// indexed by columns, so just transpose before roundkey->encryption key transformation.
	for(int r=1; r<rounds2hack; r++){
		int roundBase = roundStart+r-1;
		for(int col=0; col<4; col++){
			for (int row=0; row<4; row++){
				int i  = 4*row + col;
				roundKeys[r-1][generator.shiftRows[i]] = prop3structVector[r*4 + col]->P[row].affineMap[ Qaffine->qj[roundBase][ generator.shiftRows[i] ] ];
			}
		}
	}

	// just simply printout
	for(int r=0; doCout && r<rounds2hack-1; r++){
		cout << "* Round keys extracted from the process, r=" << (r+roundStart+1) << endl;
		dumpVectorT(roundKeys[r], 16);
	}

	//
	// Recover encryption key from recovered round keys
	//
	if (doCout) cout << "Recovering cipher key from round keys..." << endl;
	recoverCipherKey(defAES, roundKeys, encKey);

	if (doCout) cout << "Final result: " << endl;
	dumpVector(encKey);

	// Just debug piece of code - was needed in proof that BGE attack works also
	// on Dual AES scheme.
	if (0){
		// DEBUG: dump determined relations
		for(int r=1; r<rounds2hack; r++){
			int roundBase = roundStart+r-1;
			for(int i=0;i<16;i++){
				cout << "## r="<<r<<"; i="<<i
						<<"; Qj="<<CHEX(Qaffine->qj[roundBase][i])
						<<"; Aj="<<hashMatrix(Qaffine->Aj[roundBase][i])
						<<"; Matrix="<<dumpMatrix2str(Qaffine->Aj[roundBase][i], false)
						<<endl;
			}
		}

		GenericAES dualAES;
		dualAES.initFromIndex(AES_IRRED_POLYNOMIALS-1, AES_GENERATORS-1);

		if (doCout) {
			cout << "T: " << endl; dumpMatrix(dualAES.T);
			cout << "Tinv: " << endl; dumpMatrixN(dualAES.Tinv);
			cout << "T*20: " << endl; dumpMatrixN(dualAES.T * Qaffine->Aj[2][0]);
			cout << "T*21: " << endl; dumpMatrixN(dualAES.T * Qaffine->Aj[2][1]);
			cout << "20*T: " << endl; dumpMatrixN(Qaffine->Aj[2][0]*dualAES.T);
			cout << "21*T: " << endl; dumpMatrixN(Qaffine->Aj[2][1]*dualAES.T);

			cout << "Tinv*20: " << endl; dumpMatrixN(dualAES.Tinv * Qaffine->Aj[2][0]);
			cout << "Tinv*21: " << endl; dumpMatrixN(dualAES.Tinv * Qaffine->Aj[2][1]);
			cout << "20*Tinv: " << endl; dumpMatrixN(Qaffine->Aj[2][0]*dualAES.Tinv);
			cout << "21*Tinv: " << endl; dumpMatrixN(Qaffine->Aj[2][1]*dualAES.Tinv);
		}
	}

	//
	// Free memory part
	//
	for(int x=0; x<rounds2hack*4; x++){
		delete prop3structVector[x];
	}

	delete[] Sr;
	delete Qaffine;
	return 0;
}


int BGEAttack::invertCipherTest(){
	GenericAES defAES;
	defAES.init(0x11B, 0x03);

#ifndef AES_BGE_ATTACK
	cerr << "Cannot proceed with attack if \"AES_BGE_ATTACK\" is not defined, we are missing required additions"<<endl;
	exit(1);
#endif

	WBAESGenerator generator;
	ExtEncoding    coding;

	if (doCout) cout << "Generating AES..." << endl;
	//bool encrypt = true;
	generator.useDualAESARelationsIdentity=true;	// this attack works only on basic form
	generator.useDualAESIdentity=true;
	generator.useDualAESSimpeAlternate=false;
	generator.useIO04x04Identity=false;
	generator.useIO08x08Identity=false;
	generator.useMB08x08Identity=false;
	generator.useMB32x32Identity=false;
	generator.generateExtEncoding(&coding, WBAESGEN_EXTGEN_ID);

	this->wbaes = new WBAES;
	generator.generateTables(GenericAES::testVect128_key, KEY_SIZE_16, this->wbaes, &coding, true);  if (doCout) cout << "AES ENC generated" << endl;
	generator.generateTables(GenericAES::testVect128_key, KEY_SIZE_16, this->wbaes, &coding, false); if (doCout) cout << "AES DEC generated" << endl;
	//int (&nextTbox)[N_BYTES]     = encrypt ? (shiftT2) : (shiftT2);  // attack is not yet implemented for decryption

	// WBAES changed to state with affine matching bijections at round boundaries.
	if (doCout) cout << "Going to test WBAES before modifying tables" << endl;
	generator.testComputedVectors(true, this->wbaes, &coding);

	// Encrypt test vector and then try to decrypt it by inverting encryption tables
	W128b plain, cipher, state;
	arr_to_W128b(GenericAES::testVect128_plain[0], 0, plain);
	arr_to_W128b(GenericAES::testVect128_plain[0], 0, state);
	arr_to_W128b(GenericAES::testVect128_cipher[0], 0, cipher);

	// encryption
	this->wbaes->encrypt(state);
	if (doCout) {
		cout << "=====================" << endl;
		cout << "InvertTest plaintext: " << endl;
		dumpW128b(plain);
		cout << endl;
		cout << "InvertTest ciphertext: " << endl;
		dumpW128b(cipher);
		cout << endl;
		cout << "Enc(plaintext_test): " << endl;
		dumpW128b(state);
		cout << endl;
	}

	time_t lastRound = time(NULL);
	time_t lastTm = time(NULL);
	time_t curTm;
	// Invertion with use of Rbox function
	// Iterate over each round, finding inversion manually.
	//   Goal is to find X such that Rbox(X) = Y where Y is current state array.
	//
	//   Search is done simultaneously by running each column over GF(2^8)^4, since each column
	//	 is independent in one round of each other.
	//
	for(int r=9; r>=0; r--){
		W128b prevState2, prevStateFinal;
		int foundCols=0;
		unsigned long long int col = 0;	// represents whole state array column
		int colMask2compute = 15; // start with full mask

		if (doCout) cout << "Inverting cipher; round=" << r << endl;
		curTm     = time(NULL);
		lastRound = time(NULL);

		// run over GF(2^8)^4, one column, finding inverse now
		int cnt=0;
		for(col=0; col < 4294967296L && foundCols < 4; col++, cnt++){
			// Initialize prevState from col iterator - each column is the same - simultaneous search
			INIT_W128B_COL(prevState2, 0, col);
			INIT_W128B_COL(prevState2, 1, col);
			INIT_W128B_COL(prevState2, 2, col);
			INIT_W128B_COL(prevState2, 3, col);

			// current progress monitoring
			if (cnt >= 65536){
				cnt = 0;
				curTm = time(NULL);
				if ((curTm - lastTm) > 30){
					if (doCout) cout << "...progress=" <<((((double)col) / 4294967296.0)*100) << " %" << endl;
					lastTm = curTm;
				}
			}

			// Run encryption - during invertion we are using only one type of tables, finding inverses
			// Caution: prevState2 is changed in this function, so we don't have source, but in majority
			// of cases it is not needed, we can obtain it again from col variable if needed (match will be found)
			this->Rbox(prevState2, true, r, true, colMask2compute);

			// find if one of columns holds state = Rbox(prevState2)
			for(int i=0; i<4; i++){
				if ((colMask2compute & (1<<i)) > 0 && ISEQ_W128B_COL(prevState2, i, state, i)){
					if (doCout) cout << "..found inverse for col="<<i<<"; foundCols="<<foundCols
							<<"; col=" << col
							<<"; done=" << ((((double)col) / 4294967296.0)*100) << " %" << endl;

					// copy matched inverse to appropriate field in prevStateFinal
					INIT_W128B_COL(prevStateFinal, i, col);
					foundCols+=1;

					colMask2compute &= ~(1<<i); // remove col from mask
				}
			}
		}

		// If reached this point, all inverses should have been found
		if (foundCols!=4){
			if (doCout) cout << "ERROR: foundCols!=4; " << foundCols << endl;
			return -1;
		}

		// Do not invert shiftRows if this is first round - we already have a result
		if (r>0){
			// Invert ShiftRows -> apply inversion of shiftrows operation on state array
			// InvShiftRows is shifting each row to the right, each column by its number
			for(int i=1; i<4; i++){
				for(int j=0; j<=i; j++){
					BYTE tmp = prevStateFinal.B[3 + 4*i];
					prevStateFinal.B[3 + 4*i] = prevStateFinal.B[2 + 4*i];;
					prevStateFinal.B[2 + 4*i] = prevStateFinal.B[1 + 4*i];;
					prevStateFinal.B[1 + 4*i] = prevStateFinal.B[0 + 4*i];;
					prevStateFinal.B[0 + 4*i] = tmp;
				}
			}
		}

		time_t curTime = time(NULL);
		if (doCout) cout << "Inverse found; r="<<r
				<<"; time elapsed=" << (curTime - lastRound)
				<<" s; state dump: " << endl;

		dumpW128b(prevStateFinal);

		// copy back to state
		W128CP(state, prevStateFinal);
	}


	delete this->wbaes;
	return 0;
}

const vec_GF2E &BGEAttack::getEncKey() const {
    return encKey;
}

void BGEAttack::setEncKey(const vec_GF2E &encKey) {
    BGEAttack::encKey = encKey;
}

bool BGEAttack::isDoCout() const {
    return doCout;
}

void BGEAttack::setDoCout(bool doCout) {
    BGEAttack::doCout = doCout;
}

} /* namespace attack */
} /* namespace wbacr */
