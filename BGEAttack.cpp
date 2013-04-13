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
	;

}

BGEAttack::~BGEAttack() {
	;
}

using namespace std;
using namespace NTL;
using namespace boost;
using namespace wbacr::attack;
using namespace wbacr::laeqv;

int BGEAttack::shiftIdentity[16] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
int BGEAttack::shiftT2[16] =
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

// Rbox operation - one round of AES. State is indexed by rows
// First row: 0 1 2 3; It is processed by columns, first column is computed from
// state array indexes: 0,4,8,12; second from 1,5,9,13, and so on...
void BGEAttack::Rbox(W128b& state, bool encrypt, int r, bool noShift){
	int i=0;
	W32b ires[N_BYTES];				// intermediate result for T-boxes

	// encryption/decryption dependent operations and tables
	int (&shiftOp)[N_BYTES]  = noShift ? this->shiftIdentity : (encrypt ? (this->wbaes.shiftRows)   : (this->wbaes.shiftRowsInv));
	W32XTB (&edXTab)[N_ROUNDS][N_SECTIONS][N_XOR_GROUPS] = encrypt ? (this->wbaes.eXTab) 	   : (this->wbaes.dXTab);
	AES_TB_TYPE2 (&edTab2)[N_ROUNDS][N_BYTES]			 = encrypt ? (this->wbaes.eTab2) 	   : (this->wbaes.dTab2);
	AES_TB_TYPE3 (&edTab3)[N_ROUNDS][N_BYTES]			 = encrypt ? (this->wbaes.eTab3) 	   : (this->wbaes.dTab3);
#ifdef AES_BGE_ATTACK
	GF256_func_t (&edOutputBijection)[N_ROUNDS][N_BYTES] = encrypt ? (this->wbaes.eOutputBijection) : (this->wbaes.dOutputBijection);
#endif

	// Perform rest of the operations on 4 tuples.
	for(i=0; i<N_BYTES; i+=4){
		// Apply type 2 tables to all bytes, counting also shift rows selector.
		// One section ~ 1 column of state array, so select 1 column, first will
		// have indexes 0,4,8,12. Also take ShiftRows() into consideration.
		ires[i+0].l = edTab2[r][i+0][state.B[shiftOp[i/4+0*4]]].l;
		ires[i+1].l = edTab2[r][i+1][state.B[shiftOp[i/4+1*4]]].l;
		ires[i+2].l = edTab2[r][i+2][state.B[shiftOp[i/4+2*4]]].l;
		ires[i+3].l = edTab2[r][i+3][state.B[shiftOp[i/4+3*4]]].l;

		// In the last round, result is directly in T2 boxes
		if (r==(N_ROUNDS-1)){
			continue;
		}

		// XOR results of T2 boxes
		op8xor(ires[i+0], ires[i+1], edXTab[r][i/4][0], ires[i+0]);  // 1 xor 2
		op8xor(ires[i+2], ires[i+3], edXTab[r][i/4][1], ires[i+2]);  // 3 xor 4
		op8xor(ires[i+0], ires[i+2], edXTab[r][i/4][2], ires[i+0]);  // (1 xor 2) xor (3 xor 4) - next XOR stage

		// Apply T3 boxes, valid XOR results are in ires[0], ires[4], ires[8], ires[12]
		// Start from the end, because in ires[i] is our XORing result.
		//
		//                    ________________________ ROUND
		//                   |    ____________________ T3 box for 1 section
		//                   |   |      ______________ (1 xor 2) xor (3 xor 4)
		//                   |   |     |        ______ 8bit parts of 32 bit result
		//                   |   |     |       |
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
		state.B[i/4+ 0] = r<(N_ROUNDS-1) ? ires[i].B[0] : ires[i+0].B[0];
		state.B[i/4+ 4] = r<(N_ROUNDS-1) ? ires[i].B[1] : ires[i+1].B[0];
		state.B[i/4+ 8] = r<(N_ROUNDS-1) ? ires[i].B[2] : ires[i+2].B[0];
		state.B[i/4+12] = r<(N_ROUNDS-1) ? ires[i].B[3] : ires[i+3].B[0];
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
				cout << "i="<<i<<"; Beta=["<<CHEX(intBeta)
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
			cout << "i="<<i<<"; Beta=["<<CHEX(intBeta)
					<<"]; a="<<CHEX(coefs[a].intCoef)
					<<"; b="<<CHEX(coefs[b].intCoef)
					<<"; c="<<CHEX(coefs[c].intCoef) << ", inv="<<CHEX(coefs[c].intCoefInv)
					<<"; d="<<CHEX(coefs[d].intCoef) << ", inv="<<CHEX(coefs[d].intCoefInv)<<endl;

			if (bset.count(intBeta)==0) bset.insert(intBeta);
		}
	}

	return 0;
}

NTL::GF2X BGEAttack::characteristicPolynomial(mat_GF2X_t m){
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

NTL::GF2X BGEAttack::characteristicPolynomial(mat_GF2 m){
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
	for(c=0; c<0xff; c++){
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
			colVector(lei, L.L, i);
		}

		// If transformation L is linear, then it has to have proper matrix inverse and determinant!=0
		GF2 determinant;
		NTL::inv(determinant, L.Linv, L.L);
		if (determinant==GF2::zero()){
			continue;	// matrix is not invertible -> not linear
		}

		// Now check whether already determined mapping is linear for whole space, just check.
		bool works = true;
		for(int x=0; x<0xff; x++){
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
			for(int i=63-freeVariables, c=0; i>=0; i--, c++){
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
				outM.put(i/8, i%8, x[i]);
			}

			// Final verification of given solution - is equation on the top correct now?
			mat_GF2 lhs = L * outM;
			mat_GF2 rhs = outM * beta;
			if(lhs != rhs){
				cout << "Something wrong with the result, equation does not hold in self-test..." << endl;
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
		cout << "Something wrong with the result, equation does not hold..." << endl;
		return -2;
	}

	out.push_back(x);
	return 0;
}

int BGEAttack::run(void) {
	GenericAES defAES;
	defAES.init(0x11B, 0x03);

#ifndef AES_BGE_ATTACK
	cerr << "Cannot proceed with attack if \"AES_BGE_ATTACK\" is not defined, we are missing required additions"<<endl;
	exit(1);
#endif

	WBAESGenerator generator;
	CODING8X8_TABLE coding[16];
	W128b state;
	cout << "Generating AES..." << endl;
	bool encrypt = true;
	generator.useDualAESARelationsIdentity=true;	// this attack works only on basic form
	generator.useDualAESIdentity=true;
	generator.generateIO128Coding(coding, true);
	generator.generateTables(GenericAES::testVect128_key, KEY_SIZE_16, this->wbaes, coding, true);  cout << "AES ENC generated" << endl;
	generator.generateTables(GenericAES::testVect128_key, KEY_SIZE_16, this->wbaes, coding, false); cout << "AES DEC generated" << endl;
	int (&nextTbox)[N_BYTES]     = encrypt ? (generator.shiftRowsLBijection) : (generator.shiftRowsLBijectionInv);

	// WBAES changed to state with affine matching bijections at round boundaries.
	cout << "Going to test WBAES before modifying tables" << endl;
	generator.testComputedVectors(true, this->wbaes, coding);

	//
	//
	// Attack below...
	//
	//

	// Recover affine parts of Q for round r=0
	int r = 0;
	int row=0;
	int x = 0;
	int i = 0;
	int j = 0;
	int c1;

	// At first compute base function f_{00}, we will need it for computing all next functions,
	// to be exact, its inverse, f^{-1}_{00}
	cout << "Allocating memory for the attack" << endl;
	SsetPerRound_t * Sr = new SsetPerRound_t[10];
	Qaffine_t * Qaffine = new Qaffine_t;
	cout << "Memory allocated; Sr=" << dec << (sizeof(SsetPerRound_t)*10) << "; Qaffine=" << dec << sizeof(Qaffine_t) << endl;
	cout << "Memory allocated totally: " << (((sizeof(SsetPerRound_t)*10) + sizeof(Qaffine_t)) / 1024.0 / 1024.0) << " MB" << endl;

	cout << "Starting attack phase 1 ..." << endl;
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
		cout << "Generating f_00 for round r="<<r<<endl;
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

		// f(x,0,0,0) finalization - compute hash
		for(i=0; i<AES_BYTES; i++){
			Sr[r].S[i%4][i/4].f_00.initHash();
		}

		// now generate f(x,0,0,0) .. f(x,0xff,0,0) functions, generate then corresponding sets and whole Sr for round r
		cout << "Generating set S..." << endl;
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

				// By iterating c1 over GF256 we constructed vector space, with 256 elements, thus every
				// should be different. We are checking, whether this assumption is OK and vector space is correct.
				if (Sr[r].S[i%4][i/4].fmap.count(hash)!=0){
					cerr << "Sr["<<r<<"].["<<(i%4)<<"]["<<(i/4)<<"] set already contains hash[" << hash << "]; for j=" << j << std::endl;
					assert(Sr[r].S[i%4][i/4].fmap.count(hash)!=0);
				}

				Sr[r].S[i%4][i/4].fmap.insert(fctionMap_elem_t(hash, j));
			}
		}



		// Now we have Sr[r] generated.
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
			cout << "Recovering psi for Sidx="<<i<<endl;
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

				// correctness
				if (keySet.count(x)>0){
					cerr << "Key x=["<<x<<"] already defined in Q["<<r<<"]["<<i<<"].f["<<x<<"]"<< endl;
				}
				if (valSet.count(y)>0){
					cerr << "Val y=["<<y<<"] already defined in Q["<<r<<"]["<<i<<"].f["<<x<<"]"<< endl;
				}

				Qaffine->Q[r][i].f[x]    = y;
				Qaffine->Q[r][i].finv[y] = x;
				keySet.insert(x);
				valSet.insert(y);
				//cout << "Q["<<r<<"]["<<i<<"].f["<<x<<"] = " << y << endl;
			}

			Qaffine->Q[r][i].initHash();
			cout << "Q~ recovered; hash=" << Qaffine->Q[r][i].hash << endl;
		}

		cout << "PSI recovered for all sets in given round" << endl;
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
	// Q will be on transfered to AES output bijection layer - additional one
	// P will be on this->wbaes.eTab2[r][?][0..256]
	cout << "Going to transform AES input/output matching bijection to affine" << endl;
	for(r=0; r<N_ROUNDS-1; r++){
		for(i=0; i<AES_BYTES; i++){
			// Choice of indexes for eTab2 explained below:
			// Since we are interested in output encoding, we don't really care about input to Rbox (shiftrows is ignored).
			// We have Q[r][i \ in [1..16]] which corresponds to state array indexed by rows.
			// In WBAES implementation, shift rows will follow, thus first 4 T2boxes will be feeded by indexes of state array:
			// state[0,7,10,13]; next 4 T2 boxes will be feeded by state[4,1,14,11]
			// The point is that nextTbox \ocirc shiftRows (during evaluation) = identity
			AES_TB_TYPE2 & curT2     = this->wbaes.eTab2[(r+1)][shiftT2[i]];     // T2 table in next round connected to particular HILO(xtb1, xtb2)
			XTB & curXtb1            = this->wbaes.eXTab[r][i%4][5][2*(i/4)+0];  // 3.rd index - XOR table number 5, last in round
			XTB & curXtb2            = this->wbaes.eXTab[r][i%4][5][2*(i/4)+1];
			GF256_func_t & outBiject = this->wbaes.eOutputBijection[r][i];

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
				for(x=0; x<GF256; x++){
					BYTE newXtbX = HILO(curXtb1[x], curXtb2[y]);
					outBiject[newXtbX] = bijectionFinv[newXtbX];	// this will be repeated 256 times, but with same value, so OK
				}

				// For updating T2 boxes we only need 1 for loop iterating over field since function is 8->8
				curT2[y].l = tmpT2[bijectionF[y]].l;
			}
		}
	}

	// WBAES changed to state with affine matching bijections at round boundaries.
	cout << "Going to test WBAES after modifying tables" << endl;
	generator.testComputedVectors(true, this->wbaes, coding);

	//
	// Attack, proceeding to phase 2 - removing affine parts for output encoding
	//

	// Derive B set
	Bset bset2;
	deriveBset(bset2, defAES, true);

	cout << "Bset: { ";
	for(Bset::const_iterator it = bset2.begin(); it!=bset2.end(); ++it){
		cout << CHEX(*it) << ", ";
	} cout << "} " << endl;

	// Now we can generate multiplication matrices from given constants and then compute char. polynomials for them.
	// L = A0 * beta * A0^{-1} is similar to beta matrix, so they have same characteristic polynomials. From
	// this we can determine possible beta values for given L as matrix.
	charPolynomialMultimap_t charPolyMap;
	for(Bset::const_iterator it = bset2.begin(); it!=bset2.end(); ++it){
		mat_GF2 m = defAES.makeMultAMatrix(*it);
		cout << "Computing characteristic polynomial for: " << CHEX(*it) << "; ....";

		GF2X poly = characteristicPolynomial(m);
		cout << "Characteristic polynomial=" << poly << "; hex=" << GF2EHEX(poly) << endl;

		// Insert characteristic polynomial to multi-map. Later it will be used to determine B
		charPolyMap.insert(charPolynomialMultimap_elem_t(getLong(poly), *it));
	}

	// Now we can compute L0, L1 according to the paper from proposition1
	affineEquiv_t L0, L1;
	int resL0 = proposition1(L0, 2, 0, 0, 1, 0);
	int resL1 = proposition1(L1, 2, 0, 0, 1, 1);
	if (resL0==0 || resL1==0){
		cout << "One of relations is not affine! l0="<<resL0 << "; l1="<<resL1 << endl;
		return -1;
	}

	// Compute L = L0 * L1^{-1}
	mat_GF2 L = L0.L * L1.Linv;
	GF2X Lpoly = characteristicPolynomial(L);
	int LpolyInt = getLong(Lpoly);

	// determine beta values for L
	if (charPolyMap.count(LpolyInt)==0){
		cout << "Error, cannot find beta for characteristic polynomial: " << Lpoly << endl;
	} else {
		bool betaValid=false;
		auto its = charPolyMap.equal_range(LpolyInt);
		for (auto it = its.first; it != its.second; ++it) {
			BYTE possBeta = it->second;

		//for(auto it = bset2.begin(); it!=bset2.end(); ++it){
			//BYTE possBeta = *it;

			cout << "Possible beta for L is " << CHEX(possBeta) << "; poly(L)=" << Lpoly << endl;

			// Try this value of beta, proceed with proposition 2 to obtain A~ for L
			baseVectors_t bVectors;
			int AtildResult = proposition2(L, bVectors, defAES.makeMultAMatrix(possBeta));
			if (AtildResult < 0){
				cout << "Something bad happen, attack does not work here - Atild computation problem" << endl;
				continue;
			}

			betaValid = true;
			cout << "We have solutions: " << endl;
			for(int k=0, ln=bVectors.size(); k<ln; k++){
				dumpVector(bVectors[k]);
			}
		}

		// Things didn't go as expected.
		if (!betaValid){
			cout << "Problem with beta, cannot continue with attack." << endl;
			return -2;
		}
	}

	delete[] Sr;
	delete Qaffine;


	return 0;
}

} /* namespace attack */
} /* namespace wbacr */
