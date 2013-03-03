/*
 * GenericAES.cpp
 *
 *  Created on: Mar 1, 2013
 *      Author: ph4r05
 */

#include "GenericAES.h"
#include <iomanip>
using namespace std;
using namespace NTL;

GenericAES::GenericAES() {
	// TODO Auto-generated constructor stub

}

GenericAES::~GenericAES() {
	// TODO Auto-generated destructor stub
}


mat_GF2 GenericAES::getDefaultAffineMatrix (void){
	mat_GF2 ret(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);
	VectorCopy(ret[0], GF2XFromLong(0x8F, AES_FIELD_DIM), AES_FIELD_DIM); // 1 0 0 0 1 1 1 1
	VectorCopy(ret[1], GF2XFromLong(0xC7, AES_FIELD_DIM), AES_FIELD_DIM); // 1 1 0 0 0 1 1 1
	VectorCopy(ret[2], GF2XFromLong(0xE3, AES_FIELD_DIM), AES_FIELD_DIM); // 1 1 1 0 0 0 1 1
	VectorCopy(ret[3], GF2XFromLong(0xF1, AES_FIELD_DIM), AES_FIELD_DIM); // 1 1 1 1 0 0 0 1
	VectorCopy(ret[4], GF2XFromLong(0xF8, AES_FIELD_DIM), AES_FIELD_DIM); // 1 1 1 1 1 0 0 0
	VectorCopy(ret[5], GF2XFromLong(0x7C, AES_FIELD_DIM), AES_FIELD_DIM); // 0 1 1 1 1 1 0 0
	VectorCopy(ret[6], GF2XFromLong(0x3E, AES_FIELD_DIM), AES_FIELD_DIM); // 0 0 1 1 1 1 1 0
	VectorCopy(ret[7], GF2XFromLong(0x1F, AES_FIELD_DIM), AES_FIELD_DIM); // 0 0 0 1 1 1 1 1

	// fix: since NTL uses different order of elements (1+x^2+x^3) vs. (x^3+x^2+1),
	// we have to reverse now. 0x8F === 11110001, not 10001111 what we would expect,
	// but matrix multiplication works as usually.
	reverse(ret[0], ret[0]); reverse(ret[1], ret[1]);
	reverse(ret[2], ret[2]); reverse(ret[3], ret[3]);
	reverse(ret[4], ret[4]); reverse(ret[5], ret[5]);
	reverse(ret[6], ret[6]); reverse(ret[7], ret[7]);
	return ret;
}

/**
 * Affine constant using during encryption
 */
vec_GF2 GenericAES::getDefaultAffineConst (void){
	vec_GF2 ret(INIT_SIZE, AES_FIELD_DIM);
	VectorCopy(ret, GF2XFromLong(0x63, AES_FIELD_DIM), AES_FIELD_DIM);
	//reverse(ret, ret);
	return ret;
}

mat_GF2 GenericAES::getDefaultAffineMatrixDec (void){
	mat_GF2 ret(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);
	VectorCopy(ret[0], GF2XFromLong(0x25, AES_FIELD_DIM), AES_FIELD_DIM); // 0 0 1 0 0 1 0 1
	VectorCopy(ret[1], GF2XFromLong(0x92, AES_FIELD_DIM), AES_FIELD_DIM); // 1 0 0 1 0 0 1 0
	VectorCopy(ret[2], GF2XFromLong(0x49, AES_FIELD_DIM), AES_FIELD_DIM); // 0 1 0 0 1 0 0 1
	VectorCopy(ret[3], GF2XFromLong(0xA4, AES_FIELD_DIM), AES_FIELD_DIM); // 1 0 1 0 0 1 0 0
	VectorCopy(ret[4], GF2XFromLong(0x52, AES_FIELD_DIM), AES_FIELD_DIM); // 0 1 0 1 0 0 1 0
	VectorCopy(ret[5], GF2XFromLong(0x29, AES_FIELD_DIM), AES_FIELD_DIM); // 0 0 1 0 1 0 0 1
	VectorCopy(ret[6], GF2XFromLong(0x94, AES_FIELD_DIM), AES_FIELD_DIM); // 1 0 0 1 0 1 0 0
	VectorCopy(ret[7], GF2XFromLong(0x4A, AES_FIELD_DIM), AES_FIELD_DIM); // 0 1 0 0 1 0 1 0

	// fix: since NTL uses different order of elements, we have to reverse now
	// 0x8F === 11110001, not 10001111 what we would expect, but matrix multiplication
	// works as usually.
	reverse(ret[0], ret[0]); reverse(ret[1], ret[1]);
	reverse(ret[2], ret[2]); reverse(ret[3], ret[3]);
	reverse(ret[4], ret[4]); reverse(ret[5], ret[5]);
	reverse(ret[6], ret[6]); reverse(ret[7], ret[7]);
	return ret;
}

/**
 * Affine constant using during decryption
 */
vec_GF2 GenericAES::getDefaultAffineConstDec (void){
	vec_GF2 ret(INIT_SIZE, AES_FIELD_DIM);
	VectorCopy(ret, GF2XFromLong(0x05, AES_FIELD_DIM), AES_FIELD_DIM);
	//reverse(ret, ret);
	return ret;
}

/**
 * Generates whole AES table with generator
 */
void GenericAES::build() {
	modulusContext.save();

	// 1. generate whole field with given generator
	int i=0,c=0;
	generator.LoopHole().SetLength(8);
	generator.LoopHole().SetMaxLength(8);

	// current element of the field
	GF2E cur;
	cur.LoopHole().SetMaxLength(8);
	cur.LoopHole().SetLength(8);
	//cur.init(this->modulus);

	cur=1;
	this->gInv[0] = -1;
	for(i=0; i<AES_FIELD_SIZE; i++){
		this->g[i] = cur;
		this->gInv[getLong(cur)] = i;
		cur = cur * this->generator;

		this->g[i].LoopHole().SetLength(AES_FIELD_DIM);
		this->g[i].LoopHole().SetMaxLength(AES_FIELD_DIM);
	}

	// 2. compute inverses in terms of generator exponent
	this->sbox[0] = -1;
	for(i=1; i<AES_FIELD_SIZE; i++){
		this->sbox[i] = 255-i;
	}

	// 3. generate base change matrix [g^0, g^25, g^50, g^75, g^100, g^125, g^150, g^175] == polynomial
	// each element is
	T = mat_GF2(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);
	for(c=0, i=0; c<8; i+=25, c++){
		vec_GF2 tmp = VectorCopy(g[i].LoopHole(), AES_FIELD_DIM);
		for(int j=0; j<AES_FIELD_DIM; j++){
			T.put(j,c, tmp.get(j));
		}
	}

	// 4. invert base change matrix
	Tinv = inv(T);

	// 5. S-BOX with affine mappings
	mat_GF2 tmpSboxAffMatrix(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);	// T * Affine * Tinv
	mat_GF2 tmpSboxAffConst(INIT_SIZE, AES_FIELD_DIM, 1);				// column vector
	mat_GF2 tmpSboxAffMatrixDec(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);	// T * AffineDec * Tinv
	mat_GF2 tmpSboxAffConstDec(INIT_SIZE, AES_FIELD_DIM, 1);				// column dec vector

	tmpSboxAffMatrix = T * getDefaultAffineMatrix() * Tinv;
	tmpSboxAffConst = T * colVector(getDefaultAffineConst());
	tmpSboxAffMatrixDec = T * getDefaultAffineMatrixDec() * Tinv;
	tmpSboxAffConstDec = T * colVector(getDefaultAffineConstDec());


	cout << "Encryption const: " << tmpSboxAffConst << "; encMat: " << endl << tmpSboxAffMatrix << endl;
	cout << "Decryption const: " << tmpSboxAffConstDec << "; decMat: " << endl << tmpSboxAffMatrixDec << endl;
	cout << "Enc * Dec : " << endl << (tmpSboxAffMatrix * tmpSboxAffMatrixDec) << endl;
	cout << "Dec * Enc : " << endl << (tmpSboxAffMatrixDec * tmpSboxAffMatrix) << endl;
	cout << "EncInv: " << endl << (inv(tmpSboxAffMatrix)) << endl;
	cout << "DecInv: " << endl << (inv(tmpSboxAffMatrixDec)) << endl;

	// default cases for zero;
	for(i=0; i<AES_FIELD_SIZE; i++){
		GF2E transValue = i==0 ? GF2E::zero() : g[255-i];
		int target = getLong(transValue);
		mat_GF2 resMatrix = tmpSboxAffMatrix * colVector(rep(transValue), AES_FIELD_DIM) + tmpSboxAffConst;
		colVector(transValue, resMatrix, 0);
		this->sboxAffine[target] = getLong(transValue);

		GF2E transValue2 = i==0 ? GF2E::zero() : g[255-i];
		resMatrix = tmpSboxAffMatrixDec * colVector(rep(transValue2), AES_FIELD_DIM) + tmpSboxAffConstDec;
		colVector(transValue2, resMatrix, 0);
		this->sboxAffineInv[target] = getLong(transValue2);
	}

	// 6. MixColumn operations
	mixColModulus.SetLength(5);
	mixColModulus[0] = g[0];
	mixColModulus[4] = g[0];

	mixColMultiply.SetLength(4);
	mixColMultiply[0] = g[1];
	mixColMultiply[3] = g[25];

	mixColMultiplyInv.SetLength(4);
	mixColMultiplyInv[0] = g[104];
	mixColMultiplyInv[1] = g[228];
	mixColMultiplyInv[2] = g[199];
	mixColMultiplyInv[3] = g[223];
	return;
}

void GenericAES::printAll() {
	int i=0,c=0;
	cout << "Generic AES; " \
			<< "modulus[" << this->modulus << "] "
			<< "generator[" << this->generator << "]"<<endl;

	// dump fields
	cout << "Dumping field generated: " << endl;
	for(i=0; i<AES_FIELD_SIZE; i++){
		// TODO: refactor this formating to use cout
		printf("  g[%03d] -> [%02X]   repr: ", i, (int)getLong(this->g[i]));
		cout << this->g[i] << endl;
	}

	// dump inverse
	cout << "Dumping field generated: " << endl;
	for(i=1; i<AES_FIELD_SIZE; i++){
		// TODO: refactor this formating to use cout
		printf("  S[%03d] -> g[%03d] = [%02X] repr: ", i, (int) sbox[i], (int)getLong(g[sbox[i]]));
		cout << g[sbox[i]] << " test: " << (g[i] * g[sbox[i]]) << endl;
	}

	cout << "Dumping Sboxes: " << endl;
	for(i=0; i<AES_FIELD_SIZE; i++){
		printf("  Sbox  (0x%02X)=%03d -> %02X repr: ", i, i, (int)sboxAffine[i]);
		cout << sboxAffineInv[sboxAffine[i]] << endl;
	}

	cout << "Base change matrix: " << endl << T << endl;
	cout << "Base change ivn matrix: " << endl << Tinv << endl;

	cout << "MixCol mult " << mixColMultiply << endl;
	cout << "MixCol multInv " << mixColMultiplyInv << endl;
	cout << "MixCol mod " << mixColModulus << endl;

	return;
}


