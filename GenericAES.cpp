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

//#define GENERIC_AES_DEBUG 1


GenericAES::GenericAES() {
	// TODO Auto-generated constructor stub

}

GenericAES::~GenericAES() {
	// TODO Auto-generated destructor stub
}

/**
 * Initializes AES with given modulus and generator in long representation
 */
void GenericAES::init(long modulus, long generator){
	// get previous context and save it - will be later restored
	//GF2EContext oldContext;
	//oldContext.save();

	// set modulus to internal attribute
	this->modulus = GF2XFromLong(modulus, AES_FIELD_DIM+1);
	// initialize modulus to GF2E class
	GF2E::init(this->modulus);
	// store
	this->modulusContext.save();
	this->generator = GF2EFromLong(generator, AES_FIELD_DIM);
	this->build();

	// restore old context
	//oldContext.restore();
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
	return ret;
}

/**
 * Expand key to round keys
 */
void GenericAES::expandKey(vec_GF2E& expandedKey, vec_GF2E& key, enum keySize size){
	/* current expanded keySize, in bytes */
	unsigned int currentSize = 0;
	unsigned int rconIteration = 0;
	unsigned int i,j;
	unsigned int Nr = GenericAES::getNumberOfRounds(size);
	unsigned int expandedKeySize = (4 * AES_NB * (Nr + 1));
	restoreModulus();

	GF2E tmp;
	vec_GF2E t(INIT_SIZE, 4);

	// result expanded key size = NB * (NK + 1)
	#ifdef GENERIC_AES_DEBUG
	cout << "Expanded key size will be: " << expandedKeySize << endl;
	#endif
	expandedKey.SetLength(expandedKeySize);

	/* set the 16,24,32 bytes of the expanded key to the input key */
	for (i = 0; i < size; i++)
		expandedKey[i] = key[i];

	currentSize += size;
	while (currentSize < expandedKeySize){
		#ifdef GENERIC_AES_DEBUG
		cout << "CurrentSize: " << currentSize <<  "; expandedKeySize: " << expandedKeySize <<  endl;
		#endif
		/* assign the previous 4 bytes to the temporary value t */
		for (i = 0; i < 4; i++) {
			t[i] = expandedKey[(currentSize - 4) + i];
		}

		/**
		 * every 16,24,32 bytes we apply the core schedule to t
		 * and increment rconIteration afterwards
		 */
		if(currentSize % size == 0) {
			//core(t, rconIteration++);
			/* rotate the 32-bit word 8 bits to the left */
			tmp=t[0]; 	t[0]=t[1];
			t[1]=t[2];	t[2]=t[3];
			t[3]=tmp;
			/* apply S-Box substitution on all 4 parts of the 32-bit word */
			for (j = 0; j < 4; ++j){
				#ifdef GENERIC_AES_DEBUG
				cout << "Sboxing key t[" << j << "]=" << t[j] << "=" << GF2EHEX(t[j]) << "; sboxval: " << CHEX(sboxAffine[getLong(t[j])]) ;
				#endif
				t[j] = GF2EFromLong(sboxAffine[getLong(t[j])], AES_FIELD_DIM);
				#ifdef GENERIC_AES_DEBUG
				cout << " after Sbox = " << t[j] << "="  << GF2EHEX(t[j]) << endl;
				#endif
			}
			/* XOR the output of the rcon operation with i to the first part (leftmost) only */
			t[0] = t[0] + RC[rconIteration++];

			#ifdef GENERIC_AES_DEBUG
			cout << "; after XOR with RC[" << GF2EHEX(RC[rconIteration-1]) << "] = " << t[0] << " = " << GF2EHEX(t[0]) << endl;
			#endif
		}

		/* For 256-bit keys, we add an extra sbox to the calculation */
		if(size == KEY_SIZE_32 && ((currentSize % size) == 16)) {
			for(i = 0; i < 4; i++)
				t[i] = GF2EFromLong(sboxAffine[getLong(t[i])], AES_FIELD_DIM);
		}
		/* We XOR t with the four-byte block 16,24,32 bytes before the new expanded key.
		* This becomes the next four bytes in the expanded key.
		*/
		for(i = 0; i < 4; i++) {
			expandedKey[currentSize] = expandedKey[currentSize - size] + t[i];

			#ifdef GENERIC_AES_DEBUG
			cout << "t[" << i << "] = " << GF2EHEX(t[i]) << endl;
			#endif

			currentSize++;
		}
	}
}

void GenericAES::applyT(mat_GF2E& state){
	int i,j,m=state.NumRows(),n=state.NumCols();
	for(i=0; i<m; i++){
		for(j=0;j<n;j++){
			applyT(state[i][j]);
		}
	}
}

void GenericAES::applyT(vec_GF2E& state){
	int i, ln=state.length();
	for(i=0; i<ln; i++){
		applyT(state[i]);
	}
}

void GenericAES::applyT(GF2E& state){
	mat_GF2 tmpMat = colVector(state, AES_FIELD_DIM);
	mat_GF2 resMat = T * tmpMat;
	colVector(state, resMat, 0);
}

void GenericAES::applyTinv(mat_GF2E& state){
	int i,j,m=state.NumRows(),n=state.NumCols();
	for(i=0; i<m; i++){
		for(j=0;j<n;j++){
			applyTinv(state[i][j]);
		}
	}
}
void GenericAES::applyTinv(vec_GF2E& state){
	int i, ln=state.length();
	for(i=0; i<ln; i++){
		applyTinv(state[i]);
	}
}

void GenericAES::applyTinv(GF2E& state){
	mat_GF2 tmpMat = colVector(state, AES_FIELD_DIM);
	mat_GF2 resMat = Tinv * tmpMat;
	colVector(state, resMat, 0);
}

/**
 * Generates whole AES table with generator
 */
void GenericAES::build() {
	modulusContext.save();

	// 1. generate whole field with given generator under given modulus
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
	}

	// 2. compute GF(256) element inverses in terms of generator exponent
	this->sbox[0] = -1;
	for(i=1; i<AES_FIELD_SIZE; i++){
		this->sbox[i] = 255-i;
	}

	// 3. generate base change matrix [g^0, g^25, g^50, g^75, g^100, g^125, g^150, g^175] == polynomial
	// what corresponds to [x^0 x^1 x^2 x^3 x^4 x^5 x^6 x^7]. Base is expressed in terms of generator.
	T = mat_GF2(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);
	for(c=0, i=0; c<8; i+=25, c++){
		vec_GF2 tmp = VectorCopy(g[i].LoopHole(), AES_FIELD_DIM);
		for(int j=0; j<AES_FIELD_DIM; j++){
			T.put(j,c, tmp.get(j));
		}
	}

	// 4. invert base change matrix
	Tinv = inv(T);

	// 5. S-BOX with affine mappings. At first obtain default form of affine transformation for normal AES
	mat_GF2 tmpSboxAffMatrix(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);		// T * Affine * Tinv
	mat_GF2 tmpSboxAffConst(INIT_SIZE, AES_FIELD_DIM, 1);					// column vector
	mat_GF2 tmpSboxAffMatrixDec(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);	// T * AffineDec * Tinv
	mat_GF2 tmpSboxAffConstDec(INIT_SIZE, AES_FIELD_DIM, 1);				// column dec vector
	// Now transform affine transformation with T, Tinv matrices to work in generic AES
	tmpSboxAffMatrix = T * getDefaultAffineMatrix() * Tinv;
	tmpSboxAffConst = T * colVector(getDefaultAffineConst());
	tmpSboxAffMatrixDec = T * getDefaultAffineMatrixDec() * Tinv;
	tmpSboxAffConstDec = T * colVector(getDefaultAffineConstDec());

	// DEBUG
#ifdef GENERIC_AES_DEBUG
	cout << "Encryption const: " << tmpSboxAffConst << "; encMat: " << endl << tmpSboxAffMatrix << endl;
	cout << "Decryption const: " << tmpSboxAffConstDec << "; decMat: " << endl << tmpSboxAffMatrixDec << endl;
	cout << "Enc * Dec : " << endl << (tmpSboxAffMatrix * tmpSboxAffMatrixDec) << endl;
	cout << "Dec * Enc : " << endl << (tmpSboxAffMatrixDec * tmpSboxAffMatrix) << endl;
	cout << "EncInv: " << endl << (inv(tmpSboxAffMatrix)) << endl;
	cout << "DecInv: " << endl << (inv(tmpSboxAffMatrixDec)) << endl;
	cout << "Dec constant: " << endl <<  (tmpSboxAffMatrixDec * tmpSboxAffConst) << endl;
#endif

	// Computing whole Sboxes with inversion + affine transformation in generic AES
	for(i=0; i<AES_FIELD_SIZE; i++){
		long tmpLong;

		// i is now long representation, gInv transforms it to exponent power to obtain inverse.
		// Also getLong(g[gInv[i]]) == i
		GF2E transValue = i==0 ? GF2E::zero() : g[255-gInv[i]];

		// DEBUG
		//cout << "1 Inversion check["<<i<<"]"<<rep(transValue)<<": [" << g[gInv[i]] << " x " << transValue << "] = " << (transValue * g[gInv[i]]) << endl;

		mat_GF2 resMatrix = tmpSboxAffMatrix * colVector(rep(transValue), AES_FIELD_DIM) + tmpSboxAffConst;
		colVector(transValue, resMatrix, 0);
		tmpLong = getLong(transValue);
		this->sboxAffineGF2E[i] = transValue;
		this->sboxAffine[i] = tmpLong;
		this->sboxAffineInvGF2E[tmpLong] = GF2EFromLong(i, AES_FIELD_DIM);
		this->sboxAffineInv[tmpLong] = i;

		// Inversion, idea is the same, i is long representation of element in GF, apply inverted affine transformation and take inverse
		// Ax^{-1} + c is input to this transformation
		// [A^{-1} * (A{x^-1} + c) + d]^{-1} is this transformation;
		// correctness: [A^{-1} * (Ax^-1 + c) + d]^{-1} =
		//				[A^{-1}Ax^{-1} + A^{-1}c + d]^{-1} =	//	A^{-1}c = d
		//				[x^{-1}        + 0]^{-1} =
		//				x
		//
		// Computation is useless, we have inversion of transformation right from transformation above
		// by simply swapping indexes. This is just for validation purposes to show, that it really works and how
		int ii = tmpLong;
		GF2E transValue2 = transValue;
		resMatrix = tmpSboxAffMatrixDec * colVector(rep(transValue2), AES_FIELD_DIM) + tmpSboxAffConstDec;
		colVector(transValue2, resMatrix, 0);
		tmpLong = getLong(transValue2);

		transValue2 = tmpLong==0 ? GF2E::zero() : g[255-gInv[tmpLong]];
		tmpLong = getLong(transValue2);
		if (this->sboxAffineInv[ii] != getLong(transValue2)){
			cout << "!!Integrity problem in Sbox value: " << CHEX(ii) << endl;
		}
	}

	// 6. MixColumn operations
	// modulus x^4 + 1
	mixColModulus.SetLength(5);
	mixColModulus[0] = g[0];
	mixColModulus[4] = g[0];

	// 03 x^3 + 01 x^2 + 01 x + 02
	mixColMultiply.SetLength(4);
	mixColMultiply[0] = g[25];
	mixColMultiply[1] = g[0];
	mixColMultiply[2] = g[0];
	mixColMultiply[3] = g[1];

	// inverse polynomial
	// TODO:::: !check!!!!
	mixColMultiplyInv.SetLength(4);
	mixColMultiplyInv[0] = g[223];
	mixColMultiplyInv[1] = g[199];
	mixColMultiplyInv[2] = g[238];
	mixColMultiplyInv[3] = g[104];

	// MixCols multiplication matrix based on mult polynomial -  see Rijndael description of this.
	// Polynomials have coefficients in GF(256).
	mixColMat.SetDims(4,4);
	mixColInvMat.SetDims(4,4);
	for(i=0; i<4; i++){
		for(c=0; c<4; c++){
			mixColMat.put(i, c, mixColMultiply[(i+4-c) % 4]);
			mixColInvMat.put(i, c, mixColMultiplyInv[(i+4-c) % 4]);
		}
	}

	// Round key constant RC (for key schedule) obeys this reccurence:
	// RC[0] = 1
	// RC[i] = '02' * RC[i-1] = x * RC[i-1] = x^{i-1} `mod` R(X)
	RC[0] = g[0];
	for(i=1; i<16; i++){
		RC[i] =  g[25] * RC[i-1];
	}

	return;
}

// input = state array
void GenericAES::encryptInternal(mat_GF2E& state, vec_GF2E& expandedKey){
	// Now strictly assume that key is 128bit long, thus perform 10 rounds of encryption
	//
	// Cipher explanation:
	// round
	// 		ByteSub, ShiftRows, MixColumn, AddRoundKey
	// final round:
	//		ByteSub, ShiftRows,     --   , AddRoundKey
	// cipher:
	// 		AddRoundKey(State, ExpandedKey)
	//		for(i=1; i<Nr; i++) Round(State, ExpandedKey + Nb*i)
	//		finalRound(State, ExpandedKey + Nb*Nr)
	int r;
	restoreModulus();

	// Add Round key:
	this->AddRoundKey(state, expandedKey, 0);
	// rounds
	for(r=1; r<=10; r++){
		this->ByteSub(state);
		this->ShiftRows(state);
		if(r<10) this->MixColumn(state);
		this->AddRoundKey(state, expandedKey, 16*r);
	}
}

// input = state array
void GenericAES::decryptInternal(mat_GF2E& state, vec_GF2E& expandedKey){
	// Now strictly assume that key is 128bit long, thus perform 10 rounds of encryption
	//
	// Cipher explanation:
	// round
	// 		ShiftRows, ByteSub, AddRoundKey, MixColumns
	// final round:
	//		ShiftRows, ByteSub, AddRoundKey, --
	// cipher:
	// 		AddRoundKey(State, ExpandedKey + Nb*Nr)
	//		for(i=Nr-1; i>0; i++) Round(State, ExpandedKey + Nb*i)
	//		finalRound(State, ExpandedKey + 0)
	int r;
	restoreModulus();

	// Add Round key:
	this->AddRoundKey(state, expandedKey, 16*10);
	// rounds
	for(r=9; r>=0; r--){
		this->ShiftRowsInv(state);
		this->ByteSubInv(state);
		this->AddRoundKey(state, expandedKey, 16*r);
		if(r!=0) this->MixColumnInv(state);
	}
}

// test routine - verify inversion
int GenericAES::testByteSub(){
	int repr = 0;
	for(repr=0; repr < AES_FIELD_SIZE; repr++){
		if (sboxAffineInv[sboxAffine[repr]] != repr) return -1;
	}
	return 0;
}

// test routine - verify inversion
int GenericAES::testMixColumn(){
	return 1;
}

void GenericAES::printAll() {
	int i=0,c=0;
	restoreModulus();

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

	cout << "========================================================" << endl;
	cout << "Powers of generator g^i: " << endl;
	dumpVector(g, 256);

	cout << "Sbox table: " << endl;
	dumpVector(sboxAffine, 256);

	cout << "Sbox inverse table: " << endl;
	dumpVector(sboxAffineInv, 256);
	cout << "Testing inverse in Sboxes " << this->testByteSub() << endl;

	cout << "Base change matrix: " << endl;
	dumpMatrix(T);

	cout << "Base change ivn matrix: " << endl;
	dumpMatrix(Tinv);

	cout << "MixCol mult " << mixColMultiply << endl;
	cout << "MixCol multInv " << mixColMultiplyInv << endl;
	cout << "MixCol mod " << mixColModulus << endl;

	cout << "Mix col:   " << endl;
	dumpMatrix(mixColMat);

	cout << "MixInccol: " << endl;
	dumpMatrix(mixColInvMat);

	cout << "MixColInverse test: " << endl;
	mat_GF2E tmpMat = mixColInvMat * mixColMat;
	dumpMatrix(tmpMat);

	cout << "Round constants: " << endl;
	dumpVector(RC, 16);

	// try round key expansion
	vec_GF2E roundKey;
	vec_GF2E key;

	key.SetLength(128);
	expandKey(roundKey, key, KEY_SIZE_16);
	cout << "Round key for ZERO key for 16B: " << endl;
	dumpVector(roundKey);

	//key.SetLength(256);
	//expandKey(roundKey, key, KEY_SIZE_32);
	//cout << "Round key for ZERO key for 32B: " << endl;
	//dumpVector(roundKey);

	mat_GF2E state(INIT_SIZE, 4, 4);
	cout << "Plaintext: " << endl;
	dumpMatrix(state);

	cout << "Testing encryption: " << endl;
	encryptInternal(state, roundKey);
	dumpMatrix(state);

	cout << "Testing backward decryption: " << endl;
	decryptInternal(state, roundKey);
	dumpMatrix(state);

	return;
}




