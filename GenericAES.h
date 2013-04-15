/*
 * GenericAES.h
 *
 *  Created on: Mar 1, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#ifndef GENERICAES_H_
#define GENERICAES_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

// NTL dependencies
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>

#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_GF2E.h>
#include <NTL/vec_long.h>
#include <NTL/new.h>
#include <math.h>
#include <assert.h>
#include "base.h"
#include "NTLUtils.h"
#define AES_FIELD_SIZE 256
#define AES_FIELD_DIM 8

// there is 30 irreducible polynomials of degree 8 in GF2[x] ring
#define AES_IRRED_POLYNOMIALS 30
#define AES_GENERATORS 8

// number of columns of state array ; = BLOCKS / 32, static for AES
#define AES_NB 4
// number of columns of key array
#define AES_NK_128 4
#define AES_NK_256 8
#define AES_TESTVECTORS 4
#define AES_BYTES 16

// AES KEY size
enum keySize {
	KEY_SIZE_16 = 16,
	KEY_SIZE_24 = 24,
	KEY_SIZE_32 = 32
};

NTL_CLIENT
class GenericAES {
public:
	/**
	 * List of all irreducible polynomials in GF2[x] of degree 8
	 */
	static long irreduciblePolynomials[AES_IRRED_POLYNOMIALS];

	/**
	 * List of all generators usable to generate AES FIELD for given irreducible polynomial (1st index)
	 */
	static long generators[AES_IRRED_POLYNOMIALS][AES_GENERATORS];

	/**
	 * AES testvectors
	 */
	static unsigned char testVect128_key[AES_BYTES];
	static unsigned char testVect128_plain[5][AES_BYTES];
	static unsigned char testVect128_cipher[5][AES_BYTES];

	static unsigned char testVect256_key[32];
	static unsigned char testVect256_plain[AES_TESTVECTORS][AES_BYTES];
	static unsigned char testVect256_cipher[AES_TESTVECTORS][AES_BYTES];

	GenericAES();
	virtual ~GenericAES();
	
	const GF2X& getModulus() { return modulus; }
	const GF2E& getGenerator() { return generator; }
    
	/**
	 * Initializes AES with given indices to modulus and generator array
	 */
	void initFromIndex(int modulus, int generator);

	/**
	 * Initializes AES with given modulus and generator in long representation
	 */
	void init(long modulus, long generator);


	void setModulus(GF2X aMod){ modulus = aMod; }
	void setGenerator(GF2E aGen){ generator = aGen; }
	inline void restoreModulus(){ modulusContext.restore(); }

	/**
	 * Main method that prepares whole object for use.
	 * Generates Sboxes, MixColumn matrices, T, Tinv matrices.
	 *
	 * Modulus and generator has to be set.
	 * Calling this method restores modulus in GF2E class that was used
	 * to construct this AES.
	 */
	void build();

	/**
	 * Debug info - dump Sboxes, MixColumns and performs basic tests
	 */
	void printAll();

	static mat_GF2 getDefaultAffineMatrix (void);
	static vec_GF2 getDefaultAffineConst (void);
	static mat_GF2 getDefaultAffineMatrixDec (void);
	static vec_GF2 getDefaultAffineConstDec (void);

	/**
	 * Number of rounds of AES depends on key size
	 */
	inline static int getNumberOfRounds(enum keySize keySize){ return keySize/4+6; };

	/**
	 * Expands encryption key <key> to round key <expandedKey> according to AES specification.
	 * This round key is used during encryption in AddRoundKey() operation
	 */
	void expandKey(vec_GF2E& expandedKey, vec_GF2E& key, enum keySize size);

	void encryptInternal(mat_GF2E& result, vec_GF2E& expandedKey);
	void decryptInternal(mat_GF2E& result, vec_GF2E& expandedKey);

	void encrypt(vec_GF2E& result, vec_GF2E& expandedKey);
	void decrypt(vec_GF2E& result, vec_GF2E& expandedKey);

	void applyT(mat_GF2E& state);
	void applyT(vec_GF2E& state);
	void applyT(GF2E& state);
	void applyTinv(mat_GF2E& state);
	void applyTinv(vec_GF2E& state);
	void applyTinv(GF2E& state);

	inline int mod4(int a){ int c = a % 4; return c<0 ? c+4 : c; }

	inline GF2E& ByteSub(GF2E& e){
		this->restoreModulus();
		return this->sboxAffineGF2E[getLong(e)];
	}

	inline long ByteSub(long e){
		this->restoreModulus();
		return this->sboxAffine[e];
	}

 	inline void ByteSub(mat_GF2E& state){
		int i,j;
		this->restoreModulus();
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
				state[i][j] = this->sboxAffineGF2E[getLong(state[i][j])];
			}
		}
	}

	inline GF2E& ByteSubInv(GF2E& e){
		this->restoreModulus();
		return this->sboxAffineInvGF2E[getLong(e)];
	}

	inline long ByteSubInv(long e){
		this->restoreModulus();
		return this->sboxAffineInv[e];
	}

	inline void ByteSubInv(mat_GF2E& state){
		int i,j;
		this->restoreModulus();
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
				state[i][j] = this->sboxAffineInvGF2E[getLong(state[i][j])];
			}
		}
	}

	inline void AddRoundKey(mat_GF2E& state, vec_GF2E& expandedKey, unsigned int offset){
		int i,j;
		this->restoreModulus();
		for(i=0; i<4; i++){
			for(j=0; j<4; j++){
				state[i][j] += expandedKey[offset + j*4+i];
			}
		}
	}

	inline void ShiftRows(mat_GF2E& state){
		// 1. row = no shift. 2. row = cyclic shift to the left by 1
		// for AES with Nb=4, left shift for rows are: 1=1, 2=2, 3=3.
		GF2E tmp;
		int i,j=0;
		this->restoreModulus();
		for(i=1; i<4;i++){
			for(j=1; j<=i; j++){
				tmp = state[i][0];
				state[i][0] = state[i][1];
				state[i][1] = state[i][2];
				state[i][2] = state[i][3];
				state[i][3] = tmp;
			}
		}
	}

	inline void ShiftRowsInv(mat_GF2E& state){
		// 1. row = no shift. 2. row = cyclic shift to the left by 1
		// for AES with Nb=4, left shift for rows are: 1=1, 2=2, 3=3.
		GF2E tmp;
		signed int i=0,j=0;
		this->restoreModulus();
		for(i=1; i<4;i++){
			for(j=1; j<=i; j++){
				tmp = state[i][3];
				state[i][3] = state[i][2];
				state[i][2] = state[i][1];
				state[i][1] = state[i][0];
				state[i][0] = tmp;
			}
		}
	}

	inline void MixColumn(mat_GF2E& state){
		int i,j;
		this->restoreModulus();
		mat_GF2E tmpMat(INIT_SIZE, 4,1);
		mat_GF2E resMat(INIT_SIZE, 4,1);

		for(i=0; i<4; i++){
			// copy i-th column to 4*1 matrix - for multiplication
			for(j=0; j<4; j++){
				tmpMat.put(j,0, state.get(j,i));
			}

			resMat = this->mixColMat * tmpMat;

			// copy result back to i-th column
			for(j=0; j<4; j++){
				state.put(j,i, resMat.get(j,0));
			}
		}
	}

	inline void MixColumnInv(mat_GF2E& state){
		int i,j;
		this->restoreModulus();
		mat_GF2E tmpMat(INIT_SIZE, 4,1);
		mat_GF2E resMat(INIT_SIZE, 4,1);

		for(i=0; i<4; i++){
			// copy i-th column to 4*1 matrix - for multiplication
			for(j=0; j<4; j++){
				tmpMat.put(j,0, state.get(j,i));
			}

			resMat = this->mixColInvMat * tmpMat;
			// copy result back to i-th column
			for(j=0; j<4; j++){
				state.put(j,i, resMat.get(j,0));
			}
		}
	}

	int testWithVectors();
	int testByteSub();
	int testMixColumn();

	/**
	 * Generates A1 and A2 affine transformations for S-box self affine equivalence.
	 * Result is returned as lookup table.
	 *
	 * @param a = affine constant term [0..255]
	 * @param q = power matrix exponent [1..8]
	 *
	 * Creates A1 affine transformation: A1(x) = [a] . Q^q . x
	 *  Parameters:
	 *       a     \in a \in GF(2^8) \ {0}
	 *       q     \in [1,8]
	 *       field     = AES field generated by generator
	 *       generator = generator of GF(2^8) field
	 *       modulo    = irreducible polynomial defining elements of field GF(2^8)
	 *
	 * Creates A2 affine transformation: A2(x) = A(Q^(8-q) . [a] . A^(-1)(x))
     *  Parameters:
     *       a     \in a \in GF(2^8) \ {0}
     *       q     \in [1,8]
     *       field     = AES field generated by generator
     *       generator = generator of GF(2^8) field
     *       modulo    = irreducible polynomial defining elements of field GF(2^8)
	 */
	void generateA1A2Relations(vec_GF2E& A1, vec_GF2E& A2, int a, int q);
	/**
	 * Same as previous, but generates a,q by random from interval
	 */
	void generateA1A2Relations(vec_GF2E& A1, vec_GF2E& A2);

	/**
	 * Generates 8x8 bit matrix representing multiplication by a \in GF(2^8) by
     * modulo <modulo> polynomial. Element of matrix is in GF(2). Columns are
     * transformed bases = i-th column = e_i * a (mod <modulo>), where e_i=x^i.
     *
     * Since multiplication by constant in G is linear operation, there exists
     * matrix [a] representing this multiplication.
     *
     * Implementation: perform this transformation on polynomial base of space.
     * See Theorem 13, page 92, prednasky_MB101_podzim2008.pdf to verify linearity
     * and expressing linear functions in vector spaces as matrices.
	 */
	mat_GF2 makeMultAMatrix(int a);

	/**
	 *  Computes matrix performing squaring in GF(2^8).
	 *  Squaring is linear operation in this field: (a+b)^2 = a^2 + b^2
	 *
	 *  @see [In How Many Ways Can You Write Rijndael? (2002) by Elad Barkan , Eli Biham]
	 */
	mat_GF2 makeSquareMatrix(int q);

	/**
	 * Simple test of A1 and A2 relations.
	 * This equations hould hold for all x: (A2 * S * A1) (x) = S (x)
	 *
	 * @returns 0 on success, -num in case of num of failed mappings
	 */
	int testA1A2Relations(vec_GF2E& A1, vec_GF2E& A2, bool encryption=true, bool output=true);

	/**
	 * Tests linearity of A1 relation: A1(x) ^ A1(y) == A1(x^y)
	 */
	int testA1XorLinearity(vec_GF2E& A1);

    /**
     * Mix column GF2E matrix computed
     */
    mat_GF2E mixColMat;
    mat_GF2E mixColInvMat;

    /**
     * 4. table: Sbox + affine transformation
     * long representation -> long representation
     */
    long sboxAffine[AES_FIELD_SIZE];
    long sboxAffineInv[AES_FIELD_SIZE];
    GF2E sboxAffineGF2E[AES_FIELD_SIZE];		// cached object version for easy manipulation
    GF2E sboxAffineInvGF2E[AES_FIELD_SIZE];	// cached object version for easy manipulation

    /**
     * Base change matrices
     */
    mat_GF2 T;
    mat_GF2 Tinv;

    /**
     * Round key constant
     */
    GF2E RC[16];
private:

    /**
     * Basic definition parameters
     */
    GF2X modulus;
    GF2E generator;
    
    /**
     * Modulus context
     */
    GF2EContext modulusContext;

    /**
     * Previous context
     */
    GF2EContext prevModulusContext;

    /**
     * Mix column modulus polynomial & multiply polynomial
     */
    GF2EX mixColModulus;
    GF2EX mixColMultiply;
    GF2EX mixColMultiplyInv;

    /**
     * AES field.
     *
     * Build lookup tables here from generator.
     * AES field has 256 elements, so we use simple static arrays
     * 
     * 1. lookup table: generator exponents. g[0]=1, g[1]=g, g[2]=g^2, ...
     * generator exponent -> field element in GF2E
     */
    GF2E g[AES_FIELD_SIZE+1];
    
    /**
     * 2. lookup table: given element from field (represented as long) to 
     * generator exponent mapping
     *
     * long representation -> generator exponent
     */
     long gInv[AES_FIELD_SIZE+1];
     
     /**
      * 3. table: Sbox (simple field element inversion)
      * generator exponent to generator exponent
      */
     long sbox[AES_FIELD_SIZE];
};

#endif /* GENERICAES_H_ */
