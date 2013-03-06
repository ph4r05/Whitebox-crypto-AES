/*
 * GenericAES.h
 *
 *  Created on: Mar 1, 2013
 *      Author: ph4r05
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
#include "NTLUtils.h"
#define AES_FIELD_SIZE 256
#define AES_FIELD_DIM 8

// number of columns of state array ; = BLOCKS / 32, static for AES
#define AES_NB 4
// number of columns of key array
#define AES_NK_128 4
#define AES_NK_256 8

// AES KEY size
enum keySize {
	KEY_SIZE_16 = 16,
	KEY_SIZE_24 = 24,
	KEY_SIZE_32 = 32
};

NTL_CLIENT
class GenericAES {
public:
	GenericAES();
	virtual ~GenericAES();
	
	const GF2X& getModulus() { return modulus; }
	const GF2E& getGenerator() { return generator; }
    
	void setModulus(GF2X aMod){ modulus = aMod; }
	void setGenerator(GF2E aGen){ generator = aGen; }
	void build();
	void printAll();

	static mat_GF2 getDefaultAffineMatrix (void);
	static vec_GF2 getDefaultAffineConst (void);
	static mat_GF2 getDefaultAffineMatrixDec (void);
	static vec_GF2 getDefaultAffineConstDec (void);

	inline static int getNumberOfRounds(enum keySize keySize){ return keySize/4+6; };
	void expandKey(vec_GF2E& expandedKey, vec_GF2E& key, enum keySize size);

	void encryptInternal(mat_GF2E& result, vec_GF2E& expandedKey);
	void decryptInternal(mat_GF2E& result, vec_GF2E& expandedKey);

	inline int mod4(int a){
		int c = a % 4; return c<0 ? c+4 : c;
	}

 	inline void ByteSub(mat_GF2E& state){
		int i,j;
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
				state[i][j] = this->sboxAffineGF2E[getLong(state[i][j])];
			}
		}
	}

	inline void ByteSubInv(mat_GF2E& state){
		int i,j;
		for(i=0;i<4;i++){
			for(j=0;j<4;j++){
				state[i][j] = this->sboxAffineInvGF2E[getLong(state[i][j])];
			}
		}
	}

	inline void AddRoundKey(mat_GF2E& state, vec_GF2E& expandedKey, unsigned int offset){
		int i,j;
		for(i=0; i<4; i++){
			for(j=0; j<4; j++){
				state[i][j] += expandedKey[offset + i*4+j];
			}
		}
	}

	inline void ShiftRows(mat_GF2E& state){
		// 1. row = no shift. 2. row = cyclic shift to the left by 1
		// for AES with Nb=4, left shift for rows are: 1=1, 2=2, 3=3.
		GF2E tmp;
		int i,j=0;
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

	int testByteSub();
	int testMixColumn();

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
     * Base change address
     */
    mat_GF2 T;
    mat_GF2 Tinv;

    /**
     * Mix column modulus polynomial & multiply polynomial
     */
    GF2EX mixColModulus;
    GF2EX mixColMultiply;
    GF2EX mixColMultiplyInv;

    /**
     * Mix column GF2E matrix computed
     */
    mat_GF2E mixColMat;
    mat_GF2E mixColInvMat;

    /**
     * Round key constant
     */
    GF2E RC[16];

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
     
     /**
      * 4. table: Sbox + affine transformation
      * long representation -> long representation
      */
     long sboxAffine[AES_FIELD_SIZE];
     long sboxAffineInv[AES_FIELD_SIZE];
     GF2E sboxAffineGF2E[AES_FIELD_SIZE];		// cached object version for easy manipulation
     GF2E sboxAffineInvGF2E[AES_FIELD_SIZE];	// cached object version for easy manipulation
};

#endif /* GENERICAES_H_ */
