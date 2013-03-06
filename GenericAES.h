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
	void keycore(unsigned char *word, int iteration);
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
    GF2E RC[20];

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
};

#endif /* GENERICAES_H_ */
