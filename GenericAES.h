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
#include <NTL/vec_long.h>
#include <NTL/new.h>
#include <math.h>
#include "NTLUtils.h"
#define AES_FIELD_SIZE 256

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
private:
    /**
     * Basic definition parameters
     */
    GF2X modulus;
    GF2E generator;
    
    /**
     * Mix column modulus polynomial & multiply polynomial
     */
    GF2EX mixColModulus;
    GF2EX mixColMultiply;    
  
    /**
     * Build lookup tables here from generator.
     * AES fiels has 256 elements, so we use simple static arrays
     * 
     * 1. lookup table: generator exponents. g[0]=1, g[1]=g, g[2]=g^2, ... 
     */
    GF2E g[AES_FIELD_SIZE];
    
    /**
     * 2. lookup table: given element from field (represented as long) to 
     * generator exponent mapping
     */
     long gInv[AES_FIELD_SIZE];
     
     /**
      * 3. table: Sbox (simple inversion)
      * generator exponent to generator exponent
      */
     long sbox[AES_FIELD_SIZE];
     
     /**
      * 4. table: Sbox + affine transformation
      * exponent -> exponent
      */
      long sboxAffine[AES_FIELD_SIZE];  
};

#endif /* GENERICAES_H_ */
