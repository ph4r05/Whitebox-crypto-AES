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

/**
 * Generates whole AES table with generator
 */
void GenericAES::build() {
	// 1. generate whole field with given generator
	int i=0;
	generator.LoopHole().SetLength(8);
	generator.LoopHole().SetMaxLength(8);

	// current element of the field
	GF2E cur;
	cur.LoopHole().SetMaxLength(8);
	cur.LoopHole().SetLength(8);
	//cur.init(this->modulus);

	cur=1;
	for(i=0; i<AES_FIELD_SIZE; i++){
		this->g[i] = cur;
		this->gInv[getLong(cur)] = i;
		cur = cur * this->generator;

		this->g[i].LoopHole().SetLength(8);
		this->g[i].LoopHole().SetMaxLength(8);
		this->g[i].LoopHole().normalize();
	}

	// 2. compute inverses in terms of generator exponent
	this->sbox[0] = -1;
	for(i=1; i<AES_FIELD_SIZE; i++){
		this->sbox[i] = 255-i;
	}

	return;
}

void GenericAES::printAll() {
	int i=0;
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

	return;
}
