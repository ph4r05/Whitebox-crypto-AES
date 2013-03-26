//============================================================================
// Name        : MGR_NTL.cpp
// Author      : Dusan Klinec (ph4r05)
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
// *  Author: Dusan Klinec (ph4r05)
// *
// *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#define WBAES_BOOTS_SERIALIZATION 1

// NTL dependencies
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_long.h>
#include <math.h>
#include <vector>
#include "GenericAES.h"
#include "NTLUtils.h"
#include "MixingBijections.h"
#include "WBAES.h"
#include "WBAESGenerator.h"
#include "md5.h"
#include <iostream>
#include <fstream>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
#include <set>
#include <map>
#include "LinearAffineEq.h"


NTL_CLIENT
using namespace std;
using namespace NTL;
using namespace boost;
using namespace wbacr;
using namespace wbacr::laeqv;

typedef struct _a1a2rec {
	unsigned long int id;
	int count;
	vector<unsigned long int> * vec;
} a1a2rec;
typedef unordered_map<std::string, a1a2rec> a1a2map;

#define GENERIC_AES_DEBUG 1

std::string hashLookupTable(vec_GF2E s);
int MBgen(void);
int dualAESTest(void);
int A1A2relationsGenerator(void);

// hashes lookup table using MD5 hash
std::string hashLookupTable(vec_GF2E s){
	std::string inputBuffer = dumpVector2str(s);
	char buff[inputBuffer.size()];					// c++0x feature
	inputBuffer.copy(buff, inputBuffer.size(), 0);	// copy serialized lookup table to char array

	MD5_CTX mdContext;
	MD5Init(&mdContext);
	MD5Update(&mdContext, (unsigned char * )buff, inputBuffer.size());
	MD5Final(&mdContext);

	// hexcoding
	ostringstream ss;
	int i=0; for(i=0; i<MD5_DIGEST_LENGTH; i++) ss << CHEX(mdContext.digest[i]);
	inputBuffer = ss.str();
	return inputBuffer;
}

// hardcoded elements
// http://stackoverflow.com/questions/2236197/c-easiest-way-to-initialize-an-stl-vector-with-hardcoded-elements


int main(void) {
	long i,j;

	// very poor PRNG seeding, but just for now
	srand((unsigned)time(0));
	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	//A1A2relationsGenerator();
	//exit(2);

	//dualAESTest();
	//exit(2);

	GenericAES defAES;
	defAES.init(0x11B, 0x03);
	defAES.printAll();
	WBAESGenerator generator;
	WBAES genAES;


	//generator.useDualAESIdentity=true;
	//generator.useIO04x04Identity=true;
	//generator.useIO08x08Identity=true;
	//generator.useMB08x08Identity=true;
	///generator.useMB32x32Identity=true;
	//int errors = generator.testWithVectors(true, genAES);
	//cout << "Testing done, errors: " << errors << endl;

/*
	cout << "===Done===" << endl << "Going to generate WBAES..." << endl;
	WBAESGenerator generator;
	WBAES genAES;
	BYTE aesKey[16]; memset(aesKey, 0, sizeof(BYTE)*16);
	CODING8X8_TABLE coding[16];
 	generator.generateIO128Coding(coding);
 	generator.generateTables(aesKey, KEY_SIZE_16, genAES, coding, true);

 	cout << "WBAES generated! ; size: " << sizeof(genAES) << endl;
 	cout << "Going to encrypt..." << endl;
 	W128b newState; memset(&newState, 0, sizeof(W128b));
 	genAES.encrypt(newState);

 	cout << "Encrypted! " << endl;
 	dumpW128b(newState);

 	cout << "Going to generate decryption tables" << endl;
 	generator.generateTables(aesKey, KEY_SIZE_16, genAES, coding, false);
 	cout << "Generated; goind to decrypt back..." << endl;
 	genAES.decrypt(newState);

 	cout << "Encrypted! " << endl;
 	dumpW128b(newState);

 	cout << endl << "Exiting..." << endl;
*/


	//
	// Linear equivalence algorithm
	//
	bset Ua, Ub, Na, Nb, Ca, Cb;
	bset Ua_back, Ub_back, Ca_back, Cb_back;	// backup in case of incorrect guess
	bsetElem S2[256];
	bsetElem S2inv[256];
	bsetElem * S1 = S2;
	bsetElem * S1inv = S2inv;
	bsetElem guesses1[255]; int guess1Idx=0;
	randomPermutationT(guesses1, 255, 1);
	for(i=0; i<256; i++) {
		S2[i] = defAES.sboxAffine[i];
		S2inv[i] = defAES.sboxAffineInv[i];
	}

	// mappings known
	smap mapA, mapB, tmpMapA, tmpMapB;

	// init Ua, Ub, Na, Nb
	Ca.insert(0); Na.insert(0);
	Cb.insert(0); Nb.insert(0);
	tmpMapA.insert( pair<bsetElem,bsetElem>(0, 0) );
	mapA.insert(    pair<bsetElem,bsetElem>(0, 0) );
	tmpMapB.insert( pair<bsetElem,bsetElem>(0, 0) );
	mapB.insert(    pair<bsetElem,bsetElem>(0, 0) );

	// random initialization of Ua, Ub
	bsetElem rndInit[255];
	randomPermutationT(rndInit, 255, 1);
	for(i=0; i<255; i++){
		Ua.insert(rndInit[i]);
		Ub.insert(rndInit[i]);
	}

	bset::const_iterator it1, it2, it3;
	while(Ua.empty()==false && Ub.empty()==false && guess1Idx < 255){
		cout << endl << "====================================================================================" << endl
				<< "Main cycle started " << endl;

		//
		// starting with new guess
		//
		if (Na.empty() && Nb.empty()){
			bsetElem x, guess;
			int rnd = rand() % Ua.size();
			//
			// At first, backup Ca, Cb, Ua, Ub - will be restored in case of incorrect guess
			//
			Ua_back = Ua;   Ub_back = Ub;
			Ca_back = Ca;   Cb_back = Cb;
			mapA = tmpMapA; mapB = tmpMapB;

			//
			// 1. If previous guess rejected, restore Ca, Cb, Ua, Ub
			// Guess A(x) for some x \in Ua
			// Set Na = {x}, Ua = Ua / {x}
			//
			cout << "New guess;" << endl;
			it1 = Ua.begin(); for(i=0; i<rnd; ++i, ++it1);
			x = *it1; Ua.erase(it1);
			Na.insert(x);
			cout << "G: x was selected as [" << x << "]" << endl;

			guess = guesses1[guess1Idx++];
			cout << "G: guessed value: " << guess << endl;
			tmpMapA.insert( smapElem(x, guess) );
		}

		//
		// Na cycle
		//
		while(Na.empty() == false){
			bsetElem x, y;
			cout << endl << "A: Cycle 1 start, |Na| = " << Na.size() << endl;

			// Pick x \in Na; Na = Na \ {x};
			it1 = Na.begin(); x = *it1; Na.erase(it1);
			cout << "A: newX is [" << x << "]" << endl;

			// Nb = S2( x + Ca ) \ Cb
			if (Ca.size() > 0){
				bset tmpSet;
				for (it1=Ca.begin(); it1!=Ca.end(); ++it1){
					bsetElem curr = *it1;
					bsetElem tmp = x ^ curr;
					tmpSet.insert(S2[tmp]);
					// Use linearity of A to build mapping for tmp
					tmpMapA.insert(smapElem(tmp, tmpMapA[x] ^ tmpMapA[curr]));

					cout << "    curr=" << setw(2) << (curr)
							<< "; tmp = " << setw(2) << tmp
							<< "; S2[tmp]    = " << setw(2) << S2[tmp] << endl;
					cout << "A:     adding mapping for A(x^curr) = A(" << (tmp) << ") = A(x) ^ A(curr) = " << tmpMapA[x] << " ^ " << tmpMapA[curr] << " = " << (tmpMapA[x] ^ tmpMapA[curr]) << endl;

					//
					// A(x) = S^{-1}_1 (B(S_2(x)))
					// S1 * A = B * S2
					if (Cb.count(S2[tmp])==0){
						const bsetElem amap = tmpMapA[tmp];

						cout << "A:     adding mapping for B(S2(tmp)) = B("<< S2[tmp] <<") = S1(A(tmp)) = S1(" << amap << ") = " << S1[amap] << endl;
						tmpMapB.insert( pair<bsetElem,bsetElem>(S2[tmp], S1[amap]));
					}
				}
				Nb = LinearAffineEq::setDiff(tmpSet, Cb);
			}

			// Ca = Ca U (x + Ca)
			cout << "A: Ca size pre: " << Ca.size() << endl;
			bset tmpSet(Ca);
			for(it1=tmpSet.begin(); it1 != tmpSet.end(); ++it1){
				Ca.insert( x ^ (*it1) );
			}
			Ca.insert(x);
			cout << "A: Ca size post: " << Ca.size() << endl;

			double vectKnown = Nb.size() + log2(Cb.size());
			cout << "A: vect knownB: " << vectKnown << endl;
			if (vectKnown>8){


			}
		}
		cout << endl;

		//
		// Nb cycle
		//
		while(Nb.empty() == false){
			bsetElem x, y;
			cout << "B: Cycle 2 start, |Nb| = " << Nb.size() << endl;

			it1 = Nb.begin(); x = *it1; Nb.erase(it1);
			cout << "B: newX is [" << x << "]" << endl;

			if (Cb.size() > 0){
				bset tmpSet;
				for (it1=Cb.begin(); it1!=Cb.end(); ++it1){
					bsetElem curr = *it1;
					bsetElem tmp = x ^ curr;
					tmpSet.insert(S2inv[tmp]);

					// Use linearity of B to build mapping for tmp
					tmpMapB.insert(smapElem(tmp, tmpMapB[x] ^ tmpMapB[curr]));

					cout << "    curr=" << setw(2) << (curr)
							<< "; tmp = " << setw(2) << tmp
							<< "; S2inv[tmp] = " << setw(2) << S2inv[tmp] << endl;
					cout << "B:     adding mapping for B(x^curr) = B(" << (tmp) << ") = B(x) ^ B(curr) = " << tmpMapB[x] << " ^ " << tmpMapB[curr] << " = " << (tmpMapB[x] ^ tmpMapB[curr]) << endl;

					//
					//       A            = S^{-1}_1 * B * S_2
					//       A * S_2^{-1} = S^{-1}_1 * B
					if (Ca.count(S2inv[tmp])==0){
						const bsetElem bmap = tmpMapB[tmp];

						cout << "B:     adding mapping for A(S2inv(tmp)) = A(" << S2inv[tmp] << ") = S1inv(B(tmp)) = S1inv(" << bmap << ") = " << (S1inv[bmap]) << endl;
						tmpMapA.insert(smapElem(S2inv[tmp], S1inv[bmap]));
					}
				}
				Na = LinearAffineEq::setDiff(tmpSet, Ca);
			}

			// Cb = Cb U (x + Cb)
			cout << "B: Cb size pre: " << Cb.size() << endl;
			bset tmpSet(Cb);
			for(it1=tmpSet.begin(); it1 != tmpSet.end(); ++it1){
				Cb.insert( x ^ *it1 );
			}
			Cb.insert(x);
			cout << "B: Cb size post: " << Cb.size() << endl;

			double vectKnown = Na.size() + log2(Ca.size());
			cout << "B: vect knownA: " << vectKnown << endl;
			if (vectKnown>8){
				cout << "B: ## check linearity of A, derive,..." << endl;
				LinearAffineEq::dumpMap(tmpMapA);


			}
		}

		cout << endl << "EEEnd of both cycles, remove Ca, Cb from Ua Ub" << endl;
		Ua = LinearAffineEq::setDiff(Ua, Ca);
		Ub = LinearAffineEq::setDiff(Ub, Cb);

		cout    << " |Ua| = " << setw(3) << Ua.size()
				<< " |Ub| = " << setw(3) << Ub.size()
				<< " |Ca| = " << setw(3) << Ca.size()
				<< " |Cb| = " << setw(3) << Cb.size()
				<< " |Na| = " << setw(3) << Na.size()
				<< " |Nb| = " << setw(3) << Nb.size() << endl;

		cout << "Dump mapA: " << endl;
		LinearAffineEq::dumpMap(tmpMapA);

		cout << "Dump mapB: " << endl;
		LinearAffineEq::dumpMap(tmpMapB);
	}
}

int A1A2relationsGenerator(void){
	// very poor PRNG seeding, but just for now
	srand((unsigned)time(0));
	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	// unordered_map
	a1a2map a1store;
	a1a2map a2store;

	ofstream dump;
	ofstream dumpA;
	dump.open("/media/share/AES_A1A2dump.txt");
	dumpA.open("/media/share/AES_signature.txt");


	dump  << "id;polynomial;generator;qq;ii;problems;A1;A2" << endl;
	dumpA << "polynomial;generator;sbox;sboxinv;mixcol;mixcolinv" << endl;

	GenericAES defAES;
	defAES.init(0x11B, 0x03);
	defAES.printAll();

	int AES_gen, AES_poly;
	for(AES_poly=0; AES_poly < AES_IRRED_POLYNOMIALS; AES_poly++){
		for(AES_gen=0; AES_gen < AES_GENERATORS; AES_gen++){
			GenericAES dualAES;
			dualAES.initFromIndex(AES_poly, AES_gen);

			// write to file
			dumpA   << CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
					<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << ";";
			dumpVector(dumpA, dualAES.sboxAffineGF2E, 256); dumpA << ";";
			dumpVector(dumpA, dualAES.sboxAffineInvGF2E, 256); dumpA << ";";
			dumpMatrix(dumpA, dualAES.mixColMat); dumpA << ";";
			dumpMatrix(dumpA, dualAES.mixColInvMat); dumpA << ";";
			dumpA << endl;
			dumpA.flush();

			cout << "+";
			int ii,qq,probAll;
			for(qq=0;qq<8; qq++){
				cout << ".";
				for(ii=1;ii<256;ii++){
					unsigned long int id = (AES_poly << 24) | (AES_gen << 16) | (qq << 8) | ii;
					int problems=0;
					vec_GF2E A1;
					vec_GF2E A2;
					dualAES.generateA1A2Relations(A1, A2, ii, qq);
					problems = dualAES.testA1A2Relations(A1, A2);

					// write to file
					dump    << CHEX(id) << ";"
							<< CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
							<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << ";"
							<< qq << ";" << ii << ";" << problems << ";";
					dumpVector(dump, A1);
					dump << ";";
					dumpVector(dump, A2);
					dump << endl;

					if (problems>0){
						cout << "!!!Current Dual AES: "
								<< CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
								<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << endl;
						cout << "!!!Problem with relations ii="<<ii<<"; qq="<<qq<<"; problems=" << problems << endl;
						probAll+=1;
					}

					//
					// A1 A2 duplicity check via unordered hashed structure
					//
					a1a2rec rec;
					rec.count=1;
					rec.id = id;
					rec.vec = NULL;

					std::string a1hash = hashLookupTable(A1);
					std::string a2hash = hashLookupTable(A2);

					// A1 hashing
					if (a1store.count(a1hash)==0){
						a1store[a1hash] = rec;		// new A1 hash
					} else {
						a1a2rec &trec = a1store.at(a1hash);	// already contains some record
						trec.count++;
						if (trec.vec==NULL) trec.vec = new vector<unsigned long int>();
						trec.vec->push_back(id);

						if (false)
						cout << "A1: Spoiler, duplicate of ["
								<< CHEX8(id) << ";"
								<< CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
								<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << ";"
								<< qq << ";" << ii << ";] collide with ["

								<< CHEX8(trec.id) << ";"
								<< CHEX(GenericAES::irreduciblePolynomials[(trec.id >> 24) & 0xff]) << ";"
								<< CHEX(GenericAES::generators[(trec.id >> 24) & 0xff][(trec.id >> 16) & 0xff]) << ";"
								<< ((trec.id >> 16) & 0xff) << ";" << ((trec.id) & 0xff) << ";] count="

								<< trec.count << endl;
					}

					// A2 hashing
					if (a2store.count(a2hash)==0){
						a2store[a2hash] = rec;		// new A2 hash
					} else {
						a1a2rec &trec = a2store.at(a2hash);	// already contains some record
						trec.count++;
						if (trec.vec==NULL) trec.vec = new vector<unsigned long int>();
						trec.vec->push_back(id);

						cout << "A2: Spoiler, duplicate of ["
							<< CHEX8(id) << ";"
							<< CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
							<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << ";"
							<< qq << ";" << ii << ";] collide with ["

							<< CHEX8(trec.id) << ";"
							<< CHEX(GenericAES::irreduciblePolynomials[(trec.id >> 24) & 0xff]) << ";"
							<< CHEX(GenericAES::generators[(trec.id >> 24) & 0xff][(trec.id >> 16) & 0xff]) << ";"
							<< ((trec.id >> 16) & 0xff) << ";" << ((trec.id) & 0xff) << ";] count="

							<< trec.count << endl;
					}
				}
			}

			// force write
			dump.flush();
		}
	}



	//
	// Hashing completed, harvest our results...
	//
	dumpA << "======================================================" << endl << "CollisionsA1" << endl;
	BOOST_FOREACH( a1a2map::value_type & v, a1store ) {
		if (v.second.count==1) continue;
		a1a2rec &trec = v.second;

		dumpA << "colisions=" << trec.count << " ["
			<< CHEX8(trec.id) << ";"
			<< CHEX(GenericAES::irreduciblePolynomials[(trec.id >> 24) & 0xff]) << ";"
			<< CHEX(GenericAES::generators[(trec.id >> 24) & 0xff][(trec.id >> 16) & 0xff]) << ";"
			<< ((trec.id >> 16) & 0xff) << ";" << ((trec.id) & 0xff) << "] ";
		vector<unsigned long int> * pvec = trec.vec;
		for(std::vector<unsigned long int>::iterator it = pvec->begin(); it != pvec->end(); ++it) {
		    unsigned long int id = *it;
		    dumpA << "["
				<< CHEX8(id) << ";"
				<< CHEX(GenericAES::irreduciblePolynomials[(id >> 24) & 0xff]) << ";"
				<< CHEX(GenericAES::generators[(id >> 24) & 0xff][(id >> 16) & 0xff]) << ";"
				<< ((id >> 16) & 0xff) << ";" << ((id) & 0xff) << "] ";
		}
		dumpA << endl;
	}

	dumpA << "======================================================" << endl << "CollisionsA2" << endl;
	BOOST_FOREACH( a1a2map::value_type & v, a2store ) {
		if (v.second.count==1) continue;
		a1a2rec &trec = v.second;

		dumpA << "colisions=" << trec.count << " ["
			<< CHEX8(trec.id) << ";"
			<< CHEX(GenericAES::irreduciblePolynomials[(trec.id >> 24) & 0xff]) << ";"
			<< CHEX(GenericAES::generators[(trec.id >> 24) & 0xff][(trec.id >> 16) & 0xff]) << ";"
			<< ((trec.id >> 16) & 0xff) << ";" << ((trec.id) & 0xff) << "] ";
		vector<unsigned long int> * pvec = trec.vec;
		for(std::vector<unsigned long int>::iterator it = pvec->begin(); it != pvec->end(); ++it) {
			unsigned long int id = *it;
			dumpA << "["
				<< CHEX8(id) << ";"
				<< CHEX(GenericAES::irreduciblePolynomials[(id >> 24) & 0xff]) << ";"
				<< CHEX(GenericAES::generators[(id >> 24) & 0xff][(id >> 16) & 0xff]) << ";"
				<< ((id >> 16) & 0xff) << ";" << ((id) & 0xff) << "] ";
		}
		dumpA << endl;
	}
	dumpA.close();
	dump.close();
	return 0;
}

int dualAESTest(void){
	// very poor PRNG seeding, but just for now
	srand((unsigned)time(0));
	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	GenericAES defAES;
	defAES.init(0x11B, 0x03);
	defAES.printAll();
	defAES.testWithVectors();

	GenericAES dualAES;
	dualAES.init(0x11D, 0x9d);
	//dualAES.initFromIndex(15,5);
	dualAES.printAll();
	dualAES.testWithVectors();

	// try round key expansion
	vec_GF2E roundKey;
	vec_GF2E key;
	key.SetLength(128);
	dualAES.expandKey(roundKey, key, KEY_SIZE_16);
	cout << "Round key for ZERO key for 16B: " << endl;
	dumpVector(roundKey);

	mat_GF2E state(INIT_SIZE, 4, 4);
	cout << "Plaintext: " << endl;
	dumpMatrix(state);

	cout << "Testing encryption: " << endl;
	dualAES.encryptInternal(state, roundKey);
	dumpMatrix(state);

	dualAES.applyTinv(state);
	cout << "Testing encryption AFTER Tinv: " << endl;
	dumpMatrix(state);

	dualAES.applyT(state);
	cout << "Testing backward decryption: " << endl;
	dualAES.decryptInternal(state, roundKey);
	dumpMatrix(state);

	cout << "Multiplication matrix: " << endl;
	mat_GF2 multA= dualAES.makeMultAMatrix(2);
	dumpMatrix(multA);

	cout << "Squaring matrix: " << endl;
	mat_GF2 sqr = dualAES.makeSquareMatrix(1);
	dumpMatrix(sqr);


	cout << "A1 and A2 relations, testing all possible" << endl;
	int ii,qq,probAll;
	for(qq=0;qq<8; qq++){
		for(ii=1;ii<256;ii++){
			int problems=0;
			vec_GF2E A1;
			vec_GF2E A2;
			dualAES.generateA1A2Relations(A1, A2, ii, qq);
			problems = dualAES.testA1A2Relations(A1, A2);

			if (problems>0){
				cout << "Problem with relations ii="<<ii<<"; qq="<<qq<<"; problems=" << problems << endl;
				probAll+=1;
			}
		}
	}

	cout << "All relations tested, problemsAll = " << probAll << endl;

	vec_GF2E A1;
	vec_GF2E A2;
	dualAES.generateA1A2Relations(A1, A2, 1+(rand() % 0xfe), rand() % 7);
	cout << "Testing relations A1 A2: Problems = " << dualAES.testA1A2Relations(A1, A2) << endl;

	cout << "A1: " << endl;
	dumpVector(A1);

	cout << "A2: " << endl;
	dumpVector(A2);

	cout << "Generating random bijections: " << endl;
	vec_GF2X rndB;
	vec_GF2X rndBinv;
	generateRandomBijection(rndB, rndBinv, AES_FIELD_SIZE, AES_FIELD_DIM);
	dumpVector(rndB);
	dumpVector(rndBinv);

	return 0;
}

int MBgen(void){
	long i,j;

	// sample matrix stuff
	mat_GF2 A;
	A.SetDims(QSIZE, QSIZE);
	cout << "I will show you a nice matrix: " << endl << A << endl << endl;

	// now generate invertible matrix
	i = generateInvertiblePM(A, QSIZE);
	if (i>=0){
		cout << "found invertible matrix in [" << i << "] iterations: " << endl << A << endl << endl;
	} else {
		cout << "Invertible matrix was not found" << endl;
	}

	// Now try bigger matrix - mixing bijection 8x8
	i = generateInvertiblePM(A, 8);
	if (i>=0){
		cout << "found invertible matrix in [" << i << "] iterations: " << endl << A << endl << endl;
	} else {
		cout << "Invertible matrix was not found" << endl;
	}

	// now try our Gauss method returning modification matrix P
	mat_GF2 B;
	B.SetDims(QSIZE, QSIZE);
	long const Bdata[][4] = {
	  {1, 0, 1, 0},
	  {1, 1, 0, 0},
	  {0, 0, 1, 1},
	  {0, 0, 1, 1}
	};

	initMatrix(B, (long *)Bdata);
	cout << "My custom init matrix: " << endl << B << endl << endl;
	cout << "Now trying to find inverse: " << endl;

	GF2 d;
	ref_GF2 dd(d);
	mat_GF2 P;
	mat_GF2 Q;
	i = invP(dd,P,Q,B);
	cout << "Determinant returned: " << d << "; Rank = " << i;
	cout << "; P matrix: " << endl << P << endl << endl;
	cout << "; Q matrix: " << endl << Q << endl << endl;
	cout << "; R matrix: " << endl << (P*B*Q) << endl << endl;

	// A matrix generation
	mat_GF2 Amat;

	generateARankMatrix(Amat, 1, 8);
	cout << "Generating A matrix, r=1, p=8" << endl << Amat << endl <<endl;

	generateARankMatrix(Amat, 2, 8);
	cout << "Generating A matrix, r=2, p=8" << endl << Amat << endl <<endl;

	generateARankMatrix(Amat, 3, 8);
	cout << "Generating A matrix, r=3, p=8" << endl << Amat << endl <<endl;

	generateARankMatrix(Amat, 4, 8);
	cout << "Generating A matrix, r=4, p=8" << endl << Amat << endl <<endl;

	generateARankMatrix(Amat, 5, 8);
	cout << "Generating A matrix, r=5, p=8" << endl << Amat << endl <<endl;

	mat_GF2 MB;

	// simple mixing bijection invertibility test, 100 iterations
	for(j=0; j<1000; j++){
		generateMixingBijection(MB, 32, 4);

		mat_GF2 MBinv;
		inv(MBinv, MB);
		cout << "## Test passed: " << j << endl;
		cout << "MB: " << MB << endl << endl;
	}
	return EXIT_SUCCESS;
}


