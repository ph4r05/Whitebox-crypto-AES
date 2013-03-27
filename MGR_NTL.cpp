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
	return hashString(inputBuffer);
}

std::string hashSmap(const smap & map){
	std::string inputBuffer = LinearAffineEq::dumpMapT(map, false);
	return hashString(inputBuffer);
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
	int total=0;
	bsetElem a,b;
	unordered_set<std::string> hashes;
	for(a=0; a<256; a++){
		cout << "+++++++++++++++++++++++++++++ @@[ " << a << "]" << endl;
		for(b=0; b<256; b++){
			cout << "........................... ##[ " << b << "]" << endl;

			bsetElem S2[256];
			bsetElem S2inv[256];
			bsetElem S1[256];
			bsetElem S1inv[256];
			for(i=0; i<256; i++) {
				S1[i]      = defAES.sboxAffine[i ^ a];
				S1inv[S1[i]] = i;
				S2[i]    = defAES.sboxAffine[i] ^ b;
				S2inv[S2[i]] = i;

				//S1[i]      = defAES.sboxAffineInv[i ^ a];
				//S1inv[S1[i]] = i;
				//S2[i]      = defAES.sboxAffineInv[i] ^ b;
				//S2inv[S2[i]] = i;
			}

			LinearAffineEq eqCheck;
			eqCheck.setDimension(8);
			eqCheck.verbosity=0;
			linearEquivalencesList resultList;

			int result = eqCheck.findLinearEquivalences(S1, S1inv, S2, S2inv, &resultList);
			total += result;

			cout << "Done, result = " << result <<  " ; count= " << eqCheck.relationsCount <<  " listSize: " << resultList.size() <<  endl;
			cout << "So far we have [" << total << "]" << endl;

			// Check them now, convert to real affine lookup tables
			smap L1, L2;
			int ch = 0;
			for(linearEquivalencesList::iterator it = resultList.begin(); it != resultList.end(); ++it){
				linearEquiv_t & el = *it;
				cout << "Checking [" << ch << "]" << endl;

				bool same=true;
				bool ok = true;
				for(i=0; i<256; i++){
					if (same && el.TbinvV[i] != el.TaV[i]) same=false;
					bsetElem desired = S2[i];
					bsetElem myVal = el.TbinvV[S1[el.TaV[i]]];

					L1[i] = el.TaV[i] ^ a;
					L2[i] = el.TbinvV[i] ^ b;

					if (desired != myVal){
						ok = false;
						cout << " ! Error, mismatch S2["<<i<<"]=" << desired
							 <<	" vs. B^{-1}[S1[A[i]]=" << myVal << endl;
					}
				}

				if (same) cout << "A1 == B^{-1}" << endl;
				else	  cout << "A1 != B^{-1}" << endl;

				if (ok)  cout << "It holds!!" << endl;
				else	 cout << "BROKEN" << endl;

				// hash it
				std::string hashL1 = hashSmap(L1);
				std::string hashL2 = hashSmap(L2);
				std::string totalHash = hashL1;
				totalHash.append(";").append(hashL2);

				cout << "Total hash: " << totalHash << endl;

				if (hashes.count(totalHash)>0){
					cout << "Already in hash set, skipping" << endl;
				} else {
					hashes.insert(totalHash);
				}
			}

			cout << "So far unique affine relations: " << hashes.size() << endl;
		}
	}

	cout << "Total affine relations: " << total << endl;
	cout << "Total unique affine relations: " << hashes.size() << endl;
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


