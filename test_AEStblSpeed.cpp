//============================================================================
// Name        : test_AEStblSpeed.cpp
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
#include <iostream>
#include <fstream>
NTL_CLIENT

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace NTL;
int main(int argc, const char * argv[]) {
	// very poor PRNG seeding, but just for now
	srand((unsigned)time(0));
	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	// input parameters processing
	po::options_description description("WBAES table implementation usage");
	description.add_options()
		("help,h", "Display this help message")
		("compression,c", po::value<int>()->default_value(5)->implicit_value(10),"Compression level")
		("input-files", po::value<std::vector<std::string>>(), "Input files")
		("create-table", po::value<std::string>(), "Create encryption/decryption tables")
		("version,v", "Display the version number");


    po::positional_options_description p;
    p.add("input-files", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).positional(p).run(), vm);
    po::notify(vm);

    if(vm.count("help")){
        std::cout << description;

        return 0;
    }

    if(vm.count("version")){
        std::cout << "WBAES table implementation version 0.1" << std::endl;
        return 0;
    }

    if(vm.count("compression")){
        std::cout << "Compression level " << vm["compression"].as<int>() << std::endl;
    }

	if(vm.count("input-files")){
		std::vector<std::string> files = vm["input-files"].as<std::vector<std::string>>();
		for(std::string file : files){
			std::cout << "Input file " << file << std::endl;
		}
	}

	//A1A2relationsGenerator();
	//dualAESTest();
	//exit(2);

	/*GenericAES defAES;
	defAES.init(0x11B, 0x03);
	defAES.printAll();

	WBAESGenerator generator;
	WBAES genAES;
	int errors = generator.testWithVectors(true, genAES);

	cout << "Testing done, errors: " << errors << endl;*/

/*
	cout << "===Done===" << endl << "Going to generate WBAES..." << endl;
	WBAESGenerator generator;
	WBAES genAES;
	BYTE aesKey[16]; memset(aesKey, 0, sizeof(BYTE)*16);
	ExtCoding coding;
 	generator.generateExtCoding(&coding, WBAESGEN_EXTGEN_ID);
 	generator.generateTables(aesKey, KEY_SIZE_16, genAES, &coding, true);

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
}

int A1A2relationsGenerator(void){
	// very poor PRNG seeding, but just for now
	srand((unsigned)time(0));
	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	ofstream dump;
	ofstream dumpA;
	dump.open("/media/share/AES_A1A2dump.txt");
	dumpA.open("/media/share/AES_signature.txt");


	dump  << "polynomial;generator;qq;ii;problems;A1;A2" << endl;
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
					int problems=0;
					vec_GF2E A1;
					vec_GF2E A2;
					dualAES.generateA1A2Relations(A1, A2, ii, qq);
					problems = dualAES.testA1A2Relations(A1, A2);

					// write to file
					dump    << CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
							<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << ";"
							<< qq << ";" << ii << ";" << problems << ";";
					dumpVector(dump, A1);
					dump << ";";
					dumpVector(dump, A2);
					dump << endl;

					if (problems>0){
						cout << "Current Dual AES: "
								<< CHEX(GenericAES::irreduciblePolynomials[AES_poly]) << ";"
								<< CHEX(GenericAES::generators[AES_poly][AES_gen]) << endl;
						cout << "Problem with relations ii="<<ii<<"; qq="<<qq<<"; problems=" << problems << endl;
						probAll+=1;
					}
				}
			}

			// force write
			dump.flush();
		}
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



