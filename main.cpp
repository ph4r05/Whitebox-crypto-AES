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
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

// NTL dependencies
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/vec_GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_long.h>

#include "GenericAES.h"
#include "NTLUtils.h"
#include "MixingBijections.h"
#include "WBAES.h"
#include "WBAESGenerator.h"
#include "BGEAttack.h"
#include "InputObjectIstream.h"
#include "InputObjectOstream.h"
#include "EncTools.h"
NTL_CLIENT

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/exception/diagnostic_information.hpp>
namespace po = boost::program_options;

using namespace std;
using namespace NTL;
using namespace boost;
using namespace wbacr;
using namespace wbacr::laeqv;
using namespace wbacr::attack;

int tryMain(int argc, const char * argv[]);
int main(int argc, const char * argv[]) {
	try {
		return tryMain(argc, argv);
	}  catch(...) {
		std::clog << boost::current_exception_diagnostic_information() << std::endl;
		return -1;
	}
}

int tryMain(int argc, const char * argv[]) {
	time_t start=0, end=0;
	bool useExternal = false;
	int benchgen=0;
	int benchbge=0;
	bool randomKey=false;
	bool decrypt=false;
	bool pkcs5Padding=false;
	bool cbc=false;
	std::string outFile;
	std::string outTables;
	std::string inTables;
	std::string aesKey;
	std::string cbcIv;
	unsigned char keyFromString[AES_BYTES];
	unsigned char ivFromString[N_BYTES] = {0};
	unsigned char * keyToUse = GenericAES::testVect128_key;

	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	// input parameters processing
	po::options_description description("WBAES table implementation usage");
	description.add_options()
		("help,h",                                                                         "Display this help message")
		("bench-gen",      po::value<int>()->default_value(0)->implicit_value(0),          "Benchmarking rounds for AES gen")
		("bench-bge",      po::value<int>()->default_value(0)->implicit_value(0),          "Benchmarking rounds for AES BGE attack")
		("extEnc,e",       po::value<bool>()->default_value(false)->implicit_value(false), "Use external encoding?")
		("out-file,o",     po::value<std::string>(),                                       "Output file to write encrypted data")
		("input-files",    po::value<std::vector<std::string>>(),                          "Input files")
		("create-table",   po::value<std::string>(),                                       "Create encryption/decryption tables")
		("create-random",  po::value<bool>()->default_value(false)->implicit_value(false), "Create tables with random key")
		("use-key",        po::value<std::string>(),                                       "Create encryption/decryption with given hex-coded key")
		("use-iv",         po::value<std::string>(),                                       "Use CBC with given hex-coded IV")
		("load-tables",    po::value<std::string>(),                                       "Loads encryption/decryption tables from given file")
		("decrypt",        po::value<bool>()->default_value(false)->implicit_value(false), "Should perform encryption or decryption")
		("pkcs5",          po::value<bool>()->default_value(false)->implicit_value(false), "Enables PKCS5 padding")
		("cbc",            po::value<bool>()->default_value(false)->implicit_value(false), "Uses CBC mode")
		("version,v",                                                                      "Display the version number");


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

    if (vm.count("out-file")){
    	outFile = vm["out-file"].as<std::string>();
    	std::cout << "Output file given: " << outFile << endl;
    }

	if (vm.count("load-tables")){
		inTables = vm["load-tables"].as<std::string>();
		std::cout << "Table input file given: " << inTables << endl;
		keyToUse = nullptr;
	}

	if (vm.count("use-key")){
		aesKey = vm["use-key"].as<std::string>();
		size_t bytes = hexstr2bin(aesKey, (char*)keyFromString, AES_BYTES);

		if (bytes != AES_BYTES){
			std::cerr << "Invalid AES key size, expected " << AES_BYTES << " bytes" << std::endl;
			return -1;
		}

		keyToUse = keyFromString;
	}

	if (vm.count("use-iv")){
		cbcIv = vm["use-iv"].as<std::string>();
		size_t bytes = hexstr2bin(cbcIv, (char*)ivFromString, N_BYTES);

		if (bytes != N_BYTES){
			std::cerr << "Invalid AES IV size, expected " << N_BYTES << " bytes" << std::endl;
			return -1;
		}

		keyToUse = keyFromString;
	}

	// use random key?
	randomKey = vm["create-random"].as<bool>();
	if (randomKey){
		for(int i=0; i<AES_BYTES; i++){
			keyFromString[i] = (unsigned char)(phrand() % 0x100);
		}

		keyToUse = keyFromString;
	}

    if (keyToUse != nullptr){
	    std::cout << "AES key to use: ";
	    dumpArray(std::cout, (char *)keyToUse, AES_BYTES);
	    std::cout << std::endl;
	}

    // use external coding ?
    useExternal = vm["extEnc"].as<bool>();
    decrypt = vm["decrypt"].as<bool>();
    pkcs5Padding = vm["pkcs5"].as<bool>();
	cbc = vm["cbc"].as<bool>();

    //
    // AES generator benchmark
    //
    benchgen = vm["bench-gen"].as<int>();
    if (benchgen > 0 && keyToUse != nullptr){
    	cout << "Benchmark of the WB-AES generator is starting..." << endl;

    	//
		// Encryption with WB AES time test
		//
		GenericAES defAES;
		defAES.init(0x11B, 0x03);

		WBAESGenerator generator;

		cout << "External coding will be used: " << useExternal << endl;
		ExtEncoding coding;
		generator.generateExtEncoding(&coding, useExternal ? 0 : WBAESGEN_EXTGEN_ID);

		cout << "Going to compute tables. Benchmark will iterate " << benchgen << " times" << endl;
		time(&start);

		auto * genAES = new WBAES;
		for(int i=0; i<benchgen; i++){
			generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, true);
			generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, false);

			time(&end);
			cout << "Done i="<<i<<"; time elapsed=" << (end-start) << endl;
		}

		delete genAES;
		time_t total = end-start;

		cout << "Benchmark finished! Total time = " << (total) << " s; on average = " << (total / (float)benchgen) << " s" << endl;
    }

    //
	//BGE benchmark
	//
    benchbge = vm["bench-bge"].as<int>();
    if (benchbge > 0 && keyToUse != nullptr){
    	cout << "Benchmark of the BGE attack is starting..." << endl;

    	//
		// Encryption with WB AES time test
		//
    	time_t total = 0;
		GenericAES defAES;
		defAES.init(0x11B, 0x03);

		WBAESGenerator generator;
		generator.useDualAESARelationsIdentity=false;
		generator.useDualAESIdentity=false;
		generator.useDualAESSimpeAlternate=false;
		generator.useIO04x04Identity=false;
		generator.useIO08x08Identity=false;
		generator.useMB08x08Identity=false;
		generator.useMB32x32Identity=false;

		ExtEncoding coding;
		generator.generateExtEncoding(&coding, WBAESGEN_EXTGEN_ID);

		cout << "Going to generate AES tables to crack ..." << endl;

		auto * genAES = new WBAES;
		auto * tmpAES = new WBAES;
		auto * atk = new BGEAttack;

		clock_t pstart, pend;
		clock_t pacc = 0;

		generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, true);
		generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, false);
		std::string serialized = genAES->save();

		cout << "Tables generated, starting with the attack; size of the AES struct = " << sizeof(WBAES) << endl;
		cout << "Serialized length: " << serialized.length() << endl;
		for(int i = 0; i < benchbge; i++){
			// Copy generated AES with serialization.
			tmpAES->loadString(serialized);
			atk->wbaes = tmpAES;

			// Self-test on the first round.
			if (i == 0 && keyToUse == GenericAES::testVect128_key){
				cout << "Going to test WBAES before modifying tables" << endl;
				generator.testComputedVectors(true, tmpAES, &coding);
			}

			cout << "Starting round="<<i<<endl;

			time(&start);
			pstart = clock();

			atk->attack();

			pend = clock();
			time(&end);

			total += end-start;
			pacc  += pend - pstart;
		}

		delete genAES;
		delete atk;
		delete tmpAES;
		cout << "Benchmark finished! Total time = " << (total) << " s; on average = " << (total / (float)benchbge) << " s"
				<< "; clocktime=" << ((float) pacc / CLOCKS_PER_SEC) << " s;" << endl;
    }

	//
	// Create tables & dump to a file
	//
	if (vm.count("create-table") && keyToUse != nullptr){
		outTables = vm["create-table"].as<std::string>();
		std::cout << "Table output file given: " << outTables << endl;

		GenericAES defAES;
		defAES.init(0x11B, 0x03);

		WBAESGenerator generator;
		auto * genAES = new WBAES;

		cout << "External coding will be used: " << useExternal << endl;
		ExtEncoding coding;
		generator.generateExtEncoding(&coding, useExternal ? 0 : WBAESGEN_EXTGEN_ID);

		cout << "Generating WB-AES instance (encryption)..." << endl;
		time(&start);
		generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, true);
		cout << "Generating WB-AES instance (decryption)..." << endl;
		generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, false);
		time(&end);
		cout << "Generating AES tables took: ["<<(end-start)<<"] seconds" << endl;

		generator.save(outTables.c_str(), genAES, &coding);
	}

    //
    // AES encryption - encrypt input files with table representation
    //
	if(vm.count("input-files")){
		std::vector<std::string>  files = vm["input-files"].as<std::vector<std::string>>();
		for(const std::string &file : files){
			std::cout << "Input file " << file << std::endl;
		}

		//
		// Encryption with WB AES time test
		//
		GenericAES defAES;
		defAES.init(0x11B, 0x03);

		WBAESGenerator generator;
		auto * genAES = new WBAES;
		ExtEncoding coding;
		// Generate new encoding.
		cout << "Generating External encoding, identity: " << useExternal << "..." << endl;
		generator.generateExtEncoding(&coding, useExternal ? 0 : WBAESGEN_EXTGEN_ID);

		if (inTables.empty() && keyToUse != nullptr) {
			cout << "Generating WB-AES instance..." << endl;
			time(&start);
			generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, true);
			generator.generateTables(keyToUse, KEY_SIZE_16, genAES, &coding, false);
			time(&end);
			cout << "Generating AES tables took: [" << (end - start) << "] seconds" << endl;
		} else {
			cout << "Loading stored AES tables: " << useExternal << endl;
			time(&start);
			generator.load(inTables.c_str(), genAES, &coding);
			time(&end);
			cout << "Loading AES tables took: [" << (end - start) << "] seconds" << endl;
		}

		// open the given file
		std::string fileName = files[0];
		cout << "Going to " << (decrypt ? "decrypt":"encrypt") << " file ["<<fileName<<"] with WBAES" << endl;
		cout << "External coding: " << useExternal << ", PKCS5 padding: " << pkcs5Padding << ", CBC: " << cbc << endl;

		bool writeOut = !outFile.empty();
		ofstream out;
		if (writeOut){
			out.open(outFile.c_str(), ios::out | ios::binary | ios::trunc);
		}

		// Open reading file
		ifstream inf(fileName.c_str(), ios::in | ios::binary);
		if (!inf.is_open()){
			cerr << "Cannot open specified input file" << endl;
			exit(3);
		}

        InputObjectIstream<BYTE> iois(&inf);
        InputObjectOstream<BYTE> ioos(&out);

        time_t cacc=0;
        clock_t pacc = 0;
        EncTools::processData(decrypt, genAES, &generator, &iois, writeOut ? &ioos : nullptr, &coding, pkcs5Padding,
                              cbc, ivFromString, &cacc, &pacc);

		time(&end);
		time_t total = end-start;
		cout << "Encryption ended in ["<<total<<"]s; Pure encryption took ["<<((float) pacc / CLOCKS_PER_SEC)
					<<"] s (clock call); time: ["<<cacc<<"] s; " << endl;

		// free allocated memory
		delete genAES;
		// close reading file
		inf.close();
		// close output writing file
		if (writeOut){
			out.flush();
			out.close();
		}
	}

	return 0;
}

