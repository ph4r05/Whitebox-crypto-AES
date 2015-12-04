//============================================================================
// Name        : MGR_NTL.cpp
// Author      : Dusan Klinec (ph4r05)
// Version     :
// Copyright   : Your copyright notice
// Description : Tests for implemented routines. Tests MixingBijection generator,
//				   A1A2 relations generator and BGE attack
// *  Author: Dusan Klinec (ph4r05)
// *
// *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
//============================================================================
#include "testing.h"

// hardcoded elements
// http://stackoverflow.com/questions/2236197/c-easiest-way-to-initialize-an-stl-vector-with-hardcoded-elements

NTL_CLIENT
using namespace std;
using namespace NTL;
using namespace boost;
using namespace wbacr;
using namespace wbacr::laeqv;
using namespace wbacr::attack;

#define GENERIC_AES_DEBUG 1

int main(void) {
	// very poor PRNG seeding, but just for now
	time_t start, end;

	GF2X defaultModulus = GF2XFromLong(0x11B, 9);
	GF2E::init(defaultModulus);

	GenericAES defAES;
	defAES.init(0x11B, 0x03);
	defAES.printAll();
	WBAESGenerator generator;
	WBAES * genAES = new WBAES;

	// Test WB AES with test vectors.
	// This test also demonstrates usage of external encodings by wrapping AES
	generator.testWithVectors(true, genAES);

	// Invert test
	BGEAttack atk;

	// BGE attack
	time(&start);
	cout << "Starting an attack! Obj: " << endl;
	atk.run();

	time(&end);
	cout << "Computation took: " << (end - start) << " seconds." << endl;
	delete genAES;

	exit(3);

	// Cipher inversion - works only without external encodings,
	//  more specifically, big diffusion matrix has to be I_{128}
	//cout << "Testing cipher inversion;" << endl;
	//atk.invertCipherTest();
	//exit(3);

	//generator.useDualAESIdentity=true;
	//generator.useIO04x04Identity=true;
	//generator.useIO08x08Identity=true;
	//generator.useMB08x08Identity=true;
	//generator.useMB32x32Identity=true;
	//int errors = generator.testWithVectors(true, genAES);
	//cout << "Testing done, errors: " << errors << endl;
	//exit(3);

	//mat_GF2 m = defAES.makeMultAMatrix(0x44);
	//mat_GF2 minv = inv(m);

	//mat_GF2 t1 = defAES.makeSquareMatrix(0x1);
	//dumpMatrix(t1);

	//mat_GF2 t1inv = inv(t1);
	//dumpMatrix(t1inv);

	//mat_GF2 t2 = defAES.makeSquareMatrix(7);
	//dumpMatrix(t2);

	//mat_GF2 res = m*t1;
	//dumpMatrix(res);

	//res = t1inv*minv;
	//dumpMatrix(res);

	//res = m*t1*t1inv*minv;
	//dumpMatrix(res);


	exit(4);

	//dualAESTest();
	//exit(3);

	//A1A2relationsGenerator();
	//exit(2);

	//AESAffineRelationsVerify(true);
	//exit(3);

	//MBgen();
	//exit(3);

	return 0;
}

int A1A2relationsGenerator(void){
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
			int ii=0,qq=0,probAll=0;
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
	int ii=0,qq=0,probAll=0;
	for(qq=0;qq<8; qq++){
		for(ii=1;ii<256;ii++){
			int problems=0;
			vec_GF2E A1;
			vec_GF2E A2;
			dualAES.generateA1A2Relations(A1, A2, ii, qq);
			problems = dualAES.testA1A2Relations(A1, A2, true, false);

			if (problems!=0){
				cout << "Problem with ENC relations ii="<<ii<<"; qq="<<qq<<"; problems=" << problems << endl;
				probAll+=1;
			}

			problems = dualAES.testA1A2Relations(A2, A1, false, false);
			if (problems!=0){
				cout << "Problem with DEC relations ii="<<ii<<"; qq="<<qq<<"; problems=" << problems << endl;
				probAll+=1;
			}
		}
	}

	cout << "All relations tested, problemsAll = " << probAll << endl;

	vec_GF2E A1;
	vec_GF2E A2;
	dualAES.generateA1A2Relations(A1, A2, 1+(phrand() % 0xfe), phrand() % 7);
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

//
// Find affine equivalences for Sboxes for default AES
//
int AESAffineRelationsVerify(bool inverseSbox){
	int i;
	unordered_set<std::string> hashes;

	GenericAES defAES;
	defAES.init(0x11B, 0x03);

	LinearAffineEq eqCheck;
	eqCheck.setDimension(8);
	eqCheck.verbosity=0;
	eqCheck.verbosityAffine=1;
	eqCheck.randomizeXGuess=false;

	bsetElem S2[256];
	bsetElem S2inv[256];
	bsetElem S1[256];
	bsetElem S1inv[256];
	for(i=0; i<256; i++) {
		S1[i] = (!inverseSbox) ? defAES.sboxAffine[i] : defAES.sboxAffineInv[i];
		S2[i] = (!inverseSbox) ? defAES.sboxAffine[i] : defAES.sboxAffineInv[i];
		S1inv[S1[i]] = i;
		S2inv[S2[i]] = i;
	}

	// Print out some usefull information
	cout << "#" << endl
	     << "# Finding affine equivalences between" << endl;

	if (inverseSbox) cout << "# Sinv(x + b) ~~ Sinv(x) + a" << endl;
	else             cout << "# S(x + a) ~~ S(x) + b" << endl;

	cout
	     << "#" << endl
		 << "# new iteration of a value (a=3):"  << endl
		 << "# +++++++++++++++++++++++++++++ @@[ 3]"  << endl
		 << "# new iteration of a value (b=8):"  << endl
		 << "# ........................... ##[ 8]"  << endl
		 << "#"  << endl
		 << "# Idx to matrix correspondence:"  << endl
		 << "# +-------------------------------------------------------+"  << endl
		 << "# |  0  |  1  |     2    |     3    |     4    |     5    |"  << endl
		 << "# | mQM | mMQ | mAQMAinv | mAinvQMA | mAMQAinv | mAinvMQA |"  << endl
		 << "# +-------------------------------------------------------+"  << endl
		 << "#" << endl << endl;

	//
	// Generate ordinary matrix representation & find equivalences, for each matrix, generate hashes
	//
	AESAffineMap amap;

	int j, k, l, m;
	mat_GF2 tmpMat1, tmpMat2, tmpMat3, tmpMat4;

	// 5. S-BOX with affine mappings. At first obtain default form of affine transformation for normal AES
	mat_GF2 tmpSboxAffMatrix(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);		// T * Affine * Tinv
	mat_GF2 tmpSboxAffConst(INIT_SIZE, AES_FIELD_DIM, 1);					// column vector
	mat_GF2 tmpSboxAffMatrixDec(INIT_SIZE, AES_FIELD_DIM, AES_FIELD_DIM);	// T * AffineDec * Tinv
	mat_GF2 tmpSboxAffConstDec(INIT_SIZE, AES_FIELD_DIM, 1);				// column dec vector
	tmpSboxAffMatrix = defAES.getDefaultAffineMatrix();
	tmpSboxAffConst = colVector(defAES.getDefaultAffineConst());
	tmpSboxAffMatrixDec = defAES.getDefaultAffineMatrixDec();
	tmpSboxAffConstDec = colVector(defAES.getDefaultAffineConstDec());

	// prepare squaring matrices
	cout << "Generating square matrices ..." << endl;
	mat_GF2 squares[8];
	for(j=0; j<8; j++){
		squares[j] = defAES.makeSquareMatrix(j);
	}

	// mult & squares
	cout << "Generating multiplication & squared matrices ..." << endl;
	for(k=1; k<256; k++){
		tmpMat2 = defAES.makeMultAMatrix(k);

		cout << '+';
		cout.flush();
		for(j=0; j<8; j++){
			cout << '.';
			cout.flush();

			smap mQM, mMQ, mAQMAinv, mAinvQMA, mAMQAinv, mAinvMQA;
			tmpMat1 = squares[j];

			mat_GF2 QM = tmpMat1 * tmpMat2;
			mat_GF2 MQ = tmpMat2 * tmpMat1;
			eqCheck.buildLookupTableAndCheck(QM, 0, mQM);
			eqCheck.buildLookupTableAndCheck(MQ, 0, mMQ);

			// Now generate transformation for L1
			for(l=0; l<AES_FIELD_SIZE; l++){
				GF2E tmpElem = GF2EFromLong(l, AES_FIELD_DIM);
				GF2E transformedElem;
				mat_GF2 resMatrix;

				// mAQMAinv
				resMatrix = tmpSboxAffMatrix*(QM * ((tmpSboxAffMatrixDec * colVector(tmpElem, AES_FIELD_DIM)) + tmpSboxAffConstDec)) + tmpSboxAffConst;
				colVector(transformedElem, resMatrix, 0);
				mAQMAinv[l] = (bsetElem) getLong(transformedElem);

				// mAMQAinv
				resMatrix = tmpSboxAffMatrix*(MQ * ((tmpSboxAffMatrixDec * colVector(tmpElem, AES_FIELD_DIM)) + tmpSboxAffConstDec)) + tmpSboxAffConst;
				colVector(transformedElem, resMatrix, 0);
				mAMQAinv[l] = (bsetElem) getLong(transformedElem);

				// mAinvQMA
				resMatrix = tmpSboxAffMatrixDec*(QM * ((tmpSboxAffMatrix * colVector(tmpElem, AES_FIELD_DIM)) + tmpSboxAffConst)) + tmpSboxAffConstDec;
				colVector(transformedElem, resMatrix, 0);
				mAinvQMA[l] = (bsetElem) getLong(transformedElem);

				// mAinvMQA
				resMatrix = tmpSboxAffMatrixDec*(MQ * ((tmpSboxAffMatrix * colVector(tmpElem, AES_FIELD_DIM)) + tmpSboxAffConst)) + tmpSboxAffConstDec;
				colVector(transformedElem, resMatrix, 0);
				mAinvMQA[l] = (bsetElem) getLong(transformedElem);
			}

			//mQM, mMQ, mAQMAinv, mAinvQMA, mAMQAinv, mAinvMQA
			std::string mhash[6] = {
					LinearAffineEq::hashSmap(mQM), LinearAffineEq::hashSmap(mMQ),
					LinearAffineEq::hashSmap(mAQMAinv), LinearAffineEq::hashSmap(mAinvQMA),
					LinearAffineEq::hashSmap(mAMQAinv), LinearAffineEq::hashSmap(mAinvMQA) };

			for (m=0; m<6; m++){
				AESAffineElement nEl;
				nEl.square = j;
				nEl.multi = k;
				nEl.type=m;
				amap.insert(AESAffineMapElem(mhash[m], nEl));
			}
		}
	}

	cout << endl << "Done, starting affine relations finder" << endl;

	// Launch main affine equivalences finding algorithm
	affineEquivalencesList list;
	return eqCheck.findAffineEquivalences(S1, S1inv, S2, S2inv, &list, inverseSbox, &AffCallbackCorrespondence, &amap);
}

//
// Find affine relations correspondence for AES Sboxes
//
int AffCallbackCorrespondence(wbacr::laeqv::affineEquiv_t * el, wbacr::laeqv::affineEquivalencesList * lish, boost::unordered_set<std::string> * hashes, wbacr::laeqv::LinearAffineEq * eqCheck, void * usrData){
	std::string hashL1 = LinearAffineEq::hashSmap(el->L1);
	std::string hashL2 = LinearAffineEq::hashSmap(el->L2);
	std::string totalHash = hashL1;
	totalHash.append(";").append(hashL2);

	// try to determine form of the matrix
	// In inversion case, S2 should be linear, try several options
	std::string lhash[2] = {hashL1, hashL2};
	std::string L1str = dumpMatrix2str(el->linPart.Ta, false);
	std::string L2str = dumpMatrix2str(el->linPart.Tbinv, false);

	// convert usrData to map
	AESAffineMap * amap = (AESAffineMap * ) usrData;

	int n=0;
	for(n=0; n<2; n++){
		if (amap->count(lhash[n])>0){
			auto its = amap->equal_range(lhash[n]);
			for (auto it = its.first; it != its.second; ++it) {
				cout << "We have match! L"<<(n+1)
						<<" ~ idx="<<(it->second.type)
						<<" ; [a]="<< setw(2) << (it->second.multi)
						<<" ; Q^i="<<(it->second.square)
						<<" ; L1 const=" << (el->a)
						<<" ; L2 const=" << (el->b) << endl;

				if (el->a==0 && el->b==0){
					cout << "L1 matrix: " << endl;
					dumpMatrix(el->linPart.Ta);

					cout << "L2 matrix: " << endl;
					dumpMatrix(el->linPart.Tbinv);
				}
			}
		}
	}

	return 0;
}
