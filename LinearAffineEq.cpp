/*
 * LinearAffineEq.cpp
 *
 *  Created on: Mar 26, 2013
 *      Author: ph4r05
 */

#include "LinearAffineEq.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

#define WBAES_BOOTS_SERIALIZATION 1

// NTL dependencies
namespace wbacr {
namespace laeqv {

NTL_CLIENT
using namespace std;
using namespace NTL;
using namespace boost;

LinearAffineEq::LinearAffineEq() {
	verbosity=0;
	relationsCount=0;
	size=256;
	dim=8;
}

LinearAffineEq::~LinearAffineEq() {
	;
}

// returns set C = A \ B
bset LinearAffineEq::setDiff(const bset &A, const bset &B){
	bset newSet(A);
	for(bset::const_iterator it = B.begin(); it != B.end(); ++it){
		newSet.erase(*it);
	}

	return newSet;
}

void LinearAffineEq::dumpMap(const smap& mp){
	for(smap::const_iterator it = mp.begin(); it != mp.end(); ++it){
		cout << "  mp[" << setw(2) << (it->first) << "] = " << setw(2) << (it->second) << endl;
	}
}

void LinearAffineEq::dumpSet(const bset s){
	for(bset::const_iterator it = s.begin(); it != s.end(); ++it){
		cout << setw(2) << (*it) << " ";
	}
}

// extract lineary independent vectors from input map - keys
bset LinearAffineEq::extractLinearlyIndependent(const smap& mp){
	bset tmpKeySet;
	bset resKeySet;
	bset spanSpace, spanSpace2;
	smap::const_iterator it = mp.begin();
	bset::const_iterator its;
	int i;
	for(i=0; it != mp.end(); ++it){
		if (it->first == 0) continue;
		if (i==0){
			resKeySet.insert( it->first );
			spanSpace.insert( it->first );
		}
		else
			tmpKeySet.insert( it->first );
		i++;
	}

	// Do while key set is exhausted
	while(tmpKeySet.empty()==false){
		// 1. take next vector in set as a base vector, is linearly independent from
		// vectors in resKeySet. Then remove it from keySet.
		its = tmpKeySet.begin();
		bsetElem newVect = *its;
		tmpKeySet.erase(its);

		// Remove all linear combinations of new vectors and old vectors in resKeySet, add to span space
		spanSpace.insert(newVect);
		resKeySet.insert(newVect);

		// generate new vectors
		spanSpace2 = spanSpace;
		bset::iterator its2 = spanSpace2.end();
		for(its = spanSpace2.begin(); its != its2; ++its){
			spanSpace.insert(newVect ^ (*its));
		}

		// tmpKeySet = tmpKeySet \ spanSpace
		for(its = spanSpace.begin(); its != spanSpace.end(); ++its){
			tmpKeySet.erase(*its);
		}
	}

	return resKeySet;
}

mat_GF2 LinearAffineEq::vectorSet2GF2matrix(const bset & s, int dim){
	mat_GF2 ret(INIT_SIZE, dim, s.size());
	int i=0,j=0;
	for(bset::const_iterator it = s.begin(); it != s.end(); ++it, j++){
		for(i=0; i<dim; i++){
			ret.put(i, j, ((*it) & (1 << i)) > 0 ? 1 : 0);
		}
	}

	return ret;
}

mat_GF2 LinearAffineEq::values2GF2matrix(const bset & s, const smap & m, int dim){
	mat_GF2 ret(INIT_SIZE, dim, s.size());
	int i=0,j=0;
	for(bset::const_iterator it = s.begin(); it != s.end(); ++it, j++){
		const bsetElem e = *it;
		const bsetElem v = m.at(e);
		for(i=0; i<dim; i++){
			ret.put(i, j, (( v ) & (1 << i)) > 0 ? 1 : 0);
		}
	}

	return ret;
}

int LinearAffineEq::buildLookupTableAndCheck(mat_GF2 & Ta, bsetElem cst, smap & mapA){
	for(unsigned int i=0; i<size; i++){
		mat_GF2 e = colVector(i);
		mat_GF2 v = Ta * e;
		GF2X res = colVector_GF2X(v, 0);
		bsetElem resE = getLong(res) + cst;

		// check mapping, if it is correct (consistent table data with new computed data)
		if (mapA.count(i)>0){
			if (mapA[i] != resE) {
				if (verbosity){
					cout << "Mapping inconsistent for i=" << i
						<< "; mapA[i] = " << mapA[i]
						<< "; resE = " << resE << endl;
					cout << "colVector: " << endl;
					dumpMatrix(e);
					cout << " Ta*e " << endl;
					dumpMatrix(v);
					cout << "res: " << res << endl;
				}
				return -4;
			}
		} else {
			mapA.insert(smapElem(i, resE));
		}
	}

	return 0;
}

int LinearAffineEq::checkInvertibleLinear(const bset & Ua,   const bset & Ub,
							  	  	  	  smap & mapA,       smap & mapB,
							  	  	  	  bsetElem * S1,     bsetElem * S1inv,
							  	  	  	  bsetElem * S2,     bsetElem * S2inv,
							  	  	  	  mat_GF2 & Ta, 	 mat_GF2 & Tb){
	// Extract linearly independent vectors, to determine mapping.
	// We need matrix consisting of Avect to be invertible in order to determine
	// matrix representation of transformation.
	bset Avect = LinearAffineEq::extractLinearlyIndependent(mapA);
	if (verbosity) {
		cout << "Size of linearly independent vectors: " << Avect.size() << endl;
		LinearAffineEq::dumpSet(Avect);
	}

	// Derive matrix representation of transformation, if possible
	//
	// Ta * Ainp = Aout => Ta = Aout * Ainp^{-1}
	mat_GF2 Ainp = LinearAffineEq::vectorSet2GF2matrix(Avect, 8), AinpInv;
	mat_GF2 Aout = LinearAffineEq::values2GF2matrix(Avect, mapA, 8);
	mat_GF2 TAinv;
	GF2 det;

	// Dimension check, each matrix has to have exactly 8 rows (dimension) and at least 8 columns
	// (number of equations, sample points)
	if (       Ainp.NumCols() < dim || Ainp.NumRows() != dim
			|| Aout.NumCols() < dim || Aout.NumRows() != dim){
		if (verbosity) cout << "Dimension mismatch for Ainp || Aout matrices " << endl;
		return -1;
	}

	if (verbosity){
		cout << "Input matrix: " << endl;
		dumpMatrix(Ainp);
		cout << "Output matrix: " << endl;
		dumpMatrix(Aout);
	}

	// invertible?
	inv(det, AinpInv, Ainp);
	if (det == 0) {
		if (verbosity) cout << "A Matrix is not invertible! " << endl;
		return -2;
	}

	if (verbosity){
		cout << "Inverse matrix: " << endl;
		dumpMatrix(AinpInv);
	}

	// obtain linear transformation
	Ta = Aout * AinpInv;
	if (verbosity) {
		cout << "Ta matrix: " << endl;
		dumpMatrix(Ta);
	}

	// invertible?
	inv(det, TAinv, Ta);
	if (det==0){
		if (verbosity) cout << "Transformation is not linear & invertible!" << endl;
		return -3;
	}

	if (verbosity){
		cout << "Transformation matrix repr. obtained!!!" << endl;
		dumpMatrix(Ta);
	}

	//
	// A is known (matrix representation), build lookup table
	//
	if (buildLookupTableAndCheck(Ta, 0, mapA)<0){
		return -4;
	}

	//
	// Deriving mapping B from mapping A that is complete now and in matrix form.
	//
	// B * S2 = S1 * A
	// B(x) = S1 * A * S2^{-1} (x)
	// From this we derive mapping for B for basis vectors directly.
	Tb.SetDims(dim,dim);
	for(unsigned int i=0; i<dim; i++){
		bsetElem base = 1 << i;
		bsetElem res = S1[mapA[S2inv[base]]];
		for(unsigned int j=0; j<dim; j++){
			Tb.put(j, i, (res & (1<<j)) > 0 ? 1 : 0);
		}
	}

	if (verbosity){
		cout << "Mapping B derived" << endl;
		dumpMatrix(Tb);
	}

	// Transformation B invertibility test.
	// Inversion is needed for final test for relations properties with S-boxes
	mat_GF2 Tbinv;
	inv(det, Tbinv, Tb);
	if (det==0){
		if (verbosity) cout << "B is not linear invertible transformation, cannot create inversion." << endl;
		return -4;
	}

	// build lookup table for B and check already precomputed values
	if (buildLookupTableAndCheck(Tb, 0, mapB)<0){
		return -4;
	}

	// Whole range test for holding desired properties with Sboxes for which
	// they were designed.
	// For this we also need B^{-1} transformation.
	//
	// We apply this equation here:
	// B^{-1} * S1 * A = S2
	//
	smap mapBinv;
	buildLookupTableAndCheck(Tbinv, 0, mapBinv);

	if (verbosity) cout << "Testing matrix representation with |Ua| = " << Ua.size() << " and |Ub| = " << Ub.size() << endl;
	for(bsetElem iter = 0; iter < size; iter++){
		bsetElem desiredResult = S2[iter];
		bsetElem myResult = mapBinv[S1[mapA[iter]]];
		if (desiredResult!=myResult){
			if (verbosity){
				cout << "Problem with relations, it does not work for: " << iter << endl;
				cout << "S2["<<iter<<"]=           " << desiredResult << endl;
				cout << "B^{-1}[S1[A["<<iter<<"]]]=" << myResult << endl;
			}

			return -5;
		}
	}

	return 0;
}

int LinearAffineEq::findLinearEquivalences(bsetElem * S1, 	bsetElem * S1inv,
										   bsetElem * S2, 		bsetElem * S2inv){
	int count = 0;
	int i;
	bset Ua, Ub, Na, Nb, Ca, Cb;
	bsetElem guesses1[size-1]; int guess1Idx=0;	// random guesses for A(x) mapping
	randomPermutationT(guesses1, size-1, 1);	// make it random - random permutation

	// recursive stack for guesses
	recStack_t recStack;
	recStack.resize(1);							// force stack to contain 1 root element - guess
	recStack.reserve(dim);						// we will need at most <dim> vectors, but 2 should be enough!
	int stackIdx=-1;							// -1 initialization -> has to be initialized in routine

	// Known mapping
	smap mapA, mapB;

	// init Ua, Ub, Na, Nb
	// Our Sboxes don't map 0 to 0, so we can save this mapping.
	Ca.insert(0); Na.insert(0);
	Cb.insert(0); Nb.insert(0);
	mapA.insert(smapElem(0, 0));
	mapB.insert(smapElem(0, 0));

	// Random initialization of Ua, Ub
	bsetElem rndInit[255];
	randomPermutationT(rndInit, 255, 1);
	for(i=0; i<255; i++){
		Ua.insert(rndInit[i]);
		Ub.insert(rndInit[i]);
	}

	bool guessRejected=false;
	bset::const_iterator it1, it2, it3;
	while(Ua.empty()==false && Ub.empty()==false && guess1Idx < 255){
		if (verbosity){
			cout << endl << "===================================================================================="
				 << endl << "Main cycle started " << endl;
		}

		//
		// Starting with new guess
		//
		if (Na.empty() && Nb.empty()){
			//
			// 1. If previous guess rejected, restore Ca, Cb, Ua, Ub
			// Guess A(x) for some x \in Ua
			// Set Na = {x}, Ua = Ua / {x}
			//

			bsetElem x;

			//
			// Guess rejected? recover last backups
			//
			if (guessRejected){
				// Is everything done?
				if (stackIdx == -1) {
					break;
				}

				// Revert backups for state variables
				linEqGuess_t & gs = recStack[stackIdx];
				Ca   = gs.Ca;   Cb   = gs.Cb;
				Ua   = gs.Ua;	Ub   = gs.Ub;
				mapA = gs.mapA; mapB = gs.mapB;

				// select same X as in previous case
				x = gs.guessKey;
				Na.insert(x);
				guessRejected=false;
				gs.idx+=1;

				if (verbosity){
					cout << "#GuessWasRejected" << endl;
				}
			} else {
				// Guess was not rejected
				//  -> descent in recursive stack
				recStack.resize( recStack.size() + 1);
				stackIdx+=1;

				linEqGuess_t & gs = recStack[stackIdx];

				// At first, backup Ca, Cb, Ua, Ub - will be restored in case of incorrect guess
				gs.Ua   = Ua;      gs.Ub   = Ub;
				gs.Ca   = Ca;      gs.Cb   = Cb;
				gs.mapA = mapA;    gs.mapB = mapB;
				gs.idx  = 0;

				// Chose new X and pick value for it
				// Keep in mind linearity of mapping, so avoid duplicities.
				int rnd = rand() % Ua.size();
				it1 = Ua.begin(); for(i=0; i<rnd; ++i, ++it1);
				x = *it1; Ua.erase(it1);
				Na.insert(x);
			}
			linEqGuess_t & gs = recStack[stackIdx];

			// Guess A(x) value, avoid duplicate with guesses in recursive stack
			gs.guessKey = x;
			gs.guessVal = 0;
			for(bool freeGuess=true; gs.idx < size; gs.idx++){
				gs.guessVal = guesses1[gs.idx];
				for(i=0; i<stackIdx; i++) {
					if (recStack[i].guessVal == gs.guessVal) { freeGuess=false; break; }
				}
				if (freeGuess) break;
			}

			// Is possible to guess ? check exhaustion
			if (gs.idx >= (size-1)){
				// terminate this stack level, invalid guess one level above
				guessRejected=true;
				Na.clear(); Nb.clear();
				stackIdx-=1;
				continue;
			}

			mapA.insert(smapElem(gs.guessKey, gs.guessVal));
			if (verbosity){
				cout << "New guess;" << endl;
				cout << "G: x=" << (recStack[stackIdx].guessKey) << endl;
				cout << "G: A(x) = " << (recStack[stackIdx].guessVal) << endl;
				cout << "G: idx = " << gs.idx << endl;
				cout << "G: stackIdx = " << stackIdx << endl;
			}
		}

		//
		// Na cycle
		//
		while(Na.empty() == false){
			bsetElem x;

			if (verbosity) cout << endl << "A: Cycle 1 start, |Na| = " << Na.size() << endl;

			// Pick x \in Na; Na = Na \ {x};
			it1 = Na.begin(); x = *it1; Na.erase(it1);
			if (verbosity) cout << "A: newX is [" << x << "]" << endl;

			// Nb = S2( x + Ca ) \ Cb
			if (Ca.size() > 0){
				bset tmpSet;
				for (it1=Ca.begin(); it1!=Ca.end(); ++it1){
					bsetElem curr = *it1;
					bsetElem tmp = x ^ curr;
					tmpSet.insert(S2[tmp]);
					// Use linearity of A to build mapping for tmp
					mapA.insert(smapElem(tmp, mapA[x] ^ mapA[curr]));

					if (verbosity){
						cout << "    curr=" << setw(2) << (curr)
								<< "; tmp = " << setw(2) << tmp
								<< "; S2[tmp]    = " << setw(2) << S2[tmp] << endl;
						cout << "A:     adding mapping for A(x^curr) = A(" << (tmp) << ") = A(x) ^ A(curr) = " << mapA[x] << " ^ " << mapA[curr] << " = " << (mapA[x] ^ mapA[curr]) << endl;
					}

					//
					// A(x) = S^{-1}_1 (B(S_2(x)))
					// S1 * A = B * S2
					if (Cb.count(S2[tmp])==0){
						const bsetElem amap = mapA[tmp];
						mapB.insert(smapElem(S2[tmp], S1[amap]));

						if (verbosity) cout << "A:     adding mapping for B(S2(tmp)) = B("<< S2[tmp] <<") = S1(A(tmp)) = S1(" << amap << ") = " << S1[amap] << endl;
					}
				}
				Nb = LinearAffineEq::setDiff(tmpSet, Cb);
			}

			// Ca = Ca U (x + Ca)
			bset tmpSet(Ca);
			for(it1=tmpSet.begin(); it1 != tmpSet.end(); ++it1){
				Ca.insert( x ^ (*it1) );
			}
			Ca.insert(x);

			// check linearity and derive
			double vectKnown = Nb.size() + log2(Cb.size());
			if (verbosity) cout << "A: vect knownB: " << vectKnown << endl;
			if (vectKnown > dim){
				mat_GF2 Ta, Tb;
				int result = checkInvertibleLinear(Ua, Ub, mapA, mapB, S1, S1inv, S2, S2inv, Ta, Tb);
				if (result==-1 || result==-2) continue;
				if (result==-3 || result==-4 || result==-5){
					guessRejected=true;
					Na.clear(); Nb.clear();
					break;
				}

				if (verbosity) cout << "Process results here!! Everything OK" << endl;

				count+=1;
				relationsCount+=1;
				guessRejected=true;
				Na.clear(); Nb.clear();
			}
		}
		cout << endl;

		//
		// Nb cycle
		//
		while(Nb.empty() == false){
			bsetElem x;
			it1 = Nb.begin(); x = *it1; Nb.erase(it1);

			if (verbosity){
				cout << "B: Cycle 2 start, |Nb| = " << Nb.size() << endl;
				cout << "B: newX is [" << x << "]" << endl;
			}

			if (Cb.size() > 0){
				bset tmpSet;
				for (it1=Cb.begin(); it1!=Cb.end(); ++it1){
					bsetElem curr = *it1;
					bsetElem tmp = x ^ curr;
					tmpSet.insert(S2inv[tmp]);

					// Use linearity of B to build mapping for tmp
					mapB.insert(smapElem(tmp, mapB[x] ^ mapB[curr]));

					if (verbosity){
						cout << "    curr=" << setw(2) << (curr)
								<< "; tmp = " << setw(2) << tmp
								<< "; S2inv[tmp] = " << setw(2) << S2inv[tmp] << endl;
						cout << "B:     adding mapping for B(x^curr) = B(" << (tmp) << ") = B(x) ^ B(curr) = " << mapB[x] << " ^ " << mapB[curr] << " = " << (mapB[x] ^ mapB[curr]) << endl;
					}

					//
					//       A            = S^{-1}_1 * B * S_2
					//       A * S_2^{-1} = S^{-1}_1 * B
					if (Ca.count(S2inv[tmp])==0){
						const bsetElem bmap = mapB[tmp];
						mapA.insert(smapElem(S2inv[tmp], S1inv[bmap]));

						if (verbosity) cout << "B:     adding mapping for A(S2inv(tmp)) = A(" << S2inv[tmp] << ") = S1inv(B(tmp)) = S1inv(" << bmap << ") = " << (S1inv[bmap]) << endl;
					}
				}
				Na = LinearAffineEq::setDiff(tmpSet, Ca);
			}

			// Cb = Cb U (x + Cb)
			bset tmpSet(Cb);
			for(it1=tmpSet.begin(); it1 != tmpSet.end(); ++it1){
				Cb.insert( x ^ *it1 );
			}
			Cb.insert(x);

			double vectKnown = Na.size() + log2(Ca.size());
			if (verbosity) cout << "B: vect knownA: " << vectKnown << endl;
			if (vectKnown>dim){
				if (verbosity) cout << "B: ## check linearity of A, derive,..." << endl;
				mat_GF2 Ta, Tb;
				int result = checkInvertibleLinear(Ub, Ua, mapB, mapA, S1, S1inv, S2, S2inv, Tb, Ta);
				if (result==-1 || result==-2) continue;
				if (result==-3 || result==-4 || result==-5){
					guessRejected=true;
					Na.clear(); Nb.clear();
					break;
				}

				if (verbosity) cout << "Process results here!! Everything OK" << endl;

				count += 1;
				relationsCount+=1;
				guessRejected=true;
				Na.clear(); Nb.clear();
			}
		}

		Ua = LinearAffineEq::setDiff(Ua, Ca);
		Ub = LinearAffineEq::setDiff(Ub, Cb);

		if (verbosity){
			cout << endl << "EEEnd of both cycles, remove Ca, Cb from Ua Ub" << endl;
			cout    << " |Ua| = " << setw(3) << Ua.size()
					<< " |Ub| = " << setw(3) << Ub.size()
					<< " |Ca| = " << setw(3) << Ca.size()
					<< " |Cb| = " << setw(3) << Cb.size()
					<< " |Na| = " << setw(3) << Na.size()
					<< " |Nb| = " << setw(3) << Nb.size() << endl;
			cout << "Dump mapA: " << endl;
			LinearAffineEq::dumpMap(mapA);
			cout << "Dump mapB: " << endl;
			LinearAffineEq::dumpMap(mapB);
		}
	}

	return count;
}

} /* namespace laeqv */
} /* namespace wbacr */
