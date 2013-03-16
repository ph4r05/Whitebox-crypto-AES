/*
 * NTLUtils.c
 *
 *  Created on: Mar 2, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include <ctime>


// NTL dependencies
#include "NTLUtils.h"
#include <iomanip>
#include <iostream>
#include <fstream>
NTL_CLIENT

using namespace std;
using namespace NTL;

void dumpVector(NTL::vec_GF2E& a){
	unsigned int i, len = a.length();
	for (i=0; i<len; i++){
		cout << " " << GF2EHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpVector(NTL::vec_GF2X& a){
	unsigned int i, len = a.length();
	for (i=0; i<len; i++){
		cout << " " << GF2EHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpVector(NTL::vec_GF2& a){
	unsigned int i, len = a.length();
	for (i=0; i<len; i++){
		cout << " " << a[i] << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpVector(long * a, size_t len){
	unsigned int i;
	for (i=0; i<len; i++){
		cout << " " << CHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpVector(GF2E * a, size_t len){
	unsigned int i;
	for (i=0; i<len; i++){
		cout << " " << GF2EHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpMatrix(NTL::mat_GF2E& a){
	unsigned int i,j, n = a.NumRows(), m=a.NumCols();
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			cout << " " << GF2EHEX(a[i][j]) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void dumpMatrix(NTL::mat_GF2& a){
	int i,j, n = a.NumRows(), m=a.NumCols();
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			cout << " " << a.get(i,j) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

/**
 * Initializes matrix from specified array.
 * Data must be at least dimension of given matrix.
 */
int initMatrix(mat_GF2& M, long *data){
	long i,j,n,m;
	for(i=0, n=M.NumRows(); i<n; i++){
		for(j=0, m=M.NumCols(); j<m; j++){
			M.put(i,j,data[n*i+j]);
		}
	}

	return 0;
}

void dumpMatrix(ofstream& out, NTL::mat_GF2E& a){
	unsigned int i,j, n = a.NumRows(), m=a.NumCols();
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			out << GF2EHEX(a[i][j]) << ",";
		}
		out << "|";
	}
}

void dumpVector(ofstream& out, NTL::vec_GF2E& a){
	unsigned int i, len = a.length();
	for (i=0; i<len; i++){
		out << GF2EHEX(a[i]) << ",";
	}
}

void dumpVector(ofstream& out, NTL::GF2E * a, size_t len){
	unsigned int i;
	for (i=0; i<len; i++){
		out << GF2EHEX(a[i]) << ",";
	}
}

