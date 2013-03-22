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
	dumpVectorT<long>(a, len);
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

void dumpMatrix(ostream& out, NTL::mat_GF2E& a){
	unsigned int i,j, n = a.NumRows(), m=a.NumCols();
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			out << GF2EHEX(a[i][j]) << ",";
		}
		out << "|";
	}
}

void dumpVector(ostream& out, NTL::vec_GF2E& a){
	unsigned int i, len = a.length();
	for (i=0; i<len; i++){
		out << GF2EHEX(a[i]) << ",";
	}
}

std::string dumpVector2str(NTL::vec_GF2E& a){
	 std::ostringstream out;
	 dumpVector(out, a);
	 return out.str();
}

void dumpVector(ostream& out, NTL::GF2E * a, size_t len){
	unsigned int i;
	for (i=0; i<len; i++){
		out << GF2EHEX(a[i]) << ",";
	}
}

void matrix2vector(const NTL::mat_GF2E& src, NTL::vec_GF2E& dst, bool byRows){
	int i,j, n = byRows ? src.NumRows() : src.NumCols(), m = byRows ? src.NumCols() : src.NumRows();

	dst.SetLength(n*m);
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			dst.put(i*m + j, byRows ? src.get(i, j) : src.get(j, i));
		}
	}
}

void vector2matrix(const NTL::vec_GF2E& src, NTL::mat_GF2E& dst, int rowLen, bool byRows){
	int i,n=src.length();
	dst.SetDims(n/rowLen, rowLen);
	for(i=0; i<n; i++){
		if (byRows)
			dst.put(i / rowLen, i % rowLen, src.get(i));
		else
			dst.put(i % rowLen, i / rowLen, src.get(i));
	}
}

void charArr_to_vec_GF2E(const unsigned char * arr, size_t len, NTL::vec_GF2E& dst){
	unsigned int j;
	dst.SetLength(len);
	for(j=0; j<len; j++){
		dst.put(j, GF2EFromLong((unsigned long)arr[j], 8));
	}
}
