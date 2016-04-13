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
#include <stdexcept>
#include <boost/io/ios_state.hpp>

NTL_CLIENT

using namespace std;
using namespace NTL;

void dumpArray(ostream& out, const char * input, size_t len){
	size_t i;
	boost::io::ios_flags_saver ifs(out);
	for (i=0; i<len; i++) {
		out << setw(2) << setfill('0') << hex << ((unsigned int)input[i]&0xff) << " ";
	}
}

void dumpVector(NTL::vec_GF2E& a){
	unsigned int i, len = a.length();
	boost::io::ios_flags_saver ifs(cout);
	for (i=0; i<len; i++){
		cout << " " << GF2EHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpVector(NTL::vec_GF2X& a){
	unsigned int i, len = a.length();
	boost::io::ios_flags_saver ifs(cout);
	for (i=0; i<len; i++){
		cout << " " << GF2EHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpVector(NTL::vec_GF2& a){
	unsigned int i, len = a.length();
	boost::io::ios_flags_saver ifs(cout);
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
	boost::io::ios_flags_saver ifs(cout);
	for (i=0; i<len; i++){
		cout << " " << GF2EHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

void dumpMatrix(NTL::mat_GF2E& a){
	unsigned int i,j, n = a.NumRows(), m=a.NumCols();
	boost::io::ios_flags_saver ifs(cout);
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			cout << " " << GF2EHEX(a[i][j]) << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void dumpMatrix(ostream& out, NTL::mat_GF2& a, bool newLine){
	boost::io::ios_flags_saver ifs(out);
	int i,j, n = a.NumRows(), m=a.NumCols();
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			if (newLine) out << " ";
			out << a.get(i,j);
			if (newLine) out << " ";
			else out << ", ";
		}
		if (newLine) out << endl;
		else out << " | ";
	}
	if (newLine) out << endl;
}

void dumpMatrix(NTL::mat_GF2& a){
	dumpMatrix(cout, a, true);
}

void dumpMatrixN(NTL::mat_GF2 a){
	dumpMatrix(cout, a, true);
}

std::string dumpMatrix2str(NTL::mat_GF2& a, bool newLine){
	std::ostringstream out;
	dumpMatrix(out, a, newLine);
	return out.str();
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
	boost::io::ios_flags_saver ifs(out);
	for (i=0; i<n; i++){
		for(j=0; j<m; j++){
			out << GF2EHEX(a[i][j]) << ",";
		}
		out << "|";
	}
}

void dumpVector(ostream& out, NTL::vec_GF2E& a){
	unsigned int i, len = a.length();
	boost::io::ios_flags_saver ifs(out);
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
	boost::io::ios_flags_saver ifs(out);
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

void applyLookupTable(vec_GF2E& ltable, GF2E& tgt){
	tgt = ltable[getLong(tgt)];
}

void applyLookupTable(vec_GF2E& ltable, vec_GF2E& tgt){
	int i, n=tgt.length();
	for(i=0; i<n; i++){
		tgt[i] = ltable[getLong(tgt[i])];
	}
}

void applyLookupTable(vec_GF2E& ltable, mat_GF2E& tgt){
	int i,j,n=tgt.NumRows(),m=tgt.NumCols();
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			tgt[i][j] = ltable[getLong(tgt[i][j])];
		}
	}
}

std::string hashString(std::string inputBuffer){
	char buff[inputBuffer.size()];					// c++0x feature
	inputBuffer.copy(buff, inputBuffer.size(), 0);	// copy serialized lookup table to char array

	MD5_CTX mdContext;
	MD5Init(&mdContext);
	MD5Update(&mdContext, (unsigned char * )buff, inputBuffer.size());
	MD5Final(&mdContext);

	// hexcoding
	ostringstream ss; ss << "0x";
	int i=0; for(i=0; i<MD5_DIGEST_LENGTH; i++) ss << setw(2) << setfill('0') << hex << ((unsigned long int)mdContext.digest[i]);
	inputBuffer = ss.str();
	return inputBuffer;
}

// hashes lookup table using MD5 hash
std::string hashLookupTable(vec_GF2E s){
	std::string inputBuffer = dumpVector2str(s);
	return hashString(inputBuffer);
}

std::string hashMatrix(mat_GF2 m){
	std::string inputBuffer = dumpMatrix2str(m, false);
	return hashString(inputBuffer);
}

int char2int(char input) {
	if(input >= '0' && input <= '9')
		return input - '0';
	if(input >= 'A' && input <= 'F')
		return input - 'A' + 10;
	if(input >= 'a' && input <= 'f')
		return input - 'a' + 10;
	throw std::invalid_argument("Invalid input string");
}

// This function assumes src to be a zero terminated sanitized string with
// an even number of [0-9a-f] characters, and target to be sufficiently large
size_t hex2bin(const char* src, char* target, size_t maxLen) {
	size_t i;
	for(i = 0; *src && src[1] && i < maxLen;){
		if (*src==' '){
			src+=1;
			continue;
		}
		*(target++) = (char2int(*src)<<4) + char2int(src[1]);
		src += 2;
		i++;
	}

	return i;
}

size_t hexstr2bin(std::string hex, char * target, size_t maxLen){
	return hex2bin(hex.c_str(), target, maxLen);
}
