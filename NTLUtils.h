/*
 * NTLUtils.h
 *
 *  Created on: Mar 2, 2013
 *  Author: Dusan Klinec (ph4r05)
 *
 *  License: GPLv3 [http://www.gnu.org/licenses/gpl-3.0.html]
 */

#ifndef NTLUTILS_H_
#define NTLUTILS_H_

// NTL dependencies
#include "base.h"
#include <NTL/vector.h>
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_GF2E.h>
#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <math.h>
#include <assert.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/io/ios_state.hpp>
#include "md5.h"
NTL_CLIENT


#define CHEX8(x) "0x" << setw(8) << setfill('0') << hex << ((unsigned long int)x)
#define CHEX(x) "0x" << setw(2) << setfill('0') << hex << ((unsigned long int)x)
#define GF2EHEX(x) "0x" << setw(2) << setfill('0') << hex << (getLong(x))

/**
 * Init GF2 polynomial from long.
 * Original conv() version used only LSB from long producing not
 * what we exactly expected.
 */
inline void GF2XFromLong(NTL::GF2X& x, const long a, const long len){
    x.SetLength(len);
    x.xrep[0] = a;
}

inline NTL::GF2X GF2XFromLong(const long a, const long len){
	if (a==0) return GF2X::zero();
    NTL::GF2X x; x.SetLength(len); x.xrep[0] = a;  x.SetLength(len); return x;
}

inline NTL::GF2E GF2EFromLong(const long a, const long len){
    if (a==0) return GF2E::zero();
    NTL::GF2X x; x.SetLength(len); x.xrep[0] = a;  x.SetLength(len); return to_GF2E(x);
}

inline long getLong(const NTL::GF2X& x) {
    return deg(x) < 0 ? 0L : x.xrep[0];
}

inline long getLong(NTL::GF2E& x) {
    return x.LoopHole().xrep.length() == 0 ? 0L : x.LoopHole().xrep[0];
}

inline long getLong(const NTL::GF2E& x) {
	NTL::GF2X tmp; conv(tmp, x);
    return tmp.xrep.length() == 0 ? 0L : tmp.xrep[0];
}

template<typename T>
inline NTL::mat_GF2 colVector(const T c){
	NTL::mat_GF2 ret(INIT_SIZE, 8, 1);
	ret[0][0] = (c & (1<<0)) > 0 ? 1 : 0; 	ret[1][0] = (c & (1<<1)) > 0 ? 1 : 0;
	ret[2][0] = (c & (1<<2)) > 0 ? 1 : 0; 	ret[3][0] = (c & (1<<3)) > 0 ? 1 : 0;
	ret[4][0] = (c & (1<<4)) > 0 ? 1 : 0; 	ret[5][0] = (c & (1<<5)) > 0 ? 1 : 0;
	ret[6][0] = (c & (1<<6)) > 0 ? 1 : 0; 	ret[7][0] = (c & (1<<7)) > 0 ? 1 : 0;
	return ret;
}

template<typename T>
inline void colVectorT(const T c, NTL::mat_GF2 & src, int col){
	src[0][col] = (c & (1<<0)) > 0 ? 1 : 0; 	src[1][col] = (c & (1<<1)) > 0 ? 1 : 0;
	src[2][col] = (c & (1<<2)) > 0 ? 1 : 0; 	src[3][col] = (c & (1<<3)) > 0 ? 1 : 0;
	src[4][col] = (c & (1<<4)) > 0 ? 1 : 0; 	src[5][col] = (c & (1<<5)) > 0 ? 1 : 0;
	src[6][col] = (c & (1<<6)) > 0 ? 1 : 0; 	src[7][col] = (c & (1<<7)) > 0 ? 1 : 0;
}

inline NTL::mat_GF2 colVector(const NTL::vec_GF2& v){
	int i=0, ln = v.length();
	NTL::mat_GF2 ret(INIT_SIZE, ln, 1);
	for(i=0; i < ln; i++) ret.put(i, 0, v[i]); 
	return ret;
}

inline NTL::mat_GF2 colVector(const NTL::GF2X& v, int ln){
    int i=0, realLn = deg(v);
    NTL::mat_GF2 ret(INIT_SIZE, ln, 1);
    for(i=0; i <= realLn; i++) ret.put(i, 0, v[i]); 
    return ret;
}

inline NTL::mat_GF2 colVector(const NTL::GF2E& v, int ln){
    int i=0, realLn;
    NTL::mat_GF2 ret(INIT_SIZE, ln, 1);
    NTL::GF2X tV; conv(tV, v); realLn = deg(tV);
    for(i=0; i <= realLn; i++) ret.put(i, 0, tV[i]); 
    return ret;
}

inline void colVector(NTL::vec_GF2& x, const NTL::mat_GF2& m, int which){
	int i=0, ln=m.NumRows();
	for(i=0; i<ln; i++) x.put(i, m.get(i, which));
}

inline NTL::vec_GF2 colVector_vecGF2(const NTL::mat_GF2& m, int which){
    NTL::vec_GF2 ret(INIT_SIZE, m.NumRows());
    colVector(ret, m, which);
    return ret;
}

inline void colVector(NTL::GF2X& x, const NTL::mat_GF2& m, int which){
    NTL::vec_GF2 vec = colVector_vecGF2(m, which); conv(x, vec);
}

inline NTL::GF2X colVector_GF2X(const NTL::mat_GF2& m, int which){
    NTL::GF2X ret; colVector(ret, m, which); return ret;
}

inline void colVector(NTL::GF2E& x, const NTL::mat_GF2& m, int which){
    NTL::vec_GF2 vec = colVector_vecGF2(m, which);
    GF2X gf2x; conv(gf2x, vec); conv(x, gf2x);
}

inline NTL::GF2E colVector_GF2E(const NTL::mat_GF2& m, int which){
    NTL::GF2E ret;
    colVector(ret, m, which);
    return ret;
}

/**
 * Converts matrix consisting of GF2E elements to binary matrix from
 * element representation, coding binary elements to columns. LSB is first in the row,
 * what is consistent with GenericAES.
 */
inline void mat_GF2E_to_mat_GF2_col(NTL::mat_GF2& dst, const NTL::mat_GF2E& src, int elemLen){
	int i,j,k,n = src.NumRows(),m = src.NumCols();
	dst.SetDims(elemLen * n, m);
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			GF2E curElem = src.get(i,j);
			GF2X curX = curElem.LoopHole();
			int xdeg = deg(curX);
			for(k=0; k<elemLen; k++){
				dst.put(i*elemLen + k, j, k<=xdeg ? curX[k] : GF2::zero());
			}
		}
	}
}

/**
 * Converts 8bit column representation to array of characters
 */
inline void mat_GF2_to_char(NTL::mat_GF2& src, unsigned char * dst){
	int i,j,n=src.NumRows(),m=src.NumCols();
	assert((n%8)!=0);

	for(i=0; i<n/8; i++){
		for(j=0; j<m; j++){
			dst[i*m + j] = 0;
			dst[i*m + j] |=
					  ((src[i*8+0][j] == 1) ? 1    : 0)
					| ((src[i*8+1][j] == 1) ? 1<<1 : 0)
					| ((src[i*8+2][j] == 1) ? 1<<2 : 0)
					| ((src[i*8+3][j] == 1) ? 1<<3 : 0)
					| ((src[i*8+4][j] == 1) ? 1<<4 : 0)
					| ((src[i*8+5][j] == 1) ? 1<<5 : 0)
					| ((src[i*8+6][j] == 1) ? 1<<6 : 0)
					| ((src[i*8+7][j] == 1) ? 1<<7 : 0);
		}
	}
}

void dumpArray(ostream& out, const char * input, size_t len);
void dumpVector(NTL::vec_GF2E& a);
void dumpVector(NTL::vec_GF2X& a);
void dumpVector(NTL::vec_GF2& a);
void dumpVector(NTL::GF2E * a, size_t len);
void dumpVector(long * a, size_t len);
void dumpMatrix(NTL::mat_GF2E& a);
void dumpMatrix(NTL::mat_GF2& a);
void dumpMatrixN(NTL::mat_GF2 a);
void dumpMatrix(ostream& out, NTL::mat_GF2E& a);
void dumpMatrix(ostream& out, NTL::mat_GF2& a, bool newLine);
std::string dumpMatrix2str(NTL::mat_GF2& a, bool newLine);
void dumpVector(ostream& out, NTL::vec_GF2E& a);
void dumpVector(ostream& out, NTL::GF2E * a, size_t len);
std::string dumpVector2str(NTL::vec_GF2E& a);

template<typename T> void dumpVectorT(T * a, size_t len){
	unsigned int i;
	boost::io::ios_flags_saver ifs(cout);
	for (i=0; i<len; i++){
		cout << " " << CHEX(a[i]) << " ";
		if (((i+1) % 16) == 0) cout << endl;
	}
	cout << endl;
}

// Compares two GF2E vectors - long values has to be equal => return true otherwise fase
inline bool compare_vec_GF2E(const NTL::vec_GF2E& a, const NTL::vec_GF2E& b){
	int i,n = a.length();
	if (n!=b.length()) return false;
	for(i=0; i<n; i++){
		if (getLong(a[i]) != getLong(b[i])) return false;
	}

	return true;
}

void matrix2vector(const NTL::mat_GF2E& src, NTL::vec_GF2E& dst, bool byRows);
void vector2matrix(const NTL::vec_GF2E& src, NTL::mat_GF2E& dst, int rowLen, bool byRows);
void charArr_to_vec_GF2E(const unsigned char * arr, size_t len, NTL::vec_GF2E& dst);

void applyLookupTable(vec_GF2E& ltable, GF2E& tgt);
void applyLookupTable(vec_GF2E& ltable, vec_GF2E& tgt);
void applyLookupTable(vec_GF2E& ltable, mat_GF2E& tgt);

/**
 * Initializes matrix from specified array.
 * Data must be at least dimension of given matrix.
 */
int initMatrix(mat_GF2& M, long *data);

size_t hex2bin(const char* src, char* target, size_t maxLen);
size_t hexstr2bin(std::string hex, char * target, size_t maxLen);

/**
 * Hashes input string with MD5 hash
 */
std::string hashString(std::string inputBuffer);

// hashes lookup table using MD5 hash
std::string hashLookupTable(vec_GF2E s);

std::string hashMatrix(mat_GF2 m);

template <class T> class vector_inserter{
public:
    NTL::Vec<T>& v;
    vector_inserter(NTL::Vec<T>& v):v(v){}
    vector_inserter& operator,(const T& val){v.push_back(val);return *this;}
};
template <class T> vector_inserter<T>& operator+=(NTL::Vec<T>& v,const T& x){
    return vector_inserter<T>(v),x;
}

#endif /* NTLUTILS_H_ */
