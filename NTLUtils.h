/*
 * NTLUtils.h
 *
 *  Created on: Mar 2, 2013
 *      Author: ph4r05
 */

#ifndef NTLUTILS_H_
#define NTLUTILS_H_

// NTL dependencies
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
NTL_CLIENT


#define CHEX(x) "0x" << setw(2) << setfill('0') << hex << (x)
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

inline long getLong(NTL::GF2X& x) {
    return x.xrep[0];
}

inline long getLong(NTL::GF2E& x) {
    return x.LoopHole().xrep.length() == 0 ? 0L : x.LoopHole().xrep[0];
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

void dumpVector(NTL::vec_GF2E& a);
void dumpVector(NTL::vec_GF2& a);
void dumpVector(NTL::GF2E * a, size_t len);
void dumpVector(long * a, size_t len);
void dumpMatrix(NTL::mat_GF2E& a);
void dumpMatrix(NTL::mat_GF2& a);

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
