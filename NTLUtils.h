/*
 * NTLUtils.h
 *
 *  Created on: Mar 2, 2013
 *      Author: ph4r05
 */

#ifndef NTLUTILS_H_
#define NTLUTILS_H_

// NTL dependencies
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <math.h>
NTL_CLIENT


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
    NTL::GF2X x; x.SetLength(len); x.xrep[0] = a; return x;
}

inline long getLong(NTL::GF2X& x) {
    return x.xrep[0];
}

inline long getLong(NTL::GF2E& x) {
    return getLong(x.LoopHole());
}


#endif /* NTLUTILS_H_ */
