
#ifndef NTL_ZZ_pE__H
#define NTL_ZZ_pE__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/vec_long.h>
#include <NTL/ZZ_pX.h>

NTL_OPEN_NNS


class ZZ_pEInfoT {
private:
   ZZ_pEInfoT();                       // disabled
   ZZ_pEInfoT(const ZZ_pEInfoT&);   // disabled
   void operator=(const ZZ_pEInfoT&);  // disabled
public:
   long ref_count;

   ZZ_pEInfoT(const ZZ_pX&);
   ~ZZ_pEInfoT() { }

   ZZ_pXModulus p;

   ZZ   _card;
   long _card_init;
   ZZ   _card_base;
   long _card_exp;


};

extern ZZ_pEInfoT *ZZ_pEInfo; // info for current modulus, initially null




class ZZ_pEContext {
private:
ZZ_pEInfoT *ptr;

public:
void save();
void restore() const;

ZZ_pEContext() { ptr = 0; }
ZZ_pEContext(const ZZ_pX& p);

ZZ_pEContext(const ZZ_pEContext&); 


ZZ_pEContext& operator=(const ZZ_pEContext&); 


~ZZ_pEContext();


};


class ZZ_pEBak {
private:
long MustRestore;
ZZ_pEInfoT *ptr;

ZZ_pEBak(const ZZ_pEBak&); // disabled
void operator=(const ZZ_pEBak&); // disabled

public:
void save();
void restore();

ZZ_pEBak() { MustRestore = 0; ptr = 0; }

~ZZ_pEBak();


};



struct ZZ_pE_NoAlloc_type { ZZ_pE_NoAlloc_type() { } };
const ZZ_pE_NoAlloc_type ZZ_pE_NoAlloc = ZZ_pE_NoAlloc_type();



class ZZ_pE {

public:

ZZ_pX _ZZ_pE__rep;

// static data


static long DivCross() { return 16; }
static long ModCross() { return 8; }


// ****** constructors and assignment

ZZ_pE();

ZZ_pE(const ZZ_pE& a)  { _ZZ_pE__rep.rep.SetMaxLength(ZZ_pE::degree()); _ZZ_pE__rep = a._ZZ_pE__rep; }

ZZ_pE(ZZ_pE_NoAlloc_type) { }  // allocates no space

~ZZ_pE() { } 

ZZ_pE& operator=(const ZZ_pE& a) { _ZZ_pE__rep = a._ZZ_pE__rep; return *this; }

inline ZZ_pE& operator=(long a);
inline ZZ_pE& operator=(const ZZ_p& a);

ZZ_pE(ZZ_pE& x, INIT_TRANS_TYPE) : _ZZ_pE__rep(x._ZZ_pE__rep, INIT_TRANS) { }


// You can always access the _ZZ_pE__representation directly...if you dare.
ZZ_pX& LoopHole() { return _ZZ_pE__rep; }

static const ZZ_pXModulus& modulus() { return ZZ_pEInfo->p; }

static long degree() { return deg(ZZ_pEInfo->p); }

static const ZZ& cardinality();

static const ZZ_pE& zero();

static long initialized() { return (ZZ_pEInfo != 0); }

static void init(const ZZ_pX&);


};




// read-only access to _ZZ_pE__representation
inline const ZZ_pX& rep(const ZZ_pE& a) { return a._ZZ_pE__rep; }

inline void clear(ZZ_pE& x)
// x = 0
   { clear(x._ZZ_pE__rep); }

inline void set(ZZ_pE& x)
// x = 1
   { set(x._ZZ_pE__rep); }

inline void swap(ZZ_pE& x, ZZ_pE& y)
// swap x and y

   { swap(x._ZZ_pE__rep, y._ZZ_pE__rep); }

// ****** addition

inline void add(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a + b

   { add(x._ZZ_pE__rep, a._ZZ_pE__rep, b._ZZ_pE__rep); }

inline void sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a - b

   { sub(x._ZZ_pE__rep, a._ZZ_pE__rep, b._ZZ_pE__rep); }


inline void negate(ZZ_pE& x, const ZZ_pE& a) 

   { negate(x._ZZ_pE__rep, a._ZZ_pE__rep); }


inline void add(ZZ_pE& x, const ZZ_pE& a, long b)
   { add(x._ZZ_pE__rep, a._ZZ_pE__rep, b); }

inline void add(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { add(x._ZZ_pE__rep, a._ZZ_pE__rep, b); }

inline void add(ZZ_pE& x, long a, const ZZ_pE& b)
   { add(x._ZZ_pE__rep, a, b._ZZ_pE__rep); }

inline void add(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { add(x._ZZ_pE__rep, a, b._ZZ_pE__rep); }





inline void sub(ZZ_pE& x, const ZZ_pE& a, long b)
   { sub(x._ZZ_pE__rep, a._ZZ_pE__rep, b); }

inline void sub(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { sub(x._ZZ_pE__rep, a._ZZ_pE__rep, b); }

inline void sub(ZZ_pE& x, long a, const ZZ_pE& b)
   { sub(x._ZZ_pE__rep, a, b._ZZ_pE__rep); }

inline void sub(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { sub(x._ZZ_pE__rep, a, b._ZZ_pE__rep); }





// ****** multiplication

inline void mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b)
// x = a*b

   { MulMod(x._ZZ_pE__rep, a._ZZ_pE__rep, b._ZZ_pE__rep, ZZ_pE::modulus()); }


inline void sqr(ZZ_pE& x, const ZZ_pE& a)
// x = a^2

   { SqrMod(x._ZZ_pE__rep, a._ZZ_pE__rep, ZZ_pE::modulus()); }

inline ZZ_pE sqr(const ZZ_pE& a)
   { ZZ_pE x; sqr(x, a); NTL_OPT_RETURN(ZZ_pE, x); }


inline void mul(ZZ_pE& x, const ZZ_pE& a, long b)
   { mul(x._ZZ_pE__rep, a._ZZ_pE__rep, b); }

inline void mul(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b)
   { mul(x._ZZ_pE__rep, a._ZZ_pE__rep, b); }

inline void mul(ZZ_pE& x, long a, const ZZ_pE& b)
   { mul(x._ZZ_pE__rep, a, b._ZZ_pE__rep); }

inline void mul(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b)
   { mul(x._ZZ_pE__rep, a, b._ZZ_pE__rep); }


// ****** division



void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_pE& b);
void div(ZZ_pE& x, const ZZ_pE& a, long b);
void div(ZZ_pE& x, const ZZ_pE& a, const ZZ_p& b);
void div(ZZ_pE& x, long a, const ZZ_pE& b);
void div(ZZ_pE& x, const ZZ_p& a, const ZZ_pE& b);

void inv(ZZ_pE& x, const ZZ_pE& a);

inline ZZ_pE inv(const ZZ_pE& a)
   { ZZ_pE x; inv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }



// ****** exponentiation

inline void power(ZZ_pE& x, const ZZ_pE& a, const ZZ& e)
// x = a^e

   { PowerMod(x._ZZ_pE__rep, a._ZZ_pE__rep, e, ZZ_pE::modulus()); }

inline ZZ_pE power(const ZZ_pE& a, const ZZ& e)
   { ZZ_pE x; power(x, a, e); NTL_OPT_RETURN(ZZ_pE, x); }

inline void power(ZZ_pE& x, const ZZ_pE& a, long e)
   { power(x, a, ZZ_expo(e)); }

inline ZZ_pE power(const ZZ_pE& a, long e)
   { ZZ_pE x; power(x, a, e); NTL_OPT_RETURN(ZZ_pE, x); }




// ****** conversion

inline void conv(ZZ_pE& x, const ZZ_pX& a)
   { rem(x._ZZ_pE__rep, a, ZZ_pE::modulus()); }

inline void conv(ZZ_pE& x, long a)
   { conv(x._ZZ_pE__rep, a); }

inline void conv(ZZ_pE& x, const ZZ_p& a)
   { conv(x._ZZ_pE__rep, a); }

inline void conv(ZZ_pE& x, const ZZ& a)
   { conv(x._ZZ_pE__rep, a); }

inline ZZ_pE to_ZZ_pE(const ZZ_pX& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE to_ZZ_pE(long a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE to_ZZ_pE(const ZZ_p& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE to_ZZ_pE(const ZZ& a) 
   { ZZ_pE x; conv(x, a); NTL_OPT_RETURN(ZZ_pE, x); }



// ****** comparison

inline long IsZero(const ZZ_pE& a)
   { return IsZero(a._ZZ_pE__rep); }

inline long IsOne(const ZZ_pE& a)
   { return IsOne(a._ZZ_pE__rep); }

inline long operator==(const ZZ_pE& a, const ZZ_pE& b)
   { return a._ZZ_pE__rep == b._ZZ_pE__rep; }
inline long operator==(const ZZ_pE& a, long b)
   { return a._ZZ_pE__rep == b; }
inline long operator==(const ZZ_pE& a, const ZZ_p& b)
   { return a._ZZ_pE__rep == b; }
inline long operator==(long a, const ZZ_pE& b)
   { return a == b._ZZ_pE__rep; }
inline long operator==(const ZZ_p& a, const ZZ_pE& b)
   { return a == b._ZZ_pE__rep; }

inline long operator!=(const ZZ_pE& a, const ZZ_pE& b)
   { return !(a == b); }
inline long operator!=(const ZZ_pE& a, long b)
   { return !(a == b); }
inline long operator!=(const ZZ_pE& a, const ZZ_p& b)
   { return !(a == b); }
inline long operator!=(long a, const ZZ_pE& b)
   { return !(a == b); }
inline long operator!=(const ZZ_p& a, const ZZ_pE& b)
   { return !(a == b); }


// ****** norm and trace

inline void trace(ZZ_p& x, const ZZ_pE& a)
   { TraceMod(x, a._ZZ_pE__rep, ZZ_pE::modulus()); }
inline ZZ_p trace(const ZZ_pE& a)
   { return TraceMod(a._ZZ_pE__rep, ZZ_pE::modulus()); }

inline void norm(ZZ_p& x, const ZZ_pE& a)
   { NormMod(x, a._ZZ_pE__rep, ZZ_pE::modulus()); }
inline ZZ_p norm(const ZZ_pE& a)
   { return NormMod(a._ZZ_pE__rep, ZZ_pE::modulus()); }


// ****** random numbers

inline void random(ZZ_pE& x)
// x = random element in ZZ_pE

   { random(x._ZZ_pE__rep, ZZ_pE::degree()); }

inline ZZ_pE random_ZZ_pE()
   { ZZ_pE x; random(x); NTL_OPT_RETURN(ZZ_pE, x); }


// ****** input/output

inline NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const ZZ_pE& a)
   { return s << a._ZZ_pE__rep; }
   
NTL_SNS istream& operator>>(NTL_SNS istream& s, ZZ_pE& x);


inline ZZ_pE& ZZ_pE::operator=(long a) { conv(*this, a); return *this; }
inline ZZ_pE& ZZ_pE::operator=(const ZZ_p& a) { conv(*this, a); return *this; }



inline ZZ_pE operator+(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(const ZZ_pE& a, long b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator+(long a, const ZZ_pE& b) 
   { ZZ_pE x; add(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }


inline ZZ_pE operator-(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_pE& a, long b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(long a, const ZZ_pE& b) 
   { ZZ_pE x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator-(const ZZ_pE& a)
   { ZZ_pE x; negate(x, a); NTL_OPT_RETURN(ZZ_pE, x); } 


inline ZZ_pE& operator+=(ZZ_pE& x, const ZZ_pE& b)
   { add(x, x, b); return x; }

inline ZZ_pE& operator+=(ZZ_pE& x, const ZZ_p& b)
   { add(x, x, b); return x; }

inline ZZ_pE& operator+=(ZZ_pE& x, long b)
   { add(x, x, b); return x; }


inline ZZ_pE& operator-=(ZZ_pE& x, const ZZ_pE& b)
   { sub(x, x, b); return x; }

inline ZZ_pE& operator-=(ZZ_pE& x, const ZZ_p& b)
   { sub(x, x, b); return x; }

inline ZZ_pE& operator-=(ZZ_pE& x, long b)
   { sub(x, x, b); return x; }


inline ZZ_pE& operator++(ZZ_pE& x) { add(x, x, 1); return x; }

inline void operator++(ZZ_pE& x, int) { add(x, x, 1); }

inline ZZ_pE& operator--(ZZ_pE& x) { sub(x, x, 1); return x; }

inline void operator--(ZZ_pE& x, int) { sub(x, x, 1); }



inline ZZ_pE operator*(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(const ZZ_pE& a, long b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator*(long a, const ZZ_pE& b) 
   { ZZ_pE x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }


inline ZZ_pE& operator*=(ZZ_pE& x, const ZZ_pE& b)
   { mul(x, x, b); return x; }

inline ZZ_pE& operator*=(ZZ_pE& x, const ZZ_p& b)
   { mul(x, x, b); return x; }

inline ZZ_pE& operator*=(ZZ_pE& x, long b)
   { mul(x, x, b); return x; }




inline ZZ_pE operator/(const ZZ_pE& a, const ZZ_pE& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(const ZZ_pE& a, const ZZ_p& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(const ZZ_pE& a, long b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(const ZZ_p& a, const ZZ_pE& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }

inline ZZ_pE operator/(long a, const ZZ_pE& b) 
   { ZZ_pE x; div(x, a, b); NTL_OPT_RETURN(ZZ_pE, x); }


inline ZZ_pE& operator/=(ZZ_pE& x, const ZZ_pE& b)
   { div(x, x, b); return x; }

inline ZZ_pE& operator/=(ZZ_pE& x, const ZZ_p& b)
   { div(x, x, b); return x; }

inline ZZ_pE& operator/=(ZZ_pE& x, long b)
   { div(x, x, b); return x; }



/* additional legacy conversions for v6 conversion regime */

inline void conv(ZZ_pX& x, const ZZ_pE& a) { x = rep(a); }
inline void conv(ZZ_pE& x, const ZZ_pE& a) { x = a; }


/* ------------------------------------- */



NTL_CLOSE_NNS

#endif
