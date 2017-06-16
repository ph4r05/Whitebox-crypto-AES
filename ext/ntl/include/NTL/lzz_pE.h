
#ifndef NTL_zz_pE__H
#define NTL_zz_pE__H

#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/vec_long.h>
#include <NTL/lzz_pX.h>

NTL_OPEN_NNS



class zz_pEInfoT {
private:
   zz_pEInfoT();                       // disabled
   zz_pEInfoT(const zz_pEInfoT&);   // disabled
   void operator=(const zz_pEInfoT&);  // disabled
public:
   long ref_count;

   zz_pEInfoT(const zz_pX&);
   ~zz_pEInfoT() { }

   zz_pXModulus p;

   ZZ   _card;
   long _card_init;
   long _card_base;
   long _card_exp;


};

extern zz_pEInfoT *zz_pEInfo; // info for current modulus, initially null




class zz_pEContext {
private:
zz_pEInfoT *ptr;

public:
void save();
void restore() const;

zz_pEContext() { ptr = 0; }
zz_pEContext(const zz_pX& p);

zz_pEContext(const zz_pEContext&); 


zz_pEContext& operator=(const zz_pEContext&); 


~zz_pEContext();


};


class zz_pEBak {
private:
long MustRestore;
zz_pEInfoT *ptr;

zz_pEBak(const zz_pEBak&); // disabled
void operator=(const zz_pEBak&); // disabled

public:
void save();
void restore();

zz_pEBak() { MustRestore = 0; ptr = 0; }

~zz_pEBak();


};



struct zz_pE_NoAlloc_type { zz_pE_NoAlloc_type() { } };
const zz_pE_NoAlloc_type zz_pE_NoAlloc = zz_pE_NoAlloc_type();



class zz_pE {

public:

zz_pX _zz_pE__rep;


static long DivCross() { return 16; }
static long ModCross() { return 8; }


// ****** constructors and assignment

zz_pE();

zz_pE(const zz_pE& a)  { _zz_pE__rep.SetMaxLength(zz_pE::degree()); _zz_pE__rep = a._zz_pE__rep; }

zz_pE(zz_pE_NoAlloc_type) { }  // allocates no space

~zz_pE() { } 

zz_pE& operator=(const zz_pE& a) { _zz_pE__rep = a._zz_pE__rep; return *this; }

zz_pE(zz_pE& x, INIT_TRANS_TYPE) : _zz_pE__rep(x._zz_pE__rep, INIT_TRANS) { }


// You can always access the _zz_pE__representation directly...if you dare.
zz_pX& LoopHole() { return _zz_pE__rep; }


static const zz_pXModulus& modulus() { return zz_pEInfo->p; }

static long degree() { return deg(zz_pEInfo->p); }

static const ZZ& cardinality();

static const zz_pE& zero();

static long initialized() { return (zz_pEInfo != 0); }

static void init(const zz_pX&);

inline zz_pE& operator=(long a);
inline zz_pE& operator=(const zz_p& a);

 
};

inline const zz_pX& _zz_pE__rep(const zz_pE& a) { return a._zz_pE__rep; }

inline void clear(zz_pE& x)
// x = 0
   { clear(x._zz_pE__rep); }

inline void set(zz_pE& x)
// x = 1
   { set(x._zz_pE__rep); }

inline void swap(zz_pE& x, zz_pE& y)
// swap x and y

   { swap(x._zz_pE__rep, y._zz_pE__rep); }

// ****** addition

inline void add(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a + b

   { add(x._zz_pE__rep, a._zz_pE__rep, b._zz_pE__rep); }

inline void sub(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a - b

   { sub(x._zz_pE__rep, a._zz_pE__rep, b._zz_pE__rep); }


inline void negate(zz_pE& x, const zz_pE& a) 

   { negate(x._zz_pE__rep, a._zz_pE__rep); }


inline void add(zz_pE& x, const zz_pE& a, long b)
   { add(x._zz_pE__rep, a._zz_pE__rep, b); }

inline void add(zz_pE& x, const zz_pE& a, const zz_p& b)
   { add(x._zz_pE__rep, a._zz_pE__rep, b); }

inline void add(zz_pE& x, long a, const zz_pE& b)
   { add(x._zz_pE__rep, a, b._zz_pE__rep); }

inline void add(zz_pE& x, const zz_p& a, const zz_pE& b)
   { add(x._zz_pE__rep, a, b._zz_pE__rep); }





inline void sub(zz_pE& x, const zz_pE& a, long b)
   { sub(x._zz_pE__rep, a._zz_pE__rep, b); }

inline void sub(zz_pE& x, const zz_pE& a, const zz_p& b)
   { sub(x._zz_pE__rep, a._zz_pE__rep, b); }

inline void sub(zz_pE& x, long a, const zz_pE& b)
   { sub(x._zz_pE__rep, a, b._zz_pE__rep); }

inline void sub(zz_pE& x, const zz_p& a, const zz_pE& b)
   { sub(x._zz_pE__rep, a, b._zz_pE__rep); }





// ****** multiplication

inline void mul(zz_pE& x, const zz_pE& a, const zz_pE& b)
// x = a*b

   { MulMod(x._zz_pE__rep, a._zz_pE__rep, b._zz_pE__rep, zz_pE::modulus()); }


inline void sqr(zz_pE& x, const zz_pE& a)
// x = a^2

   { SqrMod(x._zz_pE__rep, a._zz_pE__rep, zz_pE::modulus()); }

inline zz_pE sqr(const zz_pE& a)
   { zz_pE x; sqr(x, a); NTL_OPT_RETURN(zz_pE, x); }


inline void mul(zz_pE& x, const zz_pE& a, long b)
   { mul(x._zz_pE__rep, a._zz_pE__rep, b); }

inline void mul(zz_pE& x, const zz_pE& a, const zz_p& b)
   { mul(x._zz_pE__rep, a._zz_pE__rep, b); }

inline void mul(zz_pE& x, long a, const zz_pE& b)
   { mul(x._zz_pE__rep, a, b._zz_pE__rep); }

inline void mul(zz_pE& x, const zz_p& a, const zz_pE& b)
   { mul(x._zz_pE__rep, a, b._zz_pE__rep); }


// ****** division



void div(zz_pE& x, const zz_pE& a, const zz_pE& b);
void div(zz_pE& x, const zz_pE& a, long b);
void div(zz_pE& x, const zz_pE& a, const zz_p& b);
void div(zz_pE& x, long a, const zz_pE& b);
void div(zz_pE& x, const zz_p& a, const zz_pE& b);

void inv(zz_pE& x, const zz_pE& a);

inline zz_pE inv(const zz_pE& a)
   { zz_pE x; inv(x, a); NTL_OPT_RETURN(zz_pE, x); }



// ****** exponentiation

inline void power(zz_pE& x, const zz_pE& a, const ZZ& e)
// x = a^e

   { PowerMod(x._zz_pE__rep, a._zz_pE__rep, e, zz_pE::modulus()); }

inline zz_pE power(const zz_pE& a, const ZZ& e)
   { zz_pE x; power(x, a, e); NTL_OPT_RETURN(zz_pE, x); }

inline void power(zz_pE& x, const zz_pE& a, long e)
   { power(x, a, ZZ_expo(e)); }

inline zz_pE power(const zz_pE& a, long e)
   { zz_pE x; power(x, a, e); NTL_OPT_RETURN(zz_pE, x); }




// ****** conversion

inline void conv(zz_pE& x, const zz_pX& a)
   { rem(x._zz_pE__rep, a, zz_pE::modulus()); }

inline void conv(zz_pE& x, long a)
   { conv(x._zz_pE__rep, a); }

inline void conv(zz_pE& x, const zz_p& a)
   { conv(x._zz_pE__rep, a); }

inline void conv(zz_pE& x, const ZZ& a)
   { conv(x._zz_pE__rep, a); }

inline zz_pE to_zz_pE(const zz_pX& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE to_zz_pE(long a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE to_zz_pE(const zz_p& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE to_zz_pE(const ZZ& a) 
   { zz_pE x; conv(x, a); NTL_OPT_RETURN(zz_pE, x); }



// ****** comparison

inline long IsZero(const zz_pE& a)
   { return IsZero(a._zz_pE__rep); }

inline long IsOne(const zz_pE& a)
   { return IsOne(a._zz_pE__rep); }

inline long operator==(const zz_pE& a, const zz_pE& b)
   { return a._zz_pE__rep == b._zz_pE__rep; }
inline long operator==(const zz_pE& a, long b)
   { return a._zz_pE__rep == b; }
inline long operator==(const zz_pE& a, const zz_p& b)
   { return a._zz_pE__rep == b; }
inline long operator==(long a, const zz_pE& b)
   { return a == b._zz_pE__rep; }
inline long operator==(const zz_p& a, const zz_pE& b)
   { return a == b._zz_pE__rep; }

inline long operator!=(const zz_pE& a, const zz_pE& b)
   { return !(a == b); }
inline long operator!=(const zz_pE& a, long b)
   { return !(a == b); }
inline long operator!=(const zz_pE& a, const zz_p& b)
   { return !(a == b); }
inline long operator!=(long a, const zz_pE& b)
   { return !(a == b); }
inline long operator!=(const zz_p& a, const zz_pE& b)
   { return !(a == b); }


// ****** norm and trace

inline void trace(zz_p& x, const zz_pE& a)
   { TraceMod(x, a._zz_pE__rep, zz_pE::modulus()); }
inline zz_p trace(const zz_pE& a)
   { return TraceMod(a._zz_pE__rep, zz_pE::modulus()); }

inline void norm(zz_p& x, const zz_pE& a)
   { NormMod(x, a._zz_pE__rep, zz_pE::modulus()); }
inline zz_p norm(const zz_pE& a)
   { return NormMod(a._zz_pE__rep, zz_pE::modulus()); }


// ****** random numbers

inline void random(zz_pE& x)
// x = random element in zz_pE

   { random(x._zz_pE__rep, zz_pE::degree()); }

inline zz_pE random_zz_pE()
   { zz_pE x; random(x); NTL_OPT_RETURN(zz_pE, x); }


// ****** input/output

inline NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const zz_pE& a)
   { return s << a._zz_pE__rep; }
   
NTL_SNS istream& operator>>(NTL_SNS istream& s, zz_pE& x);

inline const zz_pX& rep(const zz_pE& a) { return a._zz_pE__rep; }



inline zz_pE& zz_pE::operator=(long a) { conv(*this, a); return *this; }
inline zz_pE& zz_pE::operator=(const zz_p& a) { conv(*this, a); return *this; }




inline zz_pE operator+(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(const zz_pE& a, const zz_p& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(const zz_pE& a, long b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(const zz_p& a, const zz_pE& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator+(long a, const zz_pE& b) 
   { zz_pE x; add(x, a, b); NTL_OPT_RETURN(zz_pE, x); }


inline zz_pE operator-(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_pE& a, const zz_p& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_pE& a, long b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_p& a, const zz_pE& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(long a, const zz_pE& b) 
   { zz_pE x; sub(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator-(const zz_pE& a)
   { zz_pE x; negate(x, a); NTL_OPT_RETURN(zz_pE, x); } 


inline zz_pE& operator+=(zz_pE& x, const zz_pE& b)
   { add(x, x, b); return x; }

inline zz_pE& operator+=(zz_pE& x, const zz_p& b)
   { add(x, x, b); return x; }

inline zz_pE& operator+=(zz_pE& x, long b)
   { add(x, x, b); return x; }


inline zz_pE& operator-=(zz_pE& x, const zz_pE& b)
   { sub(x, x, b); return x; }

inline zz_pE& operator-=(zz_pE& x, const zz_p& b)
   { sub(x, x, b); return x; }

inline zz_pE& operator-=(zz_pE& x, long b)
   { sub(x, x, b); return x; }


inline zz_pE& operator++(zz_pE& x) { add(x, x, 1); return x; }

inline void operator++(zz_pE& x, int) { add(x, x, 1); }

inline zz_pE& operator--(zz_pE& x) { sub(x, x, 1); return x; }

inline void operator--(zz_pE& x, int) { sub(x, x, 1); }



inline zz_pE operator*(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(const zz_pE& a, const zz_p& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(const zz_pE& a, long b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(const zz_p& a, const zz_pE& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator*(long a, const zz_pE& b) 
   { zz_pE x; mul(x, a, b); NTL_OPT_RETURN(zz_pE, x); }


inline zz_pE& operator*=(zz_pE& x, const zz_pE& b)
   { mul(x, x, b); return x; }

inline zz_pE& operator*=(zz_pE& x, const zz_p& b)
   { mul(x, x, b); return x; }

inline zz_pE& operator*=(zz_pE& x, long b)
   { mul(x, x, b); return x; }




inline zz_pE operator/(const zz_pE& a, const zz_pE& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(const zz_pE& a, const zz_p& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(const zz_pE& a, long b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(const zz_p& a, const zz_pE& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }

inline zz_pE operator/(long a, const zz_pE& b) 
   { zz_pE x; div(x, a, b); NTL_OPT_RETURN(zz_pE, x); }


inline zz_pE& operator/=(zz_pE& x, const zz_pE& b)
   { div(x, x, b); return x; }

inline zz_pE& operator/=(zz_pE& x, const zz_p& b)
   { div(x, x, b); return x; }

inline zz_pE& operator/=(zz_pE& x, long b)
   { div(x, x, b); return x; }



/* additional legacy conversions for v6 conversion regime */

inline void conv(zz_pX& x, const zz_pE& a) { x = rep(a); }
inline void conv(zz_pE& x, const zz_pE& a) { x = a; }


/* ------------------------------------- */


NTL_CLOSE_NNS

#endif
