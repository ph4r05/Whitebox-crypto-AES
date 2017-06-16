
#ifndef NTL_zz_p__H
#define NTL_zz_p__H

#include <NTL/ZZ.h>
#include <NTL/FFT.h>

NTL_OPEN_NNS


class zz_pInfoT {
private:
   zz_pInfoT();                      // disabled
   zz_pInfoT(const zz_pInfoT&);  // disabled
   void operator=(const zz_pInfoT&); // disabled
public:
   zz_pInfoT(long NewP, long maxroot);
   zz_pInfoT(long Index);
   ~zz_pInfoT();

   long ref_count;

   long p;
   double pinv;

   long index;        // index >= 0 means we are directly using
                     // an FFT prime

   long PrimeCnt;     // 0 for FFT prime;  otherwise same as NumPrimes
                     // used for establishing crossover points

   long NumPrimes;

   long MaxRoot;

   long MinusMModP;  //  -M mod p, M = product of primes

   // the following arrays are indexed 0..NumPrimes-1
   // q = FFTPrime[i]


   long *CoeffModP;    // coeff mod p

   double *x;          // u/q, where u = (M/q)^{-1} mod q
   long *u;            // u, as above
};

extern zz_pInfoT *zz_pInfo;  // current modulus, initially null



class zz_pContext {
private:
zz_pInfoT *ptr;

public:
void save();
void restore() const;

zz_pContext() { ptr = 0; }
zz_pContext(long p, long maxroot=NTL_FFTMaxRoot);
zz_pContext(INIT_FFT_TYPE, long index);

zz_pContext(const zz_pContext&); 

zz_pContext& operator=(const zz_pContext&); 

~zz_pContext();


};


class zz_pBak {
private:
long MustRestore;
zz_pInfoT *ptr;

zz_pBak(const zz_pBak&); // disabled
void operator=(const zz_pBak&); // disabled

public:
void save();
void restore();

zz_pBak() { MustRestore = 0; ptr = 0; }

~zz_pBak();


};

#define NTL_zz_pRegister(x) zz_p x


class zz_p {
public:

long _zz_p__rep;


static void init(long NewP, long maxroot=NTL_FFTMaxRoot);
static void FFTInit(long index);



// ****** constructors and assignment

zz_p() { _zz_p__rep = 0; }

zz_p(const zz_p& a) :  _zz_p__rep(a._zz_p__rep) { }  

~zz_p() { } 

zz_p& operator=(const zz_p& a) { _zz_p__rep = a._zz_p__rep; return *this; }

inline zz_p& operator=(long a);

// a loop-hole for direct access to _zz_p__rep
long& LoopHole() { return _zz_p__rep; }

static long modulus() { return zz_pInfo->p; }
static zz_p zero() { return zz_p(); }
static double ModulusInverse() { return zz_pInfo->pinv; }
static long PrimeCnt() { return zz_pInfo->PrimeCnt; }


static long storage() { return sizeof(long); }

zz_p(long a, INIT_LOOP_HOLE_TYPE) { _zz_p__rep = a; }

};

zz_p to_zz_p(long a);
void conv(zz_p& x, long a);

inline zz_p& zz_p::operator=(long a) { conv(*this, a); return *this; }

zz_p to_zz_p(const ZZ& a);
void conv(zz_p& x, const ZZ& a);


// read-only access to _zz_p__representation
inline long rep(zz_p a) { return a._zz_p__rep; }

inline void clear(zz_p& x)
// x = 0
   { x._zz_p__rep = 0; }

inline void set(zz_p& x)
// x = 1
   { x._zz_p__rep = 1; }

inline void swap(zz_p& x, zz_p& y)
// swap x and y

   { long t;  t = x._zz_p__rep; x._zz_p__rep = y._zz_p__rep; y._zz_p__rep = t; }

// ****** addition

inline void add(zz_p& x, zz_p a, zz_p b)
// x = a + b

   { x._zz_p__rep = AddMod(a._zz_p__rep, b._zz_p__rep, zz_p::modulus()); }

inline void sub(zz_p& x, zz_p a, zz_p b)
// x = a - b

   { x._zz_p__rep = SubMod(a._zz_p__rep, b._zz_p__rep, zz_p::modulus()); }


inline void negate(zz_p& x, zz_p a)
// x = -a

   { x._zz_p__rep = SubMod(0, a._zz_p__rep, zz_p::modulus()); }

// scalar versions

inline void add(zz_p& x, zz_p a, long b) { add(x, a, to_zz_p(b)); }
inline void add(zz_p& x, long a, zz_p b) { add(x, to_zz_p(a), b); }

inline void sub(zz_p& x, zz_p a, long b) { sub(x, a, to_zz_p(b)); }
inline void sub(zz_p& x, long a, zz_p b) { sub(x, to_zz_p(a), b); }

inline zz_p operator+(zz_p a, zz_p b)
    { zz_p x; add(x, a, b); return x; }

inline zz_p operator+(zz_p a, long b)
    { zz_p x; add(x, a, b); return x; }

inline zz_p operator+(long a, zz_p b)
    { zz_p x; add(x, a, b); return x; }

inline zz_p operator-(zz_p a, zz_p b)
    { zz_p x; sub(x, a, b); return x; }

inline zz_p operator-(zz_p a, long b)
    { zz_p x; sub(x, a, b); return x; }

inline zz_p operator-(long a, zz_p b)
    { zz_p x; sub(x, a, b); return x; }



inline zz_p operator-(zz_p a)
   { zz_p x; negate(x, a); return x; }



inline zz_p& operator+=(zz_p& x, zz_p b)
   { add(x, x, b); return x; }

inline zz_p& operator+=(zz_p& x, long b)
   { add(x, x, b); return x; }



inline zz_p& operator-=(zz_p& x, zz_p b)
   { sub(x, x, b); return x; }

inline zz_p& operator-=(zz_p& x, long b)
   { sub(x, x, b); return x; }

inline zz_p& operator++(zz_p& x) { add(x, x, 1); return x; }
inline void operator++(zz_p& x, int) { add(x, x, 1); }
inline zz_p& operator--(zz_p& x) { sub(x, x, 1); return x; }
inline void operator--(zz_p& x, int) { sub(x, x, 1); }

// ****** multiplication

inline void mul(zz_p& x, zz_p a, zz_p b)
// x = a*b

   { x._zz_p__rep = MulMod(a._zz_p__rep, b._zz_p__rep, zz_p::modulus(), zz_p::ModulusInverse()); }

inline void mul(zz_p& x, zz_p a, long b) { mul(x, a, to_zz_p(b)); }
inline void mul(zz_p& x, long a, zz_p b) { mul(x, to_zz_p(a), b); }

inline zz_p operator*(zz_p a, zz_p b)
    { zz_p x; mul(x, a, b); return x; }

inline zz_p operator*(zz_p a, long b)
    { zz_p x; mul(x, a, b); return x; }

inline zz_p operator*(long a, zz_p b)
    { zz_p x; mul(x, a, b); return x; }


inline zz_p& operator*=(zz_p& x, zz_p b)
   { mul(x, x, b); return x; }

inline zz_p& operator*=(zz_p& x, long b)
   { mul(x, x, b); return x; }



inline void sqr(zz_p& x, zz_p a)
// x = a^2

   { x._zz_p__rep = MulMod(a._zz_p__rep, a._zz_p__rep, zz_p::modulus(), zz_p::ModulusInverse()); }

inline zz_p sqr(zz_p a)
   { zz_p x; sqr(x, a); return x; }



// ****** division

inline void div(zz_p& x, zz_p a, zz_p b)
// x = a/b

   { x._zz_p__rep = MulMod(a._zz_p__rep, InvMod(b._zz_p__rep, zz_p::modulus()), zz_p::modulus(),
                    zz_p::ModulusInverse()); }

inline void inv(zz_p& x, zz_p a)
// x = 1/a

   { x._zz_p__rep = InvMod(a._zz_p__rep, zz_p::modulus()); }

inline zz_p inv(zz_p a)
   { zz_p x; inv(x, a); return x; }

inline void div(zz_p& x, zz_p a, long b) { div(x, a, to_zz_p(b)); }
inline void div(zz_p& x, long a, zz_p b) { div(x, to_zz_p(a), b); }

inline zz_p operator/(zz_p a, zz_p b)
    { zz_p x; div(x, a, b); return x; }

inline zz_p operator/(zz_p a, long b)
    { zz_p x; div(x, a, b); return x; }

inline zz_p operator/(long a, zz_p b)
    { zz_p x; div(x, a, b); return x; }


inline zz_p& operator/=(zz_p& x, zz_p b)
   { div(x, x, b); return x; }

inline zz_p& operator/=(zz_p& x, long b)
   { div(x, x, b); return x; }


// ****** exponentiation

inline void power(zz_p& x, zz_p a, long e)
// x = a^e

   { x._zz_p__rep = PowerMod(a._zz_p__rep, e, zz_p::modulus()); }

inline zz_p power(zz_p a, long e)
   { zz_p x; power(x, a, e); return x; }

// ****** comparison

inline long IsZero(zz_p a)
   { return a._zz_p__rep == 0; }

inline long IsOne(zz_p a)
   { return a._zz_p__rep == 1; }

inline long operator==(zz_p a, zz_p b)
   { return a._zz_p__rep == b._zz_p__rep; }

inline long operator!=(zz_p a, zz_p b)
   { return !(a == b); }

inline long operator==(zz_p a, long b) { return a == to_zz_p(b); }
inline long operator==(long a, zz_p b) { return to_zz_p(a) == b; }

inline long operator!=(zz_p a, long b) { return !(a == b); }
inline long operator!=(long a, zz_p b) { return !(a == b); }

// ****** random numbers

inline void random(zz_p& x)
// x = random element in zz_p

   { x._zz_p__rep = RandomBnd(zz_p::modulus()); }

inline zz_p random_zz_p()
   { zz_p x; random(x); return x; }



// ****** input/output

NTL_SNS ostream& operator<<(NTL_SNS ostream& s, zz_p a);
   
NTL_SNS istream& operator>>(NTL_SNS istream& s, zz_p& x);


/* additional legacy conversions for v6 conversion regime */

inline void conv(int& x, zz_p a) { conv(x, rep(a)); }
inline void conv(unsigned int& x, zz_p a) { conv(x, rep(a)); }
inline void conv(long& x, zz_p a) { conv(x, rep(a)); }
inline void conv(unsigned long& x, zz_p a) { conv(x, rep(a)); }
inline void conv(ZZ& x, zz_p a) { conv(x, rep(a)); }


inline void conv(zz_p& x, zz_p a) { x = a; }

/* ------------------------------------- */



NTL_CLOSE_NNS

#endif
