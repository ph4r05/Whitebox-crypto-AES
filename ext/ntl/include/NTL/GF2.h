
#ifndef NTL_GF2__H
#define NTL_GF2__H

#include <NTL/ZZ.h>
#include <NTL/vector.h>

NTL_OPEN_NNS

class GF2 {
public:

unsigned long _GF2__rep;


GF2() : _GF2__rep(0) { }
GF2(const GF2& a) : _GF2__rep(a._GF2__rep) { }

GF2(INIT_VAL_TYPE, long a) : _GF2__rep(a & 1) { }
GF2(INIT_LOOP_HOLE_TYPE, unsigned long a) : _GF2__rep(a) { }

~GF2() { }

GF2& operator=(const GF2& a) { _GF2__rep = a._GF2__rep; return *this; }
GF2& operator=(long a) { _GF2__rep = a & 1; return *this; }

static long modulus() { return 2; }
static GF2 zero() { return GF2(); }

};



class ref_GF2 {
public:

unsigned long *_ref_GF2__ptr;
long _ref_GF2__pos;

ref_GF2() : _ref_GF2__ptr(0), _ref_GF2__pos(0) { }
ref_GF2(const ref_GF2& a) : 
   _ref_GF2__ptr(a._ref_GF2__ptr), _ref_GF2__pos(a._ref_GF2__pos) { }
ref_GF2(GF2& a) :
   _ref_GF2__ptr(&a._GF2__rep),  _ref_GF2__pos(0) { }
ref_GF2(INIT_LOOP_HOLE_TYPE, unsigned long *ptr, long pos) :
   _ref_GF2__ptr(ptr), _ref_GF2__pos(pos) { }

operator const GF2() const 
{
   return GF2(INIT_LOOP_HOLE, (*_ref_GF2__ptr >> _ref_GF2__pos) & 1);
}

~ref_GF2() { }

ref_GF2 operator=(const ref_GF2& a)
{
   unsigned long rval = (*a._ref_GF2__ptr >> a._ref_GF2__pos) & 1;
   unsigned long lval = *_ref_GF2__ptr;
   lval = (lval & ~(1UL << _ref_GF2__pos)) | (rval << _ref_GF2__pos);
   *_ref_GF2__ptr = lval;
   return *this;
}

ref_GF2 operator=(const GF2& a)
{
   unsigned long rval = (a._GF2__rep) & 1;
   unsigned long lval = *_ref_GF2__ptr;
   lval = (lval & ~(1UL << _ref_GF2__pos)) | (rval << _ref_GF2__pos);
   *_ref_GF2__ptr = lval;
   return *this;
}


ref_GF2 operator=(long a)
{
   unsigned long rval = a & 1;
   unsigned long lval = *_ref_GF2__ptr;
   lval = (lval & ~(1UL << _ref_GF2__pos)) | (rval << _ref_GF2__pos);
   *_ref_GF2__ptr = lval;
   return *this;
}


};




// functions


inline long rep(GF2 a) { return a._GF2__rep; }



inline long IsZero(GF2 a)
   { return a._GF2__rep == 0; }

inline long IsOne(GF2 a)
   { return a._GF2__rep == 1; }




inline GF2 to_GF2(long a) 
   { return GF2(INIT_VAL, a); }

inline GF2 to_GF2(const ZZ& a) 
   { return GF2(INIT_LOOP_HOLE, IsOdd(a)); }




inline GF2 operator+(GF2 a, GF2 b)
   { return GF2(INIT_LOOP_HOLE, a._GF2__rep ^ b._GF2__rep); } 

inline GF2 operator+(GF2 a, long b)
   { return a + to_GF2(b); }

inline GF2 operator+(long a, GF2 b)
   { return to_GF2(a) + b; }

inline GF2 operator-(GF2 a, GF2 b)
   { return a + b; }

inline GF2 operator-(GF2 a, long b)
   { return a + b; }

inline GF2 operator-(long a, GF2 b)
   { return a + b; }

inline GF2 operator-(GF2 a)
   { return a; }


inline GF2 sqr(GF2 a)
   { return a; }

inline GF2 operator*(GF2 a, GF2 b)
   { return GF2(INIT_LOOP_HOLE,  a._GF2__rep & b._GF2__rep); }

inline GF2 operator*(GF2 a, long b)
   { return a * to_GF2(b); } 

inline GF2 operator*(long a, GF2 b)
   { return to_GF2(a) * b; }






inline GF2 operator/(GF2 a, GF2 b)
{
   if (IsZero(b)) Error("GF2: division by zero");
   return a;
}

inline GF2 operator/(GF2 a, long b)
   { return a / to_GF2(b); }

inline GF2 operator/(long a, GF2 b)
   { return to_GF2(a) / b; }


 
inline GF2 inv(GF2 a)
   { return 1 / a; }







inline long operator==(GF2 a, GF2 b)
   { return a._GF2__rep == b._GF2__rep; }


inline long operator==(GF2 a, long b) 
   { return a == to_GF2(b); }

inline long operator==(long a, GF2 b) 
   { return to_GF2(a) == b; }

inline long operator!=(GF2 a, GF2 b)  { return !(a == b); }
inline long operator!=(GF2 a, long b) { return !(a == b); }
inline long operator!=(long a, GF2 b) { return !(a == b); }



GF2 power(GF2 a, long e);



inline GF2 random_GF2()
   { return GF2(INIT_LOOP_HOLE, RandomBnd(2)); }



// procedural versions

inline GF2& operator+=(GF2& x, GF2 b)
   { return x = x + b; }

inline GF2& operator+=(GF2& x, long b)
   { return x = x + b; }

inline GF2& operator-=(GF2& x, GF2 b)
   { return x = x - b; }

inline GF2& operator-=(GF2& x, long b)
   { return x = x - b; }

inline GF2& operator++(GF2& x) { return x = x + 1; }
inline void operator++(GF2& x, int) { x = x + 1; }
inline GF2& operator--(GF2& x) { return x = x - 1; }
inline void operator--(GF2& x, int) { x = x - 1; }

inline GF2& operator*=(GF2& x, GF2 b)
   { return x = x * b; }

inline GF2& operator*=(GF2& x, long b)
   { return x = x * b; }

inline GF2& operator/=(GF2& x, GF2 b)
   { return x = x / b; }

inline GF2& operator/=(GF2& x, long b)
   { return x = x / b; }




inline void conv(GF2& x, long a) { x = to_GF2(a); }

inline void conv(GF2& x, const ZZ& a) { x = to_GF2(a); }


inline void clear(GF2& x) { x = 0; }

inline void set(GF2& x) { x = 1; }

inline void swap(GF2& x, GF2& y)
   { GF2 t; t = x; x = y; y = t; }

inline void add(GF2& x, GF2 a, GF2 b)
   { x = a + b; }

inline void sub(GF2& x, GF2 a, GF2 b)
   { x = a - b; }

inline void negate(GF2& x, GF2 a)
   { x = -a; }

inline void add(GF2& x, GF2 a, long b)
   { x = a + b; }

inline void add(GF2& x, long a, GF2 b)
   { x = a + b; }

inline void sub(GF2& x, GF2 a, long b)
   { x = a - b; }

inline void sub(GF2& x, long a, GF2 b)
   { x = a - b; }


inline void mul(GF2& x, GF2 a, GF2 b)
   { x = a * b; }

inline void mul(GF2& x, GF2 a, long b)
   { x = a * b; }

inline void mul(GF2& x, long a, GF2 b)
   { x = a * b; }

inline void sqr(GF2& x, GF2 a)
   { x = sqr(a); }


inline void div(GF2& x, GF2 a, GF2 b)   
   { x = a / b; }

inline void div(GF2& x, long a, GF2 b)
   { x = a / b; }

inline void div(GF2& x, GF2 a, long b)
   { x = a / b; }

inline void inv(GF2& x, GF2 a)
   { x = inv(a); }


inline void power(GF2& x, GF2 a, long e)
   { x = power(a, e); }



inline void random(GF2& x)
   { x = random_GF2(); }

// ref_GF2 variants...theoretically, these would
// have sufficed, because of the implicit conversion
// from GF2& to ref_GF2, but it may be a bit more efficient
// to explicitly overload everything.  Moreover, 
// the return types of the += type operators would
// not be right.


inline ref_GF2 operator+=(ref_GF2 x, GF2 b)
   { return x = x + b; }

inline ref_GF2 operator+=(ref_GF2 x, long b)
   { return x = x + b; }

inline ref_GF2 operator-=(ref_GF2 x, GF2 b)
   { return x = x - b; }

inline ref_GF2 operator-=(ref_GF2 x, long b)
   { return x = x - b; }

inline ref_GF2 operator++(ref_GF2 x) { return x = x + 1; }
inline void operator++(ref_GF2 x, int) { x = x + 1; }
inline ref_GF2 operator--(ref_GF2 x) { return x = x - 1; }
inline void operator--(ref_GF2 x, int) { x = x - 1; }

inline ref_GF2 operator*=(ref_GF2 x, GF2 b)
   { return x = x * b; }

inline ref_GF2 operator*=(ref_GF2 x, long b)
   { return x = x * b; }

inline ref_GF2 operator/=(ref_GF2 x, GF2 b)
   { return x = x / b; }

inline ref_GF2 operator/=(ref_GF2 x, long b)
   { return x = x / b; }




inline void conv(ref_GF2 x, long a) { x = to_GF2(a); }

inline void conv(ref_GF2 x, const ZZ& a) { x = to_GF2(a); }


inline void clear(ref_GF2 x) { x = 0; }

inline void set(ref_GF2 x) { x = 1; }

inline void swap(ref_GF2 x, ref_GF2 y)
   { GF2 t; t = x; x = y; y = t; }

inline void add(ref_GF2 x, GF2 a, GF2 b)
   { x = a + b; }

inline void sub(ref_GF2 x, GF2 a, GF2 b)
   { x = a - b; }

inline void negate(ref_GF2 x, GF2 a)
   { x = -a; }

inline void add(ref_GF2 x, GF2 a, long b)
   { x = a + b; }

inline void add(ref_GF2 x, long a, GF2 b)
   { x = a + b; }

inline void sub(ref_GF2 x, GF2 a, long b)
   { x = a - b; }

inline void sub(ref_GF2 x, long a, GF2 b)
   { x = a - b; }


inline void mul(ref_GF2 x, GF2 a, GF2 b)
   { x = a * b; }

inline void mul(ref_GF2 x, GF2 a, long b)
   { x = a * b; }

inline void mul(ref_GF2 x, long a, GF2 b)
   { x = a * b; }

inline void sqr(ref_GF2 x, GF2 a)
   { x = sqr(a); }


inline void div(ref_GF2 x, GF2 a, GF2 b)   
   { x = a / b; }

inline void div(ref_GF2 x, long a, GF2 b)
   { x = a / b; }

inline void div(ref_GF2 x, GF2 a, long b)
   { x = a / b; }

inline void inv(ref_GF2 x, GF2 a)
   { x = inv(a); }


inline void power(ref_GF2 x, GF2 a, long e)
   { x = power(a, e); }



inline void random(ref_GF2 x)
   { x = random_GF2(); }




// I/O...for input, we only provide the ref_GF2 variant 

NTL_SNS ostream& operator<<(NTL_SNS ostream& s, GF2 a);

NTL_SNS istream& operator>>(NTL_SNS istream& s, ref_GF2 x);

/* additional legacy conversions for v6 conversion regime */

inline void conv(int& x, GF2 a) { conv(x, rep(a)); }
inline void conv(unsigned int& x, GF2 a) { conv(x, rep(a)); }
inline void conv(long& x, GF2 a) { conv(x, rep(a)); }
inline void conv(unsigned long& x, GF2 a) { conv(x, rep(a)); }
inline void conv(ZZ& x, GF2 a) { conv(x, rep(a)); }


inline void conv(GF2& x, GF2 a) { x = a; }
inline void conv(ref_GF2 x, GF2 a) { x = a; }

/* ------------------------------------- */


// Finally, we declare an specialization Vec<GF2>:

template<> class Vec<GF2>;

NTL_CLOSE_NNS

#endif

