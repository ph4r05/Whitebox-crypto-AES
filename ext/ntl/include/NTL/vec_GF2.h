
#ifndef NTL_vec_GF2__H
#define NTL_vec_GF2__H

#include <NTL/WordVector.h>
#include <NTL/GF2.h>

NTL_OPEN_NNS


// Vec<GF2> is an explicit specialization of Vec<T>.
// Vec<GF2> is declared, but not defined, in GF2.h,
// to prevent the generic Vec from being used.

template<> 
class Vec<GF2> {

public:

// these should be private, but they are not

   WordVector rep;

   long _len;  // length (in bits)
   long _maxlen;  // (MaxLength << 1) | (fixed)

   // invariants: rep.length() "tracks" length() ( = _len)
   //             All bits in positions >= length are zero.

   // Note:       rep.MaxLength() may exceed the value
   //             indicated by MaxLength().
   

//the following are "really" public


   Vec() : _len(0), _maxlen(0) {}
   Vec(INIT_SIZE_TYPE, long n) : _len(0), _maxlen(0) { SetLength(n); }
   Vec(const Vec<GF2>& a) : _len(0), _maxlen(0) { *this = a; }

   Vec& operator=(const Vec<GF2>& a);

   ~Vec() {}

   void kill();

   void SetLength(long n);
   void SetMaxLength(long n);
   void FixLength(long n);

   long length() const { return _len; }
   long MaxLength() const { return _maxlen >> 1; }  
   long allocated() const { return rep.MaxLength() * NTL_BITS_PER_LONG; }
   long fixed() const { return _maxlen & 1; }


   Vec(Vec<GF2>& x, INIT_TRANS_TYPE) : 
      rep(x.rep, INIT_TRANS), _len(x._len), _maxlen(x._maxlen) { }

   const GF2 get(long i) const;
   void put(long i, GF2 a);
   void put(long i, long a) { put(i, to_GF2(a)); }

   ref_GF2 operator[](long i); 

   ref_GF2 operator()(long i) 
      { return (*this)[i-1]; }

   const GF2 operator[](long i) const 
      { return get(i); }

   const GF2 operator()(long i) const 
      { return get(i-1); }



// Some partial STL compatibility...also used
// to interface with the Matrix template class

   typedef GF2 value_type;
   typedef ref_GF2 reference;
   typedef const GF2 const_reference;

};

typedef Vec<GF2> vec_GF2;


// sepcialized conversion

inline void conv(vec_GF2& x, const vec_GF2& a)
{  x = a; }


void swap(vec_GF2& x, vec_GF2& y);
void append(vec_GF2& v, const GF2& a);
void append(vec_GF2& v, const vec_GF2& a);

long operator==(const vec_GF2& a, const vec_GF2& b);
inline long operator!=(const vec_GF2& a, const vec_GF2& b)
   { return !(a == b); }

NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const vec_GF2& a);
NTL_SNS istream& operator>>(NTL_SNS istream& s, vec_GF2& a);

void shift(vec_GF2& x, const vec_GF2& a, long n);
// x = a shifted n places, i.e., if l = a.length(),
//    x.length() = l, x[i] = a[i-n] for 0 <= i-n < l,
//    and x[i] = 0 for all other i such that 0 <= i < l.

inline vec_GF2 shift(const vec_GF2& a, long n)
   { vec_GF2 x; shift(x, a, n); NTL_OPT_RETURN(vec_GF2, x); }

void reverse(vec_GF2& x, const vec_GF2& a);

inline vec_GF2 reverse(const vec_GF2& a)
   { vec_GF2 x; reverse(x, a); NTL_OPT_RETURN(vec_GF2, x); }

void random(vec_GF2& x, long n);
inline vec_GF2 random_vec_GF2(long n)
   { vec_GF2 x; random(x, n); NTL_OPT_RETURN(vec_GF2, x); }

long weight(const vec_GF2& a);

void mul(vec_GF2& x, const vec_GF2& a, GF2 b);
inline void mul(vec_GF2& x, GF2 a, const vec_GF2& b)
   { mul(x, b, a); }

inline void mul(vec_GF2& x, const vec_GF2& a, long b)
   { mul(x, a, to_GF2(b)); }
inline void mul(vec_GF2& x, long a, const vec_GF2& b)
   { mul(x, b, a); }

void add(vec_GF2& x, const vec_GF2& a, const vec_GF2& b);

inline void sub(vec_GF2& x, const vec_GF2& a, const vec_GF2& b)
   { add(x, a, b); }

void clear(vec_GF2& x);

inline void negate(vec_GF2& x, const vec_GF2& a)
   { x = a; }

inline void InnerProduct(ref_GF2 x, const vec_GF2& a, const vec_GF2& b)
   { x = to_GF2(InnerProduct(a.rep, b.rep)); }

long IsZero(const vec_GF2& a);

vec_GF2 operator+(const vec_GF2& a, const vec_GF2& b);

vec_GF2 operator-(const vec_GF2& a, const vec_GF2& b);

inline vec_GF2 operator-(const vec_GF2& a)
   { return a; }

inline vec_GF2 operator*(const vec_GF2& a, GF2 b)
   { vec_GF2 x; mul(x, a, b); NTL_OPT_RETURN(vec_GF2, x); }

inline vec_GF2 operator*(const vec_GF2& a, long b)
   { vec_GF2 x; mul(x, a, b); NTL_OPT_RETURN(vec_GF2, x); }

inline vec_GF2 operator*(GF2 a, const vec_GF2& b)
   { vec_GF2 x; mul(x, a, b); NTL_OPT_RETURN(vec_GF2, x); }

inline vec_GF2 operator*(long a, const vec_GF2& b)
   { vec_GF2 x; mul(x, a, b); NTL_OPT_RETURN(vec_GF2, x); }


inline GF2 operator*(const vec_GF2& a, const vec_GF2& b)
   { return to_GF2(InnerProduct(a.rep, b.rep)); }

// assignment operator notation:

inline vec_GF2& operator+=(vec_GF2& x, const vec_GF2& a)
{ 
   add(x, x, a);
   return x;
}

inline vec_GF2& operator-=(vec_GF2& x, const vec_GF2& a)
{ 
   sub(x, x, a);
   return x;
}

inline vec_GF2& operator*=(vec_GF2& x, GF2 a)
{ 
   mul(x, x, a);
   return x;
}

inline vec_GF2& operator*=(vec_GF2& x, long a)
{ 
   mul(x, x, a);
   return x;
}

void VectorCopy(vec_GF2& x, const vec_GF2& a, long n);
inline vec_GF2 VectorCopy(const vec_GF2& a, long n)
   { vec_GF2 x; VectorCopy(x, a, n); NTL_OPT_RETURN(vec_GF2, x); }

NTL_CLOSE_NNS


#endif


