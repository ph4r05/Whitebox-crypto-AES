
#ifndef NTL_mat_zz_p__H
#define NTL_mat_zz_p__H

#include <NTL/matrix.h>
#include <NTL/vec_vec_lzz_p.h>

NTL_OPEN_NNS

typedef Mat<zz_p> mat_zz_p;


void add(mat_zz_p& X, const mat_zz_p& A, const mat_zz_p& B); 
void sub(mat_zz_p& X, const mat_zz_p& A, const mat_zz_p& B); 
void negate(mat_zz_p& X, const mat_zz_p& A);
void mul(mat_zz_p& X, const mat_zz_p& A, const mat_zz_p& B); 
void mul(vec_zz_p& x, const mat_zz_p& A, const vec_zz_p& b); 
void mul(vec_zz_p& x, const vec_zz_p& a, const mat_zz_p& B); 

void mul(mat_zz_p& X, const mat_zz_p& A, zz_p b);
void mul(mat_zz_p& X, const mat_zz_p& A, long b);

inline void mul(mat_zz_p& X, zz_p a, const mat_zz_p& B)
   { mul(X, B, a); }

inline void mul(mat_zz_p& X, long a, const mat_zz_p& B)
   { mul(X, B, a); }


void ident(mat_zz_p& X, long n); 
inline mat_zz_p ident_mat_zz_p(long n)
   { mat_zz_p X; ident(X, n); NTL_OPT_RETURN(mat_zz_p, X); }

void determinant(zz_p& d, const mat_zz_p& A);
long IsIdent(const mat_zz_p& A, long n);
void transpose(mat_zz_p& X, const mat_zz_p& A);
void solve(zz_p& d, vec_zz_p& X,
           const mat_zz_p& A, const vec_zz_p& b);
void inv(zz_p& d, mat_zz_p& X, const mat_zz_p& A);

inline void sqr(mat_zz_p& X, const mat_zz_p& A)
   { mul(X, A, A); }

inline mat_zz_p sqr(const mat_zz_p& A)
   { mat_zz_p X; sqr(X, A); NTL_OPT_RETURN(mat_zz_p, X); }

void inv(mat_zz_p& X, const mat_zz_p& A);

inline mat_zz_p inv(const mat_zz_p& A)
   { mat_zz_p X; inv(X, A); NTL_OPT_RETURN(mat_zz_p, X); }

void power(mat_zz_p& X, const mat_zz_p& A, const ZZ& e);
inline mat_zz_p power(const mat_zz_p& A, const ZZ& e)
   { mat_zz_p X; power(X, A, e); NTL_OPT_RETURN(mat_zz_p, X); }

inline void power(mat_zz_p& X, const mat_zz_p& A, long e)
   { power(X, A, ZZ_expo(e)); }
inline mat_zz_p power(const mat_zz_p& A, long e)
   { mat_zz_p X; power(X, A, e); NTL_OPT_RETURN(mat_zz_p, X); }


void diag(mat_zz_p& X, long n, zz_p d);
inline mat_zz_p diag(long n, zz_p d)
   { mat_zz_p X; diag(X, n, d); NTL_OPT_RETURN(mat_zz_p, X); }

long IsDiag(const mat_zz_p& A, long n, zz_p d);


long gauss(mat_zz_p& M);
long gauss(mat_zz_p& M, long w);
void image(mat_zz_p& X, const mat_zz_p& A);
void kernel(mat_zz_p& X, const mat_zz_p& A);



// miscellaneous:

inline zz_p determinant(const mat_zz_p& a)
   { zz_p x; determinant(x, a); return x; }
// functional variant of determinant

inline mat_zz_p transpose(const mat_zz_p& a)
   { mat_zz_p x; transpose(x, a); NTL_OPT_RETURN(mat_zz_p, x); }

void clear(mat_zz_p& a);
// x = 0 (dimension unchanged)

long IsZero(const mat_zz_p& a);
// test if a is the zero matrix (any dimension)


// operator notation:

mat_zz_p operator+(const mat_zz_p& a, const mat_zz_p& b);
mat_zz_p operator-(const mat_zz_p& a, const mat_zz_p& b);
mat_zz_p operator*(const mat_zz_p& a, const mat_zz_p& b);

mat_zz_p operator-(const mat_zz_p& a);


// matrix/scalar multiplication:

inline mat_zz_p operator*(const mat_zz_p& a, zz_p b)
   { mat_zz_p x; mul(x, a, b); NTL_OPT_RETURN(mat_zz_p, x); }

inline mat_zz_p operator*(const mat_zz_p& a, long b)
   { mat_zz_p x; mul(x, a, b); NTL_OPT_RETURN(mat_zz_p, x); }

inline mat_zz_p operator*(zz_p a, const mat_zz_p& b)
   { mat_zz_p x; mul(x, a, b); NTL_OPT_RETURN(mat_zz_p, x); }

inline mat_zz_p operator*(long a, const mat_zz_p& b)
   { mat_zz_p x; mul(x, a, b); NTL_OPT_RETURN(mat_zz_p, x); }



// matrix/vector multiplication:

vec_zz_p operator*(const mat_zz_p& a, const vec_zz_p& b);

vec_zz_p operator*(const vec_zz_p& a, const mat_zz_p& b);




// assignment operator notation:

inline mat_zz_p& operator+=(mat_zz_p& x, const mat_zz_p& a)
{
   add(x, x, a);
   return x;
}   

inline mat_zz_p& operator-=(mat_zz_p& x, const mat_zz_p& a)
{
   sub(x, x, a);
   return x;
}   


inline mat_zz_p& operator*=(mat_zz_p& x, const mat_zz_p& a)
{
   mul(x, x, a);
   return x;
}   

inline mat_zz_p& operator*=(mat_zz_p& x, zz_p a)
{
   mul(x, x, a);
   return x;
}   

inline mat_zz_p& operator*=(mat_zz_p& x, long a)
{
   mul(x, x, a);
   return x;
}   
   

inline vec_zz_p& operator*=(vec_zz_p& x, const mat_zz_p& a)
{
   mul(x, x, a);
   return x;
}   



NTL_CLOSE_NNS


#endif
