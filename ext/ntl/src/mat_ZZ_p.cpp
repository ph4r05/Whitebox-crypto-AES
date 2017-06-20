
#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZVec.h>
#include <NTL/vec_long.h>

NTL_START_IMPL

  
void add(mat_ZZ_p& X, const mat_ZZ_p& A, const mat_ZZ_p& B)  
{  
   long n = A.NumRows();  
   long m = A.NumCols();  
  
   if (B.NumRows() != n || B.NumCols() != m)   
      Error("matrix add: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i, j;  
   for (i = 1; i <= n; i++)   
      for (j = 1; j <= m; j++)  
         add(X(i,j), A(i,j), B(i,j));  
}  
  
void sub(mat_ZZ_p& X, const mat_ZZ_p& A, const mat_ZZ_p& B)  
{  
   long n = A.NumRows();  
   long m = A.NumCols();  
  
   if (B.NumRows() != n || B.NumCols() != m)  
      Error("matrix sub: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i, j;  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= m; j++)  
         sub(X(i,j), A(i,j), B(i,j));  
}  

void negate(mat_ZZ_p& X, const mat_ZZ_p& A)  
{  
   long n = A.NumRows();  
   long m = A.NumCols();  
  
  
   X.SetDims(n, m);  
  
   long i, j;  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= m; j++)  
         negate(X(i,j), A(i,j));  
}  
  
void mul_aux(mat_ZZ_p& X, const mat_ZZ_p& A, const mat_ZZ_p& B)  
{  
   long n = A.NumRows();  
   long l = A.NumCols();  
   long m = B.NumCols();  
  
   if (l != B.NumRows())  
      Error("matrix mul: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i, j, k;  
   ZZ acc, tmp;  
  
   for (i = 1; i <= n; i++) {  
      for (j = 1; j <= m; j++) {  
         clear(acc);  
         for(k = 1; k <= l; k++) {  
            mul(tmp, rep(A(i,k)), rep(B(k,j)));  
            add(acc, acc, tmp);  
         }  
         conv(X(i,j), acc);  
      }  
   }  
}  
  
  
void mul(mat_ZZ_p& X, const mat_ZZ_p& A, const mat_ZZ_p& B)  
{  
   if (&X == &A || &X == &B) {  
      mat_ZZ_p tmp;  
      mul_aux(tmp, A, B);  
      X = tmp;  
   }  
   else  
      mul_aux(X, A, B);  
}  
  
  
static
void mul_aux(vec_ZZ_p& x, const mat_ZZ_p& A, const vec_ZZ_p& b)  
{  
   long n = A.NumRows();  
   long l = A.NumCols();  
  
   if (l != b.length())  
      Error("matrix mul: dimension mismatch");  
  
   x.SetLength(n);  
  
   long i, k;  
   ZZ acc, tmp;  
  
   for (i = 1; i <= n; i++) {  
      clear(acc);  
      for (k = 1; k <= l; k++) {  
         mul(tmp, rep(A(i,k)), rep(b(k)));  
         add(acc, acc, tmp);  
      }  
      conv(x(i), acc);  
   }  
}  
  
  
void mul(vec_ZZ_p& x, const mat_ZZ_p& A, const vec_ZZ_p& b)  
{  
   if (&b == &x || A.position1(x) != -1) {
      vec_ZZ_p tmp;
      mul_aux(tmp, A, b);
      x = tmp;
   }
   else
      mul_aux(x, A, b);
}  

static
void mul_aux(vec_ZZ_p& x, const vec_ZZ_p& a, const mat_ZZ_p& B)  
{  
   long n = B.NumRows();  
   long l = B.NumCols();  
  
   if (n != a.length())  
      Error("matrix mul: dimension mismatch");  
  
   x.SetLength(l);  
  
   long i, k;  
   ZZ acc, tmp;  
  
   for (i = 1; i <= l; i++) {  
      clear(acc);  
      for (k = 1; k <= n; k++) {  
         mul(tmp, rep(a(k)), rep(B(k,i)));
         add(acc, acc, tmp);  
      }  
      conv(x(i), acc);  
   }  
}  

void mul(vec_ZZ_p& x, const vec_ZZ_p& a, const mat_ZZ_p& B)
{
   if (&a == &x) {
      vec_ZZ_p tmp;
      mul_aux(tmp, a, B);
      x = tmp;
   }
   else
      mul_aux(x, a, B);
}

     
  
void ident(mat_ZZ_p& X, long n)  
{  
   X.SetDims(n, n);  
   long i, j;  
  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= n; j++)  
         if (i == j)  
            set(X(i, j));  
         else  
            clear(X(i, j));  
} 


void determinant(ZZ_p& d, const mat_ZZ_p& M_in)
{
   long k, n;
   long i, j;
   long pos;
   ZZ t1, t2;
   ZZ *x, *y;

   const ZZ& p = ZZ_p::modulus();

   n = M_in.NumRows();

   if (M_in.NumCols() != n)
      Error("determinant: nonsquare matrix");

   if (n == 0) {
      set(d);
      return;
   }

   vec_ZZVec M;
   sqr(t1, p);
   mul(t1, t1, n);

   M.SetLength(n);
   for (i = 0; i < n; i++) {
      M[i].SetSize(n, t1.size());
      for (j = 0; j < n; j++)
         M[i][j] = rep(M_in[i][j]);
   }

   ZZ det;
   set(det);

   for (k = 0; k < n; k++) {
      pos = -1;
      for (i = k; i < n; i++) {
         rem(t1, M[i][k], p);
         M[i][k] = t1;
         if (pos == -1 && !IsZero(t1))
            pos = i;
      }

      if (pos != -1) {
         if (k != pos) {
            swap(M[pos], M[k]);
            NegateMod(det, det, p);
         }

         MulMod(det, det, M[k][k], p);

         // make M[k, k] == -1 mod p, and make row k reduced

         InvMod(t1, M[k][k], p);
         NegateMod(t1, t1, p);
         for (j = k+1; j < n; j++) {
            rem(t2, M[k][j], p);
            MulMod(M[k][j], t2, t1, p);
         }

         for (i = k+1; i < n; i++) {
            // M[i] = M[i] + M[k]*M[i,k]

            t1 = M[i][k];   // this is already reduced

            x = M[i].elts() + (k+1);
            y = M[k].elts() + (k+1);

            for (j = k+1; j < n; j++, x++, y++) {
               // *x = *x + (*y)*t1

               mul(t2, *y, t1);
               add(*x, *x, t2);
            }
         }
      }
      else {
         clear(d);
         return;
      }
   }

   conv(d, det);
}

long IsIdent(const mat_ZZ_p& A, long n)
{
   if (A.NumRows() != n || A.NumCols() != n)
      return 0;

   long i, j;

   for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
         if (i != j) {
            if (!IsZero(A(i, j))) return 0;
         }
         else {
            if (!IsOne(A(i, j))) return 0;
         }

   return 1;
}
            

void transpose(mat_ZZ_p& X, const mat_ZZ_p& A)
{
   long n = A.NumRows();
   long m = A.NumCols();

   long i, j;

   if (&X == & A) {
      if (n == m)
         for (i = 1; i <= n; i++)
            for (j = i+1; j <= n; j++)
               swap(X(i, j), X(j, i));
      else {
         mat_ZZ_p tmp;
         tmp.SetDims(m, n);
         for (i = 1; i <= n; i++)
            for (j = 1; j <= m; j++)
               tmp(j, i) = A(i, j);
         X.kill();
         X = tmp;
      }
   }
   else {
      X.SetDims(m, n);
      for (i = 1; i <= n; i++)
         for (j = 1; j <= m; j++)
            X(j, i) = A(i, j);
   }
}
   

void solve(ZZ_p& d, vec_ZZ_p& X, 
           const mat_ZZ_p& A, const vec_ZZ_p& b)

{
   long n = A.NumRows();
   if (A.NumCols() != n)
      Error("solve: nonsquare matrix");

   if (b.length() != n)
      Error("solve: dimension mismatch");

   if (n == 0) {
      set(d);
      X.SetLength(0);
      return;
   }

   long i, j, k, pos;
   ZZ t1, t2;
   ZZ *x, *y;

   const ZZ& p = ZZ_p::modulus();

   vec_ZZVec M;
   sqr(t1, p);
   mul(t1, t1, n);

   M.SetLength(n);

   for (i = 0; i < n; i++) {
      M[i].SetSize(n+1, t1.size());
      for (j = 0; j < n; j++) 
         M[i][j] = rep(A[j][i]);
      M[i][n] = rep(b[i]);
   }

   ZZ det;
   set(det);

   for (k = 0; k < n; k++) {
      pos = -1;
      for (i = k; i < n; i++) {
         rem(t1, M[i][k], p);
         M[i][k] = t1;
         if (pos == -1 && !IsZero(t1)) {
            pos = i;
         }
      }

      if (pos != -1) {
         if (k != pos) {
            swap(M[pos], M[k]);
            NegateMod(det, det, p);
         }

         MulMod(det, det, M[k][k], p);

         // make M[k, k] == -1 mod p, and make row k reduced

         InvMod(t1, M[k][k], p);
         NegateMod(t1, t1, p);
         for (j = k+1; j <= n; j++) {
            rem(t2, M[k][j], p);
            MulMod(M[k][j], t2, t1, p);
         }

         for (i = k+1; i < n; i++) {
            // M[i] = M[i] + M[k]*M[i,k]

            t1 = M[i][k];   // this is already reduced

            x = M[i].elts() + (k+1);
            y = M[k].elts() + (k+1);

            for (j = k+1; j <= n; j++, x++, y++) {
               // *x = *x + (*y)*t1

               mul(t2, *y, t1);
               add(*x, *x, t2);
            }
         }
      }
      else {
         clear(d);
         return;
      }
   }

   X.SetLength(n);
   for (i = n-1; i >= 0; i--) {
      clear(t1);
      for (j = i+1; j < n; j++) {
         mul(t2, rep(X[j]), M[i][j]);
         add(t1, t1, t2);
      }
      sub(t1, t1, M[i][n]);
      conv(X[i], t1);
   }

   conv(d, det);
}

void inv(ZZ_p& d, mat_ZZ_p& X, const mat_ZZ_p& A)
{
   long n = A.NumRows();
   if (A.NumCols() != n)
      Error("inv: nonsquare matrix");

   if (n == 0) {
      set(d);
      X.SetDims(0, 0);
      return;
   }

   long i, j, k, pos;
   ZZ t1, t2;
   ZZ *x, *y;

   const ZZ& p = ZZ_p::modulus();

   vec_ZZVec M;
   sqr(t1, p);
   mul(t1, t1, n);

   M.SetLength(n);

   for (i = 0; i < n; i++) {
      M[i].SetSize(2*n, t1.size());
      for (j = 0; j < n; j++) {
         M[i][j] = rep(A[i][j]);
         clear(M[i][n+j]);
      }
      set(M[i][n+i]);
   }

   ZZ det;
   set(det);

   for (k = 0; k < n; k++) {
      pos = -1;
      for (i = k; i < n; i++) {
         rem(t1, M[i][k], p);
         M[i][k] = t1;
         if (pos == -1 && !IsZero(t1)) {
            pos = i;
         }
      }

      if (pos != -1) {
         if (k != pos) {
            swap(M[pos], M[k]);
            NegateMod(det, det, p);
         }

         MulMod(det, det, M[k][k], p);

         // make M[k, k] == -1 mod p, and make row k reduced

         InvMod(t1, M[k][k], p);
         NegateMod(t1, t1, p);
         for (j = k+1; j < 2*n; j++) {
            rem(t2, M[k][j], p);
            MulMod(M[k][j], t2, t1, p);
         }

         for (i = k+1; i < n; i++) {
            // M[i] = M[i] + M[k]*M[i,k]

            t1 = M[i][k];   // this is already reduced

            x = M[i].elts() + (k+1);
            y = M[k].elts() + (k+1);

            for (j = k+1; j < 2*n; j++, x++, y++) {
               // *x = *x + (*y)*t1

               mul(t2, *y, t1);
               add(*x, *x, t2);
            }
         }
      }
      else {
         clear(d);
         return;
      }
   }

   X.SetDims(n, n);
   for (k = 0; k < n; k++) {
      for (i = n-1; i >= 0; i--) {
         clear(t1);
         for (j = i+1; j < n; j++) {
            mul(t2, rep(X[j][k]), M[i][j]);
            add(t1, t1, t2);
         }
         sub(t1, t1, M[i][n+k]);
         conv(X[i][k], t1);
      }
   }

   conv(d, det);
}



long gauss(mat_ZZ_p& M_in, long w)
{
   long k, l;
   long i, j;
   long pos;
   ZZ t1, t2, t3;
   ZZ *x, *y;

   long n = M_in.NumRows();
   long m = M_in.NumCols();

   if (w < 0 || w > m)
      Error("gauss: bad args");

   const ZZ& p = ZZ_p::modulus();

   vec_ZZVec M;
   sqr(t1, p);
   mul(t1, t1, n);

   M.SetLength(n);
   for (i = 0; i < n; i++) {
      M[i].SetSize(m, t1.size());
      for (j = 0; j < m; j++) {
         M[i][j] = rep(M_in[i][j]);
      }
   }

   l = 0;
   for (k = 0; k < w && l < n; k++) {

      pos = -1;
      for (i = l; i < n; i++) {
         rem(t1, M[i][k], p);
         M[i][k] = t1;
         if (pos == -1 && !IsZero(t1)) {
            pos = i;
         }
      }

      if (pos != -1) {
         swap(M[pos], M[l]);

         InvMod(t3, M[l][k], p);
         NegateMod(t3, t3, p);

         for (j = k+1; j < m; j++) {
            rem(M[l][j], M[l][j], p);
         }

         for (i = l+1; i < n; i++) {
            // M[i] = M[i] + M[l]*M[i,k]*t3

            MulMod(t1, M[i][k], t3, p);

            clear(M[i][k]);

            x = M[i].elts() + (k+1);
            y = M[l].elts() + (k+1);

            for (j = k+1; j < m; j++, x++, y++) {
               // *x = *x + (*y)*t1

               mul(t2, *y, t1);
               add(t2, t2, *x);
               *x = t2;
            }
         }

         l++;
      }
   }
   
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         conv(M_in[i][j], M[i][j]);

   return l;
}

long gauss(mat_ZZ_p& M)
{
   return gauss(M, M.NumCols());
}

void image(mat_ZZ_p& X, const mat_ZZ_p& A)
{
   mat_ZZ_p M;
   M = A;
   long r = gauss(M);
   M.SetDims(r, M.NumCols());
   X = M;
}

void kernel(mat_ZZ_p& X, const mat_ZZ_p& A)
{
   long m = A.NumRows();
   long n = A.NumCols();

   mat_ZZ_p M;
   long r;

   transpose(M, A);
   r = gauss(M);

   X.SetDims(m-r, m);

   long i, j, k, s;
   ZZ t1, t2;

   ZZ_p T3;

   vec_long D;
   D.SetLength(m);
   for (j = 0; j < m; j++) D[j] = -1;

   vec_ZZ_p inverses;
   inverses.SetLength(m);

   j = -1;
   for (i = 0; i < r; i++) {
      do {
         j++;
      } while (IsZero(M[i][j]));

      D[j] = i;
      inv(inverses[j], M[i][j]); 
   }

   for (k = 0; k < m-r; k++) {
      vec_ZZ_p& v = X[k];
      long pos = 0;
      for (j = m-1; j >= 0; j--) {
         if (D[j] == -1) {
            if (pos == k)
               set(v[j]);
            else
               clear(v[j]);
            pos++;
         }
         else {
            i = D[j];

            clear(t1);

            for (s = j+1; s < m; s++) {
               mul(t2, rep(v[s]), rep(M[i][s]));
               add(t1, t1, t2);
            }

            conv(T3, t1);
            mul(T3, T3, inverses[j]);
            negate(v[j], T3);
         }
      }
   }
}
   
void mul(mat_ZZ_p& X, const mat_ZZ_p& A, const ZZ_p& b_in)
{
   NTL_ZZ_pRegister(b);
   b = b_in;
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);

   long i, j;
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         mul(X[i][j], A[i][j], b);
}
   
void mul(mat_ZZ_p& X, const mat_ZZ_p& A, long b_in)
{
   NTL_ZZ_pRegister(b);
   b = b_in;
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);

   long i, j;
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         mul(X[i][j], A[i][j], b);
}

void diag(mat_ZZ_p& X, long n, const ZZ_p& d_in)  
{  
   ZZ_p d = d_in;
   X.SetDims(n, n);  
   long i, j;  
  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= n; j++)  
         if (i == j)  
            X(i, j) = d;  
         else  
            clear(X(i, j));  
} 

long IsDiag(const mat_ZZ_p& A, long n, const ZZ_p& d)
{
   if (A.NumRows() != n || A.NumCols() != n)
      return 0;

   long i, j;

   for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++)
         if (i != j) {
            if (!IsZero(A(i, j))) return 0;
         }
         else {
            if (A(i, j) != d) return 0;
         }

   return 1;
}


long IsZero(const mat_ZZ_p& a)
{
   long n = a.NumRows();
   long i;

   for (i = 0; i < n; i++)
      if (!IsZero(a[i]))
         return 0;

   return 1;
}

void clear(mat_ZZ_p& x)
{
   long n = x.NumRows();
   long i;
   for (i = 0; i < n; i++)
      clear(x[i]);
}


mat_ZZ_p operator+(const mat_ZZ_p& a, const mat_ZZ_p& b)
{
   mat_ZZ_p res;
   add(res, a, b);
   NTL_OPT_RETURN(mat_ZZ_p, res);
}

mat_ZZ_p operator*(const mat_ZZ_p& a, const mat_ZZ_p& b)
{
   mat_ZZ_p res;
   mul_aux(res, a, b);
   NTL_OPT_RETURN(mat_ZZ_p, res);
}

mat_ZZ_p operator-(const mat_ZZ_p& a, const mat_ZZ_p& b)
{
   mat_ZZ_p res;
   sub(res, a, b);
   NTL_OPT_RETURN(mat_ZZ_p, res);
}


mat_ZZ_p operator-(const mat_ZZ_p& a)
{
   mat_ZZ_p res;
   negate(res, a);
   NTL_OPT_RETURN(mat_ZZ_p, res);
}


vec_ZZ_p operator*(const mat_ZZ_p& a, const vec_ZZ_p& b)
{
   vec_ZZ_p res;
   mul_aux(res, a, b);
   NTL_OPT_RETURN(vec_ZZ_p, res);
}

vec_ZZ_p operator*(const vec_ZZ_p& a, const mat_ZZ_p& b)
{
   vec_ZZ_p res;
   mul_aux(res, a, b);
   NTL_OPT_RETURN(vec_ZZ_p, res);
}

void inv(mat_ZZ_p& X, const mat_ZZ_p& A)
{
   ZZ_p d;
   inv(d, X, A);
   if (d == 0) Error("inv: non-invertible matrix");
}

void power(mat_ZZ_p& X, const mat_ZZ_p& A, const ZZ& e)
{
   if (A.NumRows() != A.NumCols()) Error("power: non-square matrix");

   if (e == 0) {
      ident(X, A.NumRows());
      return;
   }

   mat_ZZ_p T1, T2;
   long i, k;

   k = NumBits(e);
   T1 = A;

   for (i = k-2; i >= 0; i--) {
      sqr(T2, T1);
      if (bit(e, i))
         mul(T1, T2, A);
      else
         T1 = T2;
   }

   if (e < 0)
      inv(X, T1);
   else
      X = T1;
}

NTL_END_IMPL
