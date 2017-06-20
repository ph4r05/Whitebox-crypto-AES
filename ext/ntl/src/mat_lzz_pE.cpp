
#include <NTL/mat_lzz_pE.h>

#include <NTL/new.h>

NTL_START_IMPL

  
void add(mat_zz_pE& X, const mat_zz_pE& A, const mat_zz_pE& B)  
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
  
void sub(mat_zz_pE& X, const mat_zz_pE& A, const mat_zz_pE& B)  
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

void negate(mat_zz_pE& X, const mat_zz_pE& A)  
{  
   long n = A.NumRows();  
   long m = A.NumCols();  
  
  
   X.SetDims(n, m);  
  
   long i, j;  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= m; j++)  
         negate(X(i,j), A(i,j));  
}  
  
void mul_aux(mat_zz_pE& X, const mat_zz_pE& A, const mat_zz_pE& B)  
{  
   long n = A.NumRows();  
   long l = A.NumCols();  
   long m = B.NumCols();  
  
   if (l != B.NumRows())  
      Error("matrix mul: dimension mismatch");  
  
   X.SetDims(n, m);  
  
   long i, j, k;  
   zz_pX acc, tmp;  
  
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
  
  
void mul(mat_zz_pE& X, const mat_zz_pE& A, const mat_zz_pE& B)  
{  
   if (&X == &A || &X == &B) {  
      mat_zz_pE tmp;  
      mul_aux(tmp, A, B);  
      X = tmp;  
   }  
   else  
      mul_aux(X, A, B);  
}  
  
  
static
void mul_aux(vec_zz_pE& x, const mat_zz_pE& A, const vec_zz_pE& b)  
{  
   long n = A.NumRows();  
   long l = A.NumCols();  
  
   if (l != b.length())  
      Error("matrix mul: dimension mismatch");  
  
   x.SetLength(n);  
  
   long i, k;  
   zz_pX acc, tmp;  
  
   for (i = 1; i <= n; i++) {  
      clear(acc);  
      for (k = 1; k <= l; k++) {  
         mul(tmp, rep(A(i,k)), rep(b(k)));  
         add(acc, acc, tmp);  
      }  
      conv(x(i), acc);  
   }  
}  
  
  
void mul(vec_zz_pE& x, const mat_zz_pE& A, const vec_zz_pE& b)  
{  
   if (&b == &x || A.position1(x) != -1) {
      vec_zz_pE tmp;
      mul_aux(tmp, A, b);
      x = tmp;
   }
   else
      mul_aux(x, A, b);
}  

static
void mul_aux(vec_zz_pE& x, const vec_zz_pE& a, const mat_zz_pE& B)  
{  
   long n = B.NumRows();  
   long l = B.NumCols();  
  
   if (n != a.length())  
      Error("matrix mul: dimension mismatch");  
  
   x.SetLength(l);  
  
   long i, k;  
   zz_pX acc, tmp;  
  
   for (i = 1; i <= l; i++) {  
      clear(acc);  
      for (k = 1; k <= n; k++) {  
         mul(tmp, rep(a(k)), rep(B(k,i)));
         add(acc, acc, tmp);  
      }  
      conv(x(i), acc);  
   }  
}  

void mul(vec_zz_pE& x, const vec_zz_pE& a, const mat_zz_pE& B)
{
   if (&a == &x) {
      vec_zz_pE tmp;
      mul_aux(tmp, a, B);
      x = tmp;
   }
   else
      mul_aux(x, a, B);

}

     
  
void ident(mat_zz_pE& X, long n)  
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


void determinant(zz_pE& d, const mat_zz_pE& M_in)
{
   long k, n;
   long i, j;
   long pos;
   zz_pX t1, t2;
   zz_pX *x, *y;

   const zz_pXModulus& p = zz_pE::modulus();

   n = M_in.NumRows();

   if (M_in.NumCols() != n)
      Error("determinant: nonsquare matrix");

   if (n == 0) {
      set(d);
      return;
   }

   vec_zz_pX *M = NTL_NEW_OP vec_zz_pX[n];

   for (i = 0; i < n; i++) {
      M[i].SetLength(n);
      for (j = 0; j < n; j++) {
         M[i][j].rep.SetMaxLength(2*deg(p)-1);
         M[i][j] = rep(M_in[i][j]);
      }
   }

   zz_pX det;
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
            negate(det, det);
         }

         MulMod(det, det, M[k][k], p);

         // make M[k, k] == -1 mod p, and make row k reduced

         InvMod(t1, M[k][k], p);
         negate(t1, t1);
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
         goto done;
      }
   }

   conv(d, det);

done:
   delete[] M;
}

long IsIdent(const mat_zz_pE& A, long n)
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
            

void transpose(mat_zz_pE& X, const mat_zz_pE& A)
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
         mat_zz_pE tmp;
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
   

void solve(zz_pE& d, vec_zz_pE& X, 
           const mat_zz_pE& A, const vec_zz_pE& b)

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
   zz_pX t1, t2;
   zz_pX *x, *y;

   const zz_pXModulus& p = zz_pE::modulus();

   vec_zz_pX *M = NTL_NEW_OP vec_zz_pX[n];

   for (i = 0; i < n; i++) {
      M[i].SetLength(n+1);
      for (j = 0; j < n; j++) {
         M[i][j].rep.SetMaxLength(2*deg(p)-1);
         M[i][j] = rep(A[j][i]);
      }
      M[i][n].rep.SetMaxLength(2*deg(p)-1);
      M[i][n] = rep(b[i]);
   }

   zz_pX det;
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
            negate(det, det);
         }

         MulMod(det, det, M[k][k], p);

         // make M[k, k] == -1 mod p, and make row k reduced

         InvMod(t1, M[k][k], p);
         negate(t1, t1);
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
         goto done;
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

done:
   delete[] M;
}

void inv(zz_pE& d, mat_zz_pE& X, const mat_zz_pE& A)
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
   zz_pX t1, t2;
   zz_pX *x, *y;

   const zz_pXModulus& p = zz_pE::modulus();


   vec_zz_pX *M = NTL_NEW_OP vec_zz_pX[n];

   for (i = 0; i < n; i++) {
      M[i].SetLength(2*n);
      for (j = 0; j < n; j++) {
         M[i][j].rep.SetMaxLength(2*deg(p)-1);
         M[i][j] = rep(A[i][j]);
         M[i][n+j].rep.SetMaxLength(2*deg(p)-1);
         clear(M[i][n+j]);
      }
      set(M[i][n+i]);
   }

   zz_pX det;
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
            negate(det, det);
         }

         MulMod(det, det, M[k][k], p);

         // make M[k, k] == -1 mod p, and make row k reduced

         InvMod(t1, M[k][k], p);
         negate(t1, t1);
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
         goto done;
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

done:
   delete[] M;
}



long gauss(mat_zz_pE& M_in, long w)
{
   long k, l;
   long i, j;
   long pos;
   zz_pX t1, t2, t3;
   zz_pX *x, *y;

   long n = M_in.NumRows();
   long m = M_in.NumCols();

   if (w < 0 || w > m)
      Error("gauss: bad args");

   const zz_pXModulus& p = zz_pE::modulus();


   vec_zz_pX *M = NTL_NEW_OP vec_zz_pX[n];

   for (i = 0; i < n; i++) {
      M[i].SetLength(m);
      for (j = 0; j < m; j++) {
         M[i][j].rep.SetMaxLength(2*deg(p)-1);
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
         negate(t3, t3);

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

   delete [] M;

   return l;
}

long gauss(mat_zz_pE& M)
{
   return gauss(M, M.NumCols());
}

void image(mat_zz_pE& X, const mat_zz_pE& A)
{
   mat_zz_pE M;
   M = A;
   long r = gauss(M);
   M.SetDims(r, M.NumCols());
   X = M;
}

void kernel(mat_zz_pE& X, const mat_zz_pE& A)
{
   long m = A.NumRows();
   long n = A.NumCols();

   mat_zz_pE M;
   long r;

   transpose(M, A);
   r = gauss(M);

   X.SetDims(m-r, m);

   long i, j, k, s;
   zz_pX t1, t2;

   zz_pE T3;

   vec_long D;
   D.SetLength(m);
   for (j = 0; j < m; j++) D[j] = -1;

   vec_zz_pE inverses;
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
      vec_zz_pE& v = X[k];
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
   
void mul(mat_zz_pE& X, const mat_zz_pE& A, const zz_pE& b_in)
{
   zz_pE b = b_in;
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);

   long i, j;
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         mul(X[i][j], A[i][j], b);
}

void mul(mat_zz_pE& X, const mat_zz_pE& A, const zz_p& b_in)
{
   NTL_zz_pRegister(b);
   b = b_in;
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);

   long i, j;
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         mul(X[i][j], A[i][j], b);
}

void mul(mat_zz_pE& X, const mat_zz_pE& A, long b_in)
{
   NTL_zz_pRegister(b);
   b = b_in;
   long n = A.NumRows();
   long m = A.NumCols();

   X.SetDims(n, m);

   long i, j;
   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         mul(X[i][j], A[i][j], b);
}

void diag(mat_zz_pE& X, long n, const zz_pE& d_in)  
{  
   zz_pE d = d_in;
   X.SetDims(n, n);  
   long i, j;  
  
   for (i = 1; i <= n; i++)  
      for (j = 1; j <= n; j++)  
         if (i == j)  
            X(i, j) = d;  
         else  
            clear(X(i, j));  
} 

long IsDiag(const mat_zz_pE& A, long n, const zz_pE& d)
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


long IsZero(const mat_zz_pE& a)
{
   long n = a.NumRows();
   long i;

   for (i = 0; i < n; i++)
      if (!IsZero(a[i]))
         return 0;

   return 1;
}

void clear(mat_zz_pE& x)
{
   long n = x.NumRows();
   long i;
   for (i = 0; i < n; i++)
      clear(x[i]);
}


mat_zz_pE operator+(const mat_zz_pE& a, const mat_zz_pE& b)
{
   mat_zz_pE res;
   add(res, a, b);
   NTL_OPT_RETURN(mat_zz_pE, res);
}

mat_zz_pE operator*(const mat_zz_pE& a, const mat_zz_pE& b)
{
   mat_zz_pE res;
   mul_aux(res, a, b);
   NTL_OPT_RETURN(mat_zz_pE, res);
}

mat_zz_pE operator-(const mat_zz_pE& a, const mat_zz_pE& b)
{
   mat_zz_pE res;
   sub(res, a, b);
   NTL_OPT_RETURN(mat_zz_pE, res);
}


mat_zz_pE operator-(const mat_zz_pE& a)
{
   mat_zz_pE res;
   negate(res, a);
   NTL_OPT_RETURN(mat_zz_pE, res);
}


vec_zz_pE operator*(const mat_zz_pE& a, const vec_zz_pE& b)
{
   vec_zz_pE res;
   mul_aux(res, a, b);
   NTL_OPT_RETURN(vec_zz_pE, res);
}

vec_zz_pE operator*(const vec_zz_pE& a, const mat_zz_pE& b)
{
   vec_zz_pE res;
   mul_aux(res, a, b);
   NTL_OPT_RETURN(vec_zz_pE, res);
}

void inv(mat_zz_pE& X, const mat_zz_pE& A)
{
   zz_pE d;
   inv(d, X, A);
   if (d == 0) Error("inv: non-invertible matrix");
}

void power(mat_zz_pE& X, const mat_zz_pE& A, const ZZ& e)
{
   if (A.NumRows() != A.NumCols()) Error("power: non-square matrix");

   if (e == 0) {
      ident(X, A.NumRows());
      return;
   }

   mat_zz_pE T1, T2;
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
