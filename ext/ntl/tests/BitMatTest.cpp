
#include <NTL/mat_GF2.h>
#include <NTL/mat_lzz_p.h>

NTL_CLIENT


void random(mat_zz_p& X, long n, long m)
{
   X.SetDims(n, m);
   long i, j;

   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         random(X[i][j]);
}

void random(vec_zz_p& X, long n)
{
   X.SetLength(n);
   long i;

   for (i = 0; i < n; i++)
      random(X[i]);
}

void cvt(mat_GF2& x, const mat_zz_p& a)
{
   long n = a.NumRows();
   long m = a.NumCols();

   x.SetDims(n, m);

   long i, j;

   for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
         x.put(i, j, rep(a[i][j]));
}


void cvt(vec_GF2& x, const vec_zz_p& a)
{
   long n = a.length();

   x.SetLength(n);

   long i;

   for (i = 0; i < n; i++)
      x.put(i, rep(a[i]));
}

int main()
{
   zz_p::init(2);

   long i;

   vec_GF2 v;
   v.SetLength(5);
   v[1] = 1;
   v[0] = v[1];

   if (v[0] != v[1]) Error("BitMatTest not OK!!");

   for (i=0; i < 8; i++) {
      mat_zz_p a, x;
      mat_GF2 A, X, X1;

      long n = RandomBnd(500) + 1;
      long m = RandomBnd(500) + 1;
      cerr << n << " " << m << "\n";

      double t;

      random(a, n, m);

      t = GetTime();
      kernel(x, a);
      t = GetTime() - t;  cerr << t << "\n";

      cvt(A, a);

      t = GetTime();
      kernel(X, A);
      t = GetTime() - t;  cerr << t << "\n";

      cerr << x.NumRows() << "\n";

      cvt(X1, x);

      if (X1 != X) Error("BitMatTest NOT OK!!");

      if (!IsZero(X*A)) Error("BitMatTest NOT OK!!");

      cerr << "\n";
   }

   cerr << "BitMatTest OK\n";

}

