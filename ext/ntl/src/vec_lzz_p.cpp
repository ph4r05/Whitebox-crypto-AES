
#include <NTL/vec_lzz_p.h>

NTL_START_IMPL

void conv(vec_zz_p& x, const vec_ZZ& a)
{
   long i, n;

   n = a.length();
   x.SetLength(n);

   zz_p* xp = x.elts();
   const ZZ* ap = a.elts();

   for (i = 0; i < n; i++)
      conv(xp[i], ap[i]);
}

void conv(vec_ZZ& x, const vec_zz_p& a)
{
   long n = a.length();
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      x[i] = rep(a[i]);
}




void InnerProduct(zz_p& x, const vec_zz_p& a, const vec_zz_p& b)
{
   long n = min(a.length(), b.length());
   long i;
   zz_p accum, t;

   clear(accum);
   for (i = 0; i < n; i++) {
      mul(t, a[i], b[i]);
      add(accum, accum, t);
   }

   x = accum;
}

void InnerProduct(zz_p& x, const vec_zz_p& a, const vec_zz_p& b,
                  long offset)
{
   if (offset < 0) Error("InnerProduct: negative offset");
   if (NTL_OVERFLOW(offset, 1, 0)) Error("InnerProduct: offset too big");

   long n = min(a.length(), b.length()+offset);
   long i;
   zz_p accum, t;

   clear(accum);
   for (i = offset; i < n; i++) {
      mul(t, a[i], b[i-offset]);
      add(accum, accum, t);
   }

   x = accum;
}

long CRT(vec_ZZ& gg, ZZ& a, const vec_zz_p& G)
{
   long n = gg.length();
   if (G.length() != n) Error("CRT: vector length mismatch");

   long p = zz_p::modulus();

   ZZ new_a;
   mul(new_a, a, p);

   long a_inv;
   a_inv = rem(a, p);
   a_inv = InvMod(a_inv, p);

   long p1;
   p1 = p >> 1;

   ZZ a1;
   RightShift(a1, a, 1);

   long p_odd = (p & 1);

   long modified = 0;

   long h;

   ZZ g;
   long i;
   for (i = 0; i < n; i++) {
      if (!CRTInRange(gg[i], a)) {
         modified = 1;
         rem(g, gg[i], a);
         if (g > a1) sub(g, g, a);
      }
      else
         g = gg[i];
   
      h = rem(g, p);
      h = SubMod(rep(G[i]), h, p);
      h = MulMod(h, a_inv, p);
      if (h > p1)
         h = h - p;
   
      if (h != 0) {
         modified = 1;
   
         if (!p_odd && g > 0 && (h == p1))
            MulSubFrom(g, a, h);
         else
            MulAddTo(g, a, h);
      }

      gg[i] = g;
   }

   a = new_a;

   return modified;
}



void mul(vec_zz_p& x, const vec_zz_p& a, zz_p b)
{
   long n = a.length();
   x.SetLength(n);

   long i;

   if (n <= 1) {

      for (i = 0; i < n; i++)
	 mul(x[i], a[i], b);

   }
   else {
 
      long p = zz_p::modulus();
      double pinv = zz_p::ModulusInverse();
      long bb = rep(b);
      mulmod_precon_t bpinv = PrepMulModPrecon(bb, p, pinv);
      
      
      const zz_p *ap = a.elts();
      zz_p *xp = x.elts();

      for (i = 0; i < n; i++)
         xp[i].LoopHole() = MulModPrecon(rep(ap[i]), bb, p, bpinv);

   }
}

void mul(vec_zz_p& x, const vec_zz_p& a, long b_in)
{
   zz_p b;
   b = b_in;
   mul(x, a, b);
}



void add(vec_zz_p& x, const vec_zz_p& a, const vec_zz_p& b)
{
   long n = a.length();
   if (b.length() != n) Error("vector add: dimension mismatch");

   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      add(x[i], a[i], b[i]);
}

void sub(vec_zz_p& x, const vec_zz_p& a, const vec_zz_p& b)
{
   long n = a.length();
   if (b.length() != n) Error("vector sub: dimension mismatch");
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      sub(x[i], a[i], b[i]);
}

void clear(vec_zz_p& x)
{
   long n = x.length();
   long i;
   for (i = 0; i < n; i++)
      clear(x[i]);
}

void negate(vec_zz_p& x, const vec_zz_p& a)
{
   long n = a.length();
   x.SetLength(n);
   long i;
   for (i = 0; i < n; i++)
      negate(x[i], a[i]);
}


long IsZero(const vec_zz_p& a)
{
   long n = a.length();
   long i;

   for (i = 0; i < n; i++)
      if (!IsZero(a[i]))
         return 0;

   return 1;
}

vec_zz_p operator+(const vec_zz_p& a, const vec_zz_p& b)
{
   vec_zz_p res;
   add(res, a, b);
   NTL_OPT_RETURN(vec_zz_p, res);
}

vec_zz_p operator-(const vec_zz_p& a, const vec_zz_p& b)
{
   vec_zz_p res;
   sub(res, a, b);
   NTL_OPT_RETURN(vec_zz_p, res);
}


vec_zz_p operator-(const vec_zz_p& a)
{
   vec_zz_p res;
   negate(res, a);
   NTL_OPT_RETURN(vec_zz_p, res);
}


zz_p operator*(const vec_zz_p& a, const vec_zz_p& b)
{
   zz_p res;
   InnerProduct(res, a, b);
   return res;
}


void VectorCopy(vec_zz_p& x, const vec_zz_p& a, long n)
{
   if (n < 0) Error("VectorCopy: negative length");
   if (NTL_OVERFLOW(n, 1, 0)) Error("overflow in VectorCopy");

   long m = min(n, a.length());

   x.SetLength(n);
  
   long i;

   for (i = 0; i < m; i++)
      x[i] = a[i];

   for (i = m; i < n; i++)
      clear(x[i]);
}

NTL_END_IMPL
