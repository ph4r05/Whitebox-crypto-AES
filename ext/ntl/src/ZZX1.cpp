

#include <NTL/ZZX.h>

#include <NTL/new.h>

NTL_START_IMPL





void conv(zz_pX& x, const ZZX& a)
{
   conv(x.rep, a.rep);
   x.normalize();
}


void conv(ZZX& x, const zz_pX& a)
{
   conv(x.rep, a.rep);
   x.normalize();
}


long CRT(ZZX& gg, ZZ& a, const zz_pX& G)
{
   long n = gg.rep.length();

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

   long m = G.rep.length();

   long max_mn = max(m, n);

   gg.rep.SetLength(max_mn);

   ZZ g;
   long i;

   for (i = 0; i < n; i++) {
      if (!CRTInRange(gg.rep[i], a)) {
         modified = 1;
         rem(g, gg.rep[i], a);
         if (g > a1) sub(g, g, a);
      }
      else
         g = gg.rep[i];
   
      h = rem(g, p);

      if (i < m)
         h = SubMod(rep(G.rep[i]), h, p);
      else
         h = NegateMod(h, p);

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

      gg.rep[i] = g;
   }


   for (; i < m; i++) {
      h = rep(G.rep[i]);
      h = MulMod(h, a_inv, p);
      if (h > p1)
         h = h - p;
   
      modified = 1;
      mul(g, a, h);
      gg.rep[i] = g;
   }

   gg.normalize();
   a = new_a;

   return modified;
}

long CRT(ZZX& gg, ZZ& a, const ZZ_pX& G)
{
   long n = gg.rep.length();

   const ZZ& p = ZZ_p::modulus();

   ZZ new_a;
   mul(new_a, a, p);

   ZZ a_inv;
   rem(a_inv, a, p);
   InvMod(a_inv, a_inv, p);

   ZZ p1;
   RightShift(p1, p, 1);

   ZZ a1;
   RightShift(a1, a, 1);

   long p_odd = IsOdd(p);

   long modified = 0;

   ZZ h;
   ZZ ah;

   long m = G.rep.length();

   long max_mn = max(m, n);

   gg.rep.SetLength(max_mn);

   ZZ g;
   long i;

   for (i = 0; i < n; i++) {
      if (!CRTInRange(gg.rep[i], a)) {
         modified = 1;
         rem(g, gg.rep[i], a);
         if (g > a1) sub(g, g, a);
      }
      else
         g = gg.rep[i];
   
      rem(h, g, p);

      if (i < m)
         SubMod(h, rep(G.rep[i]), h, p);
      else
         NegateMod(h, h, p);

      MulMod(h, h, a_inv, p);
      if (h > p1)
         sub(h, h, p);
   
      if (h != 0) {
         modified = 1;
         mul(ah, a, h);
   
         if (!p_odd && g > 0 && (h == p1))
            sub(g, g, ah);
         else
            add(g, g, ah);
      }

      gg.rep[i] = g;
   }


   for (; i < m; i++) {
      h = rep(G.rep[i]);
      MulMod(h, h, a_inv, p);
      if (h > p1)
         sub(h, h, p);
   
      modified = 1;
      mul(g, a, h);
      gg.rep[i] = g;
   }

   gg.normalize();
   a = new_a;

   return modified;
}




/* Compute a = b * 2^l mod p, where p = 2^n+1. 0<=l<=n and 0<b<p are
   assumed. */
static void LeftRotate(ZZ& a, const ZZ& b, long l, const ZZ& p, long n)
{
  if (l == 0) {
    if (&a != &b) {
      a = b;
    }
    return;
  }

  /* tmp := upper l bits of b */
  static ZZ tmp;
  RightShift(tmp, b, n - l);
  /* a := 2^l * lower n - l bits of b */
  trunc(a, b, n - l);
  LeftShift(a, a, l);
  /* a -= tmp */
  sub(a, a, tmp);
  if (sign(a) < 0) {
    add(a, a, p);
  }
}


/* Compute a = b * 2^l mod p, where p = 2^n+1. 0<=p<b is assumed. */
static void Rotate(ZZ& a, const ZZ& b, long l, const ZZ& p, long n)
{
  if (IsZero(b)) {
    clear(a);
    return;
  }

  /* l %= 2n */
  if (l >= 0) {
    l %= (n << 1);
  } else {
    l = (n << 1) - 1 - (-(l + 1) % (n << 1));
  }

  /* a = b * 2^l mod p */
  if (l < n) {
    LeftRotate(a, b, l, p, n);
  } else {
    LeftRotate(a, b, l - n, p, n);
    SubPos(a, p, a);
  }
}



/* Fast Fourier Transform. a is a vector of length 2^l, 2^l divides 2n,
   p = 2^n+1, w = 2^r mod p is a primitive (2^l)th root of
   unity. Returns a(1),a(w),...,a(w^{2^l-1}) mod p in bit-reverse
   order. */
static void fft(vec_ZZ& a, long r, long l, const ZZ& p, long n)
{
  long round;
  long off, i, j, e;
  long halfsize;
  ZZ tmp, tmp1;

  for (round = 0; round < l; round++, r <<= 1) {
    halfsize =  1L << (l - 1 - round);
    for (i = (1L << round) - 1, off = 0; i >= 0; i--, off += halfsize) {
      for (j = 0, e = 0; j < halfsize; j++, off++, e+=r) {
	/* One butterfly : 
	 ( a[off], a[off+halfsize] ) *= ( 1  w^{j2^round} )
	                                ( 1 -w^{j2^round} ) */
	/* tmp = a[off] - a[off + halfsize] mod p */
	sub(tmp, a[off], a[off + halfsize]);
	if (sign(tmp) < 0) {
	  add(tmp, tmp, p);
	}
	/* a[off] += a[off + halfsize] mod p */
	add(a[off], a[off], a[off + halfsize]);
	sub(tmp1, a[off], p);
	if (sign(tmp1) >= 0) {
	  a[off] = tmp1;
	}
	/* a[off + halfsize] = tmp * w^{j2^round} mod p */
	Rotate(a[off + halfsize], tmp, e, p, n);
      }
    }
  }
}

/* Inverse FFT. r must be the same as in the call to FFT. Result is
   by 2^l too large. */
static void ifft(vec_ZZ& a, long r, long l, const ZZ& p, long n)
{
  long round;
  long off, i, j, e;
  long halfsize;
  ZZ tmp, tmp1;

  for (round = l - 1, r <<= l - 1; round >= 0; round--, r >>= 1) {
    halfsize = 1L << (l - 1 - round);
    for (i = (1L << round) - 1, off = 0; i >= 0; i--, off += halfsize) {
      for (j = 0, e = 0; j < halfsize; j++, off++, e+=r) {
	/* One inverse butterfly : 
	 ( a[off], a[off+halfsize] ) *= ( 1               1             )
	                                ( w^{-j2^round}  -w^{-j2^round} ) */
	/* a[off + halfsize] *= w^{-j2^round} mod p */
	Rotate(a[off + halfsize], a[off + halfsize], -e, p, n);
	/* tmp = a[off] - a[off + halfsize] */
	sub(tmp, a[off], a[off + halfsize]);

	/* a[off] += a[off + halfsize] mod p */
	add(a[off], a[off], a[off + halfsize]);
	sub(tmp1, a[off], p);
	if (sign(tmp1) >= 0) {
	  a[off] = tmp1;
	}
	/* a[off+halfsize] = tmp mod p */
	if (sign(tmp) < 0) {
	  add(a[off+halfsize], tmp, p);
	} else {
	  a[off+halfsize] = tmp;
	}
      }
    }
  }
}



/* Multiplication a la Schoenhage & Strassen, modulo a "Fermat" number
   p = 2^{mr}+1, where m is a power of two and r is odd. Then w = 2^r
   is a primitive 2mth root of unity, i.e., polynomials whose product
   has degree less than 2m can be multiplied, provided that the
   coefficients of the product polynomial are at most 2^{mr-1} in
   absolute value. The algorithm is not called recursively;
   coefficient arithmetic is done directly.*/

void SSMul(ZZX& c, const ZZX& a, const ZZX& b)
{
  long na = deg(a);
  long nb = deg(b);

  if (na <= 0 || nb <= 0) {
    PlainMul(c, a, b);
    return;
  }

  long n = na + nb; /* degree of the product */


  /* Choose m and r suitably */
  long l = NextPowerOfTwo(n + 1) - 1; /* 2^l <= n < 2^{l+1} */
  long m2 = 1L << (l + 1); /* m2 = 2m = 2^{l+1} */
  /* Bitlength of the product: if the coefficients of a are absolutely less
     than 2^ka and the coefficients of b are absolutely less than 2^kb, then
     the coefficients of ab are absolutely less than
     (min(na,nb)+1)2^{ka+kb} <= 2^bound. */
  long bound = 2 + NumBits(min(na, nb)) + MaxBits(a) + MaxBits(b);
  /* Let r be minimal so that mr > bound */
  long r = (bound >> l) + 1;
  long mr = r << l;

  /* p := 2^{mr}+1 */
  ZZ p;
  set(p);
  LeftShift(p, p, mr);
  add(p, p, 1);

  /* Make coefficients of a and b positive */
  vec_ZZ aa, bb;
  aa.SetLength(m2);
  bb.SetLength(m2);

  long i;
  for (i = 0; i <= deg(a); i++) {
    if (sign(a.rep[i]) >= 0) {
      aa[i] = a.rep[i];
    } else {
      add(aa[i], a.rep[i], p);
    }
  }

  for (i = 0; i <= deg(b); i++) {
    if (sign(b.rep[i]) >= 0) {
      bb[i] = b.rep[i];
    } else {
      add(bb[i], b.rep[i], p);
    }
  }

  /* 2m-point FFT's mod p */
  fft(aa, r, l + 1, p, mr);
  fft(bb, r, l + 1, p, mr);

  /* Pointwise multiplication aa := aa * bb mod p */
  ZZ tmp, ai;
  for (i = 0; i < m2; i++) {
    mul(ai, aa[i], bb[i]);
    if (NumBits(ai) > mr) {
      RightShift(tmp, ai, mr);
      trunc(ai, ai, mr);
      sub(ai, ai, tmp);
      if (sign(ai) < 0) {
	add(ai, ai, p);
      }
    }
    aa[i] = ai;
  }
  
  ifft(aa, r, l + 1, p, mr);

  /* Retrieve c, dividing by 2m, and subtracting p where necessary */
  c.rep.SetLength(n + 1);
  for (i = 0; i <= n; i++) {
    ai = aa[i];
    ZZ& ci = c.rep[i];
    if (!IsZero(ai)) {
      /* ci = -ai * 2^{mr-l-1} = ai * 2^{-l-1} = ai / 2m mod p */
      LeftRotate(ai, ai, mr - l - 1, p, mr);
      sub(tmp, p, ai);
      if (NumBits(tmp) >= mr) { /* ci >= (p-1)/2 */
	negate(ci, ai); /* ci = -ai = ci - p */
      }
      else
        ci = tmp;
    } 
    else
       clear(ci);
  }
}


// SSRatio computes how much bigger the SS modulus must be
// to accomodate the necessary roots of unity.
// This is useful in determining algorithm crossover points.

double SSRatio(long na, long maxa, long nb, long maxb)
{
  if (na <= 0 || nb <= 0) return 0;

  long n = na + nb; /* degree of the product */


  long l = NextPowerOfTwo(n + 1) - 1; /* 2^l <= n < 2^{l+1} */
  long bound = 2 + NumBits(min(na, nb)) + maxa + maxb;
  long r = (bound >> l) + 1;
  long mr = r << l;

  return double(mr + 1)/double(bound);
}

void HomMul(ZZX& x, const ZZX& a, const ZZX& b)
{
   if (&a == &b) {
      HomSqr(x, a);
      return;
   }

   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   long bound = 2 + NumBits(min(da, db)+1) + MaxBits(a) + MaxBits(b);


   ZZ prod;
   set(prod);

   long i, nprimes;

   zz_pBak bak;
   bak.save();

   for (nprimes = 0; NumBits(prod) <= bound; nprimes++) {
      if (nprimes >= NumFFTPrimes)
         zz_p::FFTInit(nprimes);
      mul(prod, prod, FFTPrime[nprimes]);
   }


   ZZ coeff;
   ZZ t1;
   long tt;

   vec_ZZ c;

   c.SetLength(da+db+1);

   long j;

   for (i = 0; i < nprimes; i++) {
      zz_p::FFTInit(i);
      long p = zz_p::modulus();

      div(t1, prod, p);
      tt = rem(t1, p);
      tt = InvMod(tt, p);
      mul(coeff, t1, tt);

      zz_pX A, B, C;

      conv(A, a);
      conv(B, b);
      mul(C, A, B);

      long m = deg(C);

      for (j = 0; j <= m; j++) {
         /* c[j] += coeff*rep(C.rep[j]) */
         MulAddTo(c[j], coeff, rep(C.rep[j]));
         // mul(t1, coeff, rep(C.rep[j]));
         // add(c[j], c[j], t1); 
      }
   }

   x.rep.SetLength(da+db+1);

   ZZ prod2;
   RightShift(prod2, prod, 1);

   for (j = 0; j <= da+db; j++) {
      rem(t1, c[j], prod);

      if (t1 > prod2)
         sub(x.rep[j], t1, prod);
      else
         x.rep[j] = t1;
   }

   x.normalize();

   bak.restore();
}

static
long MaxSize(const ZZX& a)
{
   long res = 0;
   long n = a.rep.length();

   long i;
   for (i = 0; i < n; i++) {
      long t = a.rep[i].size();
      if (t > res)
         res = t;
   }

   return res;
}


void mul(ZZX& c, const ZZX& a, const ZZX& b)
{
   if (IsZero(a) || IsZero(b)) {
      clear(c);
      return;
   }

   if (&a == &b) {
      sqr(c, a);
      return;
   }

   long maxa = MaxSize(a);
   long maxb = MaxSize(b);

   long k = min(maxa, maxb);
   long s = min(deg(a), deg(b)) + 1;

   if (s == 1 || (k == 1 && s < 40) || (k == 2 && s < 20) || 
                 (k == 3 && s < 10)) {

      PlainMul(c, a, b);
      return;
   }

   if (s < 80 || (k < 30 && s < 150))  {
      KarMul(c, a, b);
      return;
   }


   if (maxa + maxb >= 40 && 
       SSRatio(deg(a), MaxBits(a), deg(b), MaxBits(b)) < 1.75) 
      SSMul(c, a, b);
   else
      HomMul(c, a, b);
}


void SSSqr(ZZX& c, const ZZX& a)
{
  long na = deg(a);
  if (na <= 0) {
    PlainSqr(c, a);
    return;
  }

  long n = na + na; /* degree of the product */


  long l = NextPowerOfTwo(n + 1) - 1; /* 2^l <= n < 2^{l+1} */
  long m2 = 1L << (l + 1); /* m2 = 2m = 2^{l+1} */
  long bound = 2 + NumBits(na) + 2*MaxBits(a);
  long r = (bound >> l) + 1;
  long mr = r << l;

  /* p := 2^{mr}+1 */
  ZZ p;
  set(p);
  LeftShift(p, p, mr);
  add(p, p, 1);

  vec_ZZ aa;
  aa.SetLength(m2);

  long i;
  for (i = 0; i <= deg(a); i++) {
    if (sign(a.rep[i]) >= 0) {
      aa[i] = a.rep[i];
    } else {
      add(aa[i], a.rep[i], p);
    }
  }


  /* 2m-point FFT's mod p */
  fft(aa, r, l + 1, p, mr);

  /* Pointwise multiplication aa := aa * aa mod p */
  ZZ tmp, ai;
  for (i = 0; i < m2; i++) {
    sqr(ai, aa[i]);
    if (NumBits(ai) > mr) {
      RightShift(tmp, ai, mr);
      trunc(ai, ai, mr);
      sub(ai, ai, tmp);
      if (sign(ai) < 0) {
	add(ai, ai, p);
      }
    }
    aa[i] = ai;
  }
  
  ifft(aa, r, l + 1, p, mr);

  ZZ ci;

  /* Retrieve c, dividing by 2m, and subtracting p where necessary */
  c.rep.SetLength(n + 1);

  for (i = 0; i <= n; i++) {
    ai = aa[i];
    ZZ& ci = c.rep[i];
    if (!IsZero(ai)) {
      /* ci = -ai * 2^{mr-l-1} = ai * 2^{-l-1} = ai / 2m mod p */
      LeftRotate(ai, ai, mr - l - 1, p, mr);
      sub(tmp, p, ai);
      if (NumBits(tmp) >= mr) { /* ci >= (p-1)/2 */
	negate(ci, ai); /* ci = -ai = ci - p */
      }
      else
        ci = tmp;
    } 
    else
       clear(ci);
  }
}

void HomSqr(ZZX& x, const ZZX& a)
{

   long da = deg(a);

   if (da < 0) {
      clear(x);
      return;
   }

   long bound = 2 + NumBits(da+1) + 2*MaxBits(a);


   ZZ prod;
   set(prod);

   long i, nprimes;

   zz_pBak bak;
   bak.save();

   for (nprimes = 0; NumBits(prod) <= bound; nprimes++) {
      if (nprimes >= NumFFTPrimes)
         zz_p::FFTInit(nprimes);
      mul(prod, prod, FFTPrime[nprimes]);
   }


   ZZ coeff;
   ZZ t1;
   long tt;

   vec_ZZ c;

   c.SetLength(da+da+1);

   long j;

   for (i = 0; i < nprimes; i++) {
      zz_p::FFTInit(i);
      long p = zz_p::modulus();

      div(t1, prod, p);
      tt = rem(t1, p);
      tt = InvMod(tt, p);
      mul(coeff, t1, tt);

      zz_pX A, C;

      conv(A, a);
      sqr(C, A);

      long m = deg(C);

      for (j = 0; j <= m; j++) {
         /* c[j] += coeff*rep(C.rep[j]) */
         MulAddTo(c[j], coeff, rep(C.rep[j]));
         // mul(t1, coeff, rep(C.rep[j]));
         // add(c[j], c[j], t1); 
      }
   }

   x.rep.SetLength(da+da+1);

   ZZ prod2;
   RightShift(prod2, prod, 1);

   for (j = 0; j <= da+da; j++) {
      rem(t1, c[j], prod);

      if (t1 > prod2)
         sub(x.rep[j], t1, prod);
      else
         x.rep[j] = t1;
   }

   x.normalize();

   bak.restore();
}


void sqr(ZZX& c, const ZZX& a)
{
   if (IsZero(a)) {
      clear(c);
      return;
   }

   long maxa = MaxSize(a);

   long k = maxa;
   long s = deg(a) + 1;

   if (s == 1 || (k == 1 && s < 50) || (k == 2 && s < 25) || 
                 (k == 3 && s < 25) || (k == 4 && s < 10)) {

      PlainSqr(c, a);
      return;
   }

   if (s < 80 || (k < 30 && s < 150))  {
      KarSqr(c, a);
      return;
   }

   long mba = MaxBits(a);
   
   if (2*maxa >= 40 && 
       SSRatio(deg(a), mba, deg(a), mba) < 1.75) 
      SSSqr(c, a);
   else
      HomSqr(c, a);
}


void mul(ZZX& x, const ZZX& a, const ZZ& b)
{
   ZZ t;
   long i, da;

   const ZZ *ap;
   ZZ* xp;

   if (IsZero(b)) {
      clear(x);
      return;
   }

   t = b;
   da = deg(a);
   x.rep.SetLength(da+1);
   ap = a.rep.elts();
   xp = x.rep.elts();

   for (i = 0; i <= da; i++) 
      mul(xp[i], ap[i], t);
}

void mul(ZZX& x, const ZZX& a, long b)
{
   long i, da;

   const ZZ *ap;
   ZZ* xp;

   if (b == 0) {
      clear(x);
      return;
   }

   da = deg(a);
   x.rep.SetLength(da+1);
   ap = a.rep.elts();
   xp = x.rep.elts();

   for (i = 0; i <= da; i++) 
      mul(xp[i], ap[i], b);
}




void diff(ZZX& x, const ZZX& a)
{
   long n = deg(a);
   long i;

   if (n <= 0) {
      clear(x);
      return;
   }

   if (&x != &a)
      x.rep.SetLength(n);

   for (i = 0; i <= n-1; i++) {
      mul(x.rep[i], a.rep[i+1], i+1);
   }

   if (&x == &a)
      x.rep.SetLength(n);

   x.normalize();
}

void HomPseudoDivRem(ZZX& q, ZZX& r, const ZZX& a, const ZZX& b)
{
   if (IsZero(b)) Error("division by zero");

   long da = deg(a);
   long db = deg(b);

   if (da < db) {
      r = a;
      clear(q);
      return;
   }

   ZZ LC;
   LC = LeadCoeff(b);

   ZZ LC1;

   power(LC1, LC, da-db+1);

   long a_bound = NumBits(LC1) + MaxBits(a);

   LC1.kill();

   long b_bound = MaxBits(b);

   zz_pBak bak;
   bak.save();

   ZZX qq, rr;

   ZZ prod, t;
   set(prod);

   clear(qq);
   clear(rr);

   long i;
   long Qinstable, Rinstable;

   Qinstable = 1;
   Rinstable = 1;

   for (i = 0; ; i++) {
      zz_p::FFTInit(i);
      long p = zz_p::modulus();


      if (divide(LC, p)) continue;

      zz_pX A, B, Q, R;

      conv(A, a);
      conv(B, b);
      
      if (!IsOne(LC)) {
         zz_p y;
         conv(y, LC);
         power(y, y, da-db+1);
         mul(A, A, y);
      }

      if (!Qinstable) {
         conv(Q, qq);
         mul(R, B, Q);
         sub(R, A, R);

         if (deg(R) >= db)
            Qinstable = 1;
         else
            Rinstable = CRT(rr, prod, R);
      }

      if (Qinstable) {
         DivRem(Q, R, A, B);
         t = prod;
         Qinstable = CRT(qq, t, Q);
         Rinstable =  CRT(rr, prod, R);
      }

      if (!Qinstable && !Rinstable) {
         // stabilized...check if prod is big enough

         long bound1 = b_bound + MaxBits(qq) + NumBits(min(db, da-db)+1);
         long bound2 = MaxBits(rr);
         long bound = max(bound1, bound2);

         if (a_bound > bound)
            bound = a_bound;

         bound += 4;

         if (NumBits(prod) > bound)
            break;
      }
   }

   bak.restore();

   q = qq;
   r = rr;
}




void HomPseudoDiv(ZZX& q, const ZZX& a, const ZZX& b)
{
   ZZX r;
   HomPseudoDivRem(q, r, a, b);
}

void HomPseudoRem(ZZX& r, const ZZX& a, const ZZX& b)
{
   ZZX q;
   HomPseudoDivRem(q, r, a, b);
}

void PlainPseudoDivRem(ZZX& q, ZZX& r, const ZZX& a, const ZZX& b)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ *bp;
   ZZ *qp;
   ZZ *xp;


   ZZ  s, t;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZX: division by zero");

   if (da < db) {
      r = a;
      clear(q);
      return;
   }

   ZZX lb;

   if (&q == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   ZZ LC = bp[db];
   LCIsOne = IsOne(LC);


   vec_ZZ x;

   x = a.rep;
   xp = x.elts();

   dq = da - db;
   q.rep.SetLength(dq+1);
   qp = q.rep.elts();

   if (!LCIsOne) {
      t = LC;
      for (i = dq-1; i >= 0; i--) {
         mul(xp[i], xp[i], t);
         if (i > 0) mul(t, t, LC);
      }
   }

   for (i = dq; i >= 0; i--) {
      t = xp[i+db];
      qp[i] = t;

      for (j = db-1; j >= 0; j--) {
	 mul(s, t, bp[j]);
         if (!LCIsOne) mul(xp[i+j], xp[i+j], LC);
	 sub(xp[i+j], xp[i+j], s);
      }
   }

   if (!LCIsOne) {
      t = LC;
      for (i = 1; i <= dq; i++) {
         mul(qp[i], qp[i], t);
         if (i < dq) mul(t, t, LC);
      }
   }
      

   r.rep.SetLength(db);
   for (i = 0; i < db; i++)
      r.rep[i] = xp[i];
   r.normalize();
}


void PlainPseudoDiv(ZZX& q, const ZZX& a, const ZZX& b)
{
   ZZX r;
   PlainPseudoDivRem(q, r, a, b);
}

void PlainPseudoRem(ZZX& r, const ZZX& a, const ZZX& b)
{
   ZZX q;
   PlainPseudoDivRem(q, r, a, b);
}

void div(ZZX& q, const ZZX& a, long b)
{
   if (b == 0) Error("div: division by zero");

   if (!divide(q, a, b)) Error("DivRem: quotient undefined over ZZ");
}

void div(ZZX& q, const ZZX& a, const ZZ& b)
{
   if (b == 0) Error("div: division by zero");

   if (!divide(q, a, b)) Error("DivRem: quotient undefined over ZZ");
}

static
void ConstDivRem(ZZX& q, ZZX& r, const ZZX& a, const ZZ& b)
{
   if (b == 0) Error("DivRem: division by zero");

   if (!divide(q, a, b)) Error("DivRem: quotient undefined over ZZ");

   r = 0;
}

static
void ConstRem(ZZX& r, const ZZX& a, const ZZ& b)
{
   if (b == 0) Error("rem: division by zero");

   r = 0;
}



void DivRem(ZZX& q, ZZX& r, const ZZX& a, const ZZX& b)
{
   long da = deg(a);
   long db = deg(b);

   if (db < 0) Error("DivRem: division by zero");

   if (da < db) {
      r = a;
      q = 0;
   }
   else if (db == 0) {
      ConstDivRem(q, r, a, ConstTerm(b));
   }
   else if (IsOne(LeadCoeff(b))) {
      PseudoDivRem(q, r, a, b);
   }
   else if (LeadCoeff(b) == -1) {
      ZZX b1;
      negate(b1, b);
      PseudoDivRem(q, r, a, b1);
      negate(q, q);
   }
   else if (divide(q, a, b)) {
      r = 0;
   }
   else {
      ZZX q1, r1;
      ZZ m;
      PseudoDivRem(q1, r1, a, b);
      power(m, LeadCoeff(b), da-db+1);
      if (!divide(q, q1, m)) Error("DivRem: quotient not defined over ZZ");
      if (!divide(r, r1, m)) Error("DivRem: remainder not defined over ZZ");
   }
}

void div(ZZX& q, const ZZX& a, const ZZX& b)
{
   long da = deg(a);
   long db = deg(b);

   if (db < 0) Error("div: division by zero");

   if (da < db) {
      q = 0;
   }
   else if (db == 0) {
      div(q, a, ConstTerm(b));
   }
   else if (IsOne(LeadCoeff(b))) {
      PseudoDiv(q, a, b);
   }
   else if (LeadCoeff(b) == -1) {
      ZZX b1;
      negate(b1, b);
      PseudoDiv(q, a, b1);
      negate(q, q);
   }
   else if (divide(q, a, b)) {

      // nothing to do
      
   }
   else {
      ZZX q1;
      ZZ m;
      PseudoDiv(q1, a, b);
      power(m, LeadCoeff(b), da-db+1);
      if (!divide(q, q1, m)) Error("div: quotient not defined over ZZ");
   }
}

void rem(ZZX& r, const ZZX& a, const ZZX& b)
{
   long da = deg(a);
   long db = deg(b);

   if (db < 0) Error("rem: division by zero");

   if (da < db) {
      r = a;
   }
   else if (db == 0) {
      ConstRem(r, a, ConstTerm(b));
   }
   else if (IsOne(LeadCoeff(b))) {
      PseudoRem(r, a, b);
   }
   else if (LeadCoeff(b) == -1) {
      ZZX b1;
      negate(b1, b);
      PseudoRem(r, a, b1);
   }
   else if (divide(a, b)) {
      r = 0;
   }
   else {
      ZZX r1;
      ZZ m;
      PseudoRem(r1, a, b);
      power(m, LeadCoeff(b), da-db+1);
      if (!divide(r, r1, m)) Error("rem: remainder not defined over ZZ");
   }
}

long HomDivide(ZZX& q, const ZZX& a, const ZZX& b)
{
   if (IsZero(b)) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }

   if (IsZero(a)) {
      clear(q);
      return 1;
   }

   if (deg(b) == 0) {
      return divide(q, a, ConstTerm(b));
   }

   if (deg(a) < deg(b)) return 0;

   ZZ ca, cb, cq;

   content(ca, a);
   content(cb, b);

   if (!divide(cq, ca, cb)) return 0;

   ZZX aa, bb;

   divide(aa, a, ca);
   divide(bb, b, cb);

   if (!divide(LeadCoeff(aa), LeadCoeff(bb)))
      return 0;

   if (!divide(ConstTerm(aa), ConstTerm(bb)))
      return 0;

   zz_pBak bak;
   bak.save();

   ZZX qq;

   ZZ prod;
   set(prod);

   clear(qq);
   long res = 1;
   long Qinstable = 1;


   long a_bound = MaxBits(aa);
   long b_bound = MaxBits(bb);


   long i;
   for (i = 0; ; i++) {
      zz_p::FFTInit(i);
      long p = zz_p::modulus();

      if (divide(LeadCoeff(bb), p)) continue;

      zz_pX A, B, Q, R;

      conv(A, aa);
      conv(B, bb);

      if (!Qinstable) {
         conv(Q, qq);
         mul(R, B, Q);
         sub(R, A, R);

         if (deg(R) >= deg(B))
            Qinstable = 1;
         else if (!IsZero(R)) {
            res = 0;
            break;
         }
         else
            mul(prod, prod, p);
      }

      if (Qinstable) {
         if (!divide(Q, A, B)) {
            res = 0;
            break;
         }

         Qinstable = CRT(qq, prod, Q);
      }

      if (!Qinstable) {
         // stabilized...check if prod is big enough

         long bound = b_bound + MaxBits(qq) + 
                     NumBits(min(deg(bb), deg(qq)) + 1);

         if (a_bound > bound)
            bound = a_bound;

         bound += 3;

         if (NumBits(prod) > bound) 
            break;
      }
   }

   bak.restore();

   if (res) mul(q, qq, cq);
   return res;

}


long HomDivide(const ZZX& a, const ZZX& b)
{
   if (deg(b) == 0) {
      return divide(a, ConstTerm(b));
   }
   else {
      ZZX q;
      return HomDivide(q, a, b);
   }
}

long PlainDivide(ZZX& qq, const ZZX& aa, const ZZX& bb)
{
   if (IsZero(bb)) {
      if (IsZero(aa)) {
         clear(qq);
         return 1;
      }
      else
         return 0;
   }

   if (deg(bb) == 0) {
      return divide(qq, aa, ConstTerm(bb));
   }

   long da, db, dq, i, j, LCIsOne;
   const ZZ *bp;
   ZZ *qp;
   ZZ *xp;


   ZZ  s, t;

   da = deg(aa);
   db = deg(bb);

   if (da < db) {
      return 0;
   }

   ZZ ca, cb, cq;

   content(ca, aa);
   content(cb, bb);

   if (!divide(cq, ca, cb)) {
      return 0;
   } 


   ZZX a, b, q;

   divide(a, aa, ca);
   divide(b, bb, cb);

   if (!divide(LeadCoeff(a), LeadCoeff(b)))
      return 0;

   if (!divide(ConstTerm(a), ConstTerm(b)))
      return 0;

   long coeff_bnd = MaxBits(a) + (NumBits(da+1)+1)/2 + (da-db);

   bp = b.rep.elts();

   ZZ LC;
   LC = bp[db];

   LCIsOne = IsOne(LC);

   xp = a.rep.elts();

   dq = da - db;
   q.rep.SetLength(dq+1);
   qp = q.rep.elts();

   for (i = dq; i >= 0; i--) {
      if (!LCIsOne) {
         if (!divide(t, xp[i+db], LC))
            return 0;
      }
      else
         t = xp[i+db];

      if (NumBits(t) > coeff_bnd) return 0;

      qp[i] = t;

      for (j = db-1; j >= 0; j--) {
	 mul(s, t, bp[j]);
	 sub(xp[i+j], xp[i+j], s);
      }
   }

   for (i = 0; i < db; i++)
      if (!IsZero(xp[i]))
         return 0;

   mul(qq, q, cq);
   return 1;
}

long PlainDivide(const ZZX& a, const ZZX& b)
{
   if (deg(b) == 0) 
      return divide(a, ConstTerm(b));
   else {
      ZZX q;
      return PlainDivide(q, a, b);
   }
}


long divide(ZZX& q, const ZZX& a, const ZZX& b)
{
   long da = deg(a);
   long db = deg(b);

   if (db <= 8 || da-db <= 8)
      return PlainDivide(q, a, b);
   else
      return HomDivide(q, a, b);
}

long divide(const ZZX& a, const ZZX& b)
{
   long da = deg(a);
   long db = deg(b);

   if (db <= 8 || da-db <= 8)
      return PlainDivide(a, b);
   else
      return HomDivide(a, b);
}







long divide(ZZX& q, const ZZX& a, const ZZ& b)
{
   if (IsZero(b)) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }

   if (IsOne(b)) {
      q = a;
      return 1;
   }

   if (b == -1) {
      negate(q, a);
      return 1;
   }

   long n = a.rep.length();
   vec_ZZ res(INIT_SIZE, n);
   long i;

   for (i = 0; i < n; i++) {
      if (!divide(res[i], a.rep[i], b))
         return 0;
   }

   q.rep = res;
   return 1;
}

long divide(const ZZX& a, const ZZ& b)
{
   if (IsZero(b)) return IsZero(a);

   if (IsOne(b) || b == -1) {
      return 1;
   }

   long n = a.rep.length();
   long i;

   for (i = 0; i < n; i++) {
      if (!divide(a.rep[i], b))
         return 0;
   }

   return 1;
}

long divide(ZZX& q, const ZZX& a, long b)
{
   if (b == 0) {
      if (IsZero(a)) {
         clear(q);
         return 1;
      }
      else
         return 0;
   }

   if (b == 1) {
      q = a;
      return 1;
   }

   if (b == -1) {
      negate(q, a);
      return 1;
   }

   long n = a.rep.length();
   vec_ZZ res(INIT_SIZE, n);
   long i;

   for (i = 0; i < n; i++) {
      if (!divide(res[i], a.rep[i], b))
         return 0;
   }

   q.rep = res;
   return 1;
}

long divide(const ZZX& a, long b)
{
   if (b == 0) return IsZero(a);
   if (b == 1 || b == -1) {
      return 1;
   }

   long n = a.rep.length();
   long i;

   for (i = 0; i < n; i++) {
      if (!divide(a.rep[i], b))
         return 0;
   }

   return 1;
}

   

void content(ZZ& d, const ZZX& f)
{
   ZZ res;
   long i;

   clear(res);
   for (i = 0; i <= deg(f); i++) {
      GCD(res, res, f.rep[i]);
      if (IsOne(res)) break;
   }

   if (sign(LeadCoeff(f)) < 0) negate(res, res);
   d = res;
}

void PrimitivePart(ZZX& pp, const ZZX& f)
{
   if (IsZero(f)) {
      clear(pp);
      return;
   }
 
   ZZ d;

   content(d, f);
   divide(pp, f, d);
}


static
void BalCopy(ZZX& g, const zz_pX& G)
{
   long p = zz_p::modulus();
   long p2 = p >> 1;
   long n = G.rep.length();
   long i;
   long t;

   g.rep.SetLength(n);
   for (i = 0; i < n; i++) {
      t = rep(G.rep[i]);
      if (t > p2) t = t - p;
      conv(g.rep[i], t);
   }
}


   
void GCD(ZZX& d, const ZZX& a, const ZZX& b)
{
   if (IsZero(a)) {
      d = b;
      if (sign(LeadCoeff(d)) < 0) negate(d, d);
      return;
   }

   if (IsZero(b)) {
      d = a;
      if (sign(LeadCoeff(d)) < 0) negate(d, d);
      return;
   }

   ZZ c1, c2, c;
   ZZX f1, f2;

   content(c1, a);
   divide(f1, a, c1);

   content(c2, b);
   divide(f2, b, c2);

   GCD(c, c1, c2);

   ZZ ld;
   GCD(ld, LeadCoeff(f1), LeadCoeff(f2));

   ZZX g, h, res;

   ZZ prod;
   set(prod);

   zz_pBak bak;
   bak.save();


   long FirstTime = 1;

   long i;
   for (i = 0; ;i++) {
      zz_p::FFTInit(i);
      long p = zz_p::modulus();

      if (divide(LeadCoeff(f1), p) || divide(LeadCoeff(f2), p)) continue;

      zz_pX G, F1, F2;
      zz_p  LD;

      conv(F1, f1);
      conv(F2, f2);
      conv(LD, ld);

      GCD(G, F1, F2);
      mul(G, G, LD);


      if (deg(G) == 0) { 
         set(res);
         break;
      }

      if (FirstTime || deg(G) < deg(g)) {
         FirstTime = 0;
         conv(prod, p);
         BalCopy(g, G);
      }
      else if (deg(G) > deg(g)) 
         continue;
      else if (!CRT(g, prod, G)) {
         PrimitivePart(res, g);
         if (divide(f1, res) && divide(f2, res))
            break;
      }

   }

   bak.restore();

   mul(d, res, c);
   if (sign(LeadCoeff(d)) < 0) negate(d, d);
}

void trunc(ZZX& x, const ZZX& a, long m)

// x = a % X^m, output may alias input

{
   if (m < 0) Error("trunc: bad args");

   if (&x == &a) {
      if (x.rep.length() > m) {
         x.rep.SetLength(m);
         x.normalize();
      }
   }
   else {
      long n;
      long i;
      ZZ* xp;
      const ZZ* ap;

      n = min(a.rep.length(), m);
      x.rep.SetLength(n);

      xp = x.rep.elts();
      ap = a.rep.elts();

      for (i = 0; i < n; i++) xp[i] = ap[i];

      x.normalize();
   }
}



void LeftShift(ZZX& x, const ZZX& a, long n)
{
   if (IsZero(a)) {
      clear(x);
      return;
   }

   if (n < 0) {
      if (n < -NTL_MAX_LONG)  
         clear(x);
      else
         RightShift(x, a, -n);
      return;
   }

   if (NTL_OVERFLOW(n, 1, 0))
      Error("overflow in LeftShift");

   long m = a.rep.length();

   x.rep.SetLength(m+n);

   long i;
   for (i = m-1; i >= 0; i--)
      x.rep[i+n] = a.rep[i];

   for (i = 0; i < n; i++)
      clear(x.rep[i]);
}


void RightShift(ZZX& x, const ZZX& a, long n)
{
   if (IsZero(a)) {
      clear(x);
      return;
   }

   if (n < 0) {
      if (n < -NTL_MAX_LONG) Error("overflow in RightShift");
      LeftShift(x, a, -n);
      return;
   }

   long da = deg(a);
   long i;

   if (da < n) {
      clear(x);
      return;
   }

   if (&x != &a)
      x.rep.SetLength(da-n+1);

   for (i = 0; i <= da-n; i++)
      x.rep[i] = a.rep[i+n];

   if (&x == &a)
      x.rep.SetLength(da-n+1);

   x.normalize();
}


void TraceVec(vec_ZZ& S, const ZZX& ff)
{
   if (!IsOne(LeadCoeff(ff)))
      Error("TraceVec: bad args");

   ZZX f;
   f = ff;

   long n = deg(f);

   S.SetLength(n);

   if (n == 0)
      return;

   long k, i;
   ZZ acc, t;

   S[0] = n;

   for (k = 1; k < n; k++) {
      mul(acc, f.rep[n-k], k);

      for (i = 1; i < k; i++) {
         mul(t, f.rep[n-i], S[k-i]);
         add(acc, acc, t);
      }

      negate(S[k], acc);
   }

}

static
void EuclLength(ZZ& l, const ZZX& a)
{
   long n = a.rep.length();
   long i;
 
   ZZ sum, t;

   clear(sum);
   for (i = 0; i < n; i++) {
      sqr(t, a.rep[i]);
      add(sum, sum, t);
   }

   if (sum > 1) {
      SqrRoot(l, sum);
      add(l, l, 1);
   }
   else
      l = sum;
}



static
long ResBound(const ZZX& a, const ZZX& b)
{
   if (IsZero(a) || IsZero(b)) 
      return 0;

   ZZ t1, t2, t;
   EuclLength(t1, a);
   EuclLength(t2, b);
   power(t1, t1, deg(b));
   power(t2, t2, deg(a));
   mul(t, t1, t2);
   return NumBits(t);
}



void resultant(ZZ& rres, const ZZX& a, const ZZX& b, long deterministic)
{
   if (IsZero(a) || IsZero(b)) {
      clear(rres);
      return;
   }

   zz_pBak zbak;
   zbak.save();

   ZZ_pBak Zbak;
   Zbak.save();

   long instable = 1;

   long bound = 2+ResBound(a, b);

   long gp_cnt = 0;

   ZZ res, prod;

   clear(res);
   set(prod);


   long i;
   for (i = 0; ; i++) {
      if (NumBits(prod) > bound)
         break;

      if (!deterministic &&
          !instable && bound > 1000 && NumBits(prod) < 0.25*bound) {

         ZZ P;


         long plen = 90 + NumBits(max(bound, NumBits(res)));

         do {
            GenPrime(P, plen, 90 + 2*NumBits(gp_cnt++));
         }
         while (divide(LeadCoeff(a), P) || divide(LeadCoeff(b), P));

         ZZ_p::init(P);

         ZZ_pX A, B;
         conv(A, a);
         conv(B, b);

         ZZ_p t;
         resultant(t, A, B);

         if (CRT(res, prod, rep(t), P))
            instable = 1;
         else
            break;
      }


      zz_p::FFTInit(i);
      long p = zz_p::modulus();
      if (divide(LeadCoeff(a), p) || divide(LeadCoeff(b), p))
         continue;

      zz_pX A, B;
      conv(A, a);
      conv(B, b);

      zz_p t;
      resultant(t, A, B);

      instable = CRT(res, prod, rep(t), p);
   }

   rres = res;

   zbak.restore();
   Zbak.restore();
}




void MinPolyMod(ZZX& gg, const ZZX& a, const ZZX& f)

{
   if (!IsOne(LeadCoeff(f)) || deg(f) < 1 || deg(a) >= deg(f))
      Error("MinPolyMod: bad args");

   if (IsZero(a)) {
      SetX(gg);
      return;
   }

   ZZ_pBak Zbak;
   Zbak.save();
   zz_pBak zbak;
   zbak.save();

   long n = deg(f);

   long instable = 1;

   long gp_cnt = 0;

   ZZ prod;
   ZZX g;

   clear(g);
   set(prod);

   long bound = -1;

   long i;
   for (i = 0; ; i++) {
      if (deg(g) == n) {
         if (bound < 0)
            bound = 2+CharPolyBound(a, f);

         if (NumBits(prod) > bound)
            break;
      }

      if (!instable && 
         (deg(g) < n || 
         (deg(g) == n && bound > 1000 && NumBits(prod) < 0.75*bound))) {

         // guarantees 2^{-80} error probability
         long plen = 90 + max( 2*NumBits(n) + NumBits(MaxBits(f)),
                         max( NumBits(n) + NumBits(MaxBits(a)),
                              NumBits(MaxBits(g)) ));

         ZZ P;
         GenPrime(P, plen, 90 + 2*NumBits(gp_cnt++));
         ZZ_p::init(P);


         ZZ_pX A, F, G;
         conv(A, a);
         conv(F, f);
         conv(G, g);

         ZZ_pXModulus FF;
         build(FF, F);

         ZZ_pX H;
         CompMod(H, G, A, FF);
         
         if (IsZero(H))
            break;

         instable = 1;
      } 
         
      zz_p::FFTInit(i);

      zz_pX A, F;
      conv(A, a);
      conv(F, f);

      zz_pXModulus FF;
      build(FF, F);

      zz_pX G;
      MinPolyMod(G, A, FF);

      if (deg(G) < deg(g))
         continue;

      if (deg(G) > deg(g)) {
         clear(g);
         set(prod);
      }

      instable = CRT(g, prod, G);
   }

   gg = g;

   Zbak.restore();
   zbak.restore();
}


void XGCD(ZZ& rr, ZZX& ss, ZZX& tt, const ZZX& a, const ZZX& b, 
          long deterministic)
{
   ZZ r;

   resultant(r, a, b, deterministic);

   if (IsZero(r)) {
      clear(rr);
      return;
   }

   zz_pBak bak;
   bak.save();

   long i;
   long instable = 1;

   ZZ tmp;
   ZZ prod;
   ZZX s, t;

   set(prod);
   clear(s);
   clear(t);

   for (i = 0; ; i++) {
      zz_p::FFTInit(i);
      long p = zz_p::modulus();

      if (divide(LeadCoeff(a), p) || divide(LeadCoeff(b), p) || divide(r, p))
         continue;

      zz_p R;
      conv(R, r);

      zz_pX D, S, T, A, B;
      conv(A, a);
      conv(B, b);

      if (!instable) {
         conv(S, s);
         conv(T, t);
         zz_pX t1, t2;
         mul(t1, A, S); 
         mul(t2, B, T);
         add(t1, t1, t2);

         if (deg(t1) == 0 && ConstTerm(t1) == R)
            mul(prod, prod, p);
         else
            instable = 1;
      }

      if (instable) {
         XGCD(D, S, T, A, B);
   
         mul(S, S, R);
         mul(T, T, R);
   
         tmp = prod;
         long Sinstable = CRT(s, tmp, S);
         long Tinstable = CRT(t, prod, T);
   
         instable = Sinstable || Tinstable;
      }

      if (!instable) {
         long bound1 = NumBits(min(deg(a), deg(s)) + 1) 
                      + MaxBits(a) + MaxBits(s);
         long bound2 = NumBits(min(deg(b), deg(t)) + 1) 
                      + MaxBits(b) + MaxBits(t);

         long bound = 4 + max(NumBits(r), max(bound1, bound2));

         if (NumBits(prod) > bound)
            break;
      }
   }

   rr = r;
   ss = s;
   tt = t;

   bak.restore();
}

void NormMod(ZZ& x, const ZZX& a, const ZZX& f, long deterministic)
{
   if (!IsOne(LeadCoeff(f)) || deg(a) >= deg(f) || deg(f) <= 0)
      Error("norm: bad args");

   if (IsZero(a)) {
      clear(x);
      return;
   }

   resultant(x, f, a, deterministic);
}

void TraceMod(ZZ& res, const ZZX& a, const ZZX& f)
{
   if (!IsOne(LeadCoeff(f)) || deg(a) >= deg(f) || deg(f) <= 0)
      Error("trace: bad args");

   vec_ZZ S;

   TraceVec(S, f);

   InnerProduct(res, S, a.rep);
}


void discriminant(ZZ& d, const ZZX& a, long deterministic)
{
   long m = deg(a);

   if (m < 0) {
      clear(d);
      return;
   }

   ZZX a1;
   ZZ res;

   diff(a1, a);
   resultant(res, a, a1, deterministic);
   if (!divide(res, res, LeadCoeff(a)))
      Error("discriminant: inexact division");

   m = m & 3;
   if (m >= 2)
      negate(res, res);

   d = res;
}


void MulMod(ZZX& x, const ZZX& a, const ZZX& b, const ZZX& f)
{
   if (deg(a) >= deg(f) || deg(b) >= deg(f) || deg(f) == 0 || 
       !IsOne(LeadCoeff(f)))
      Error("MulMod: bad args");

   ZZX t;
   mul(t, a, b);
   rem(x, t, f);
}

void SqrMod(ZZX& x, const ZZX& a, const ZZX& f)
{
   if (deg(a) >= deg(f) || deg(f) == 0 || !IsOne(LeadCoeff(f)))
      Error("MulMod: bad args");

   ZZX t;
   sqr(t, a);
   rem(x, t, f);
}



static
void MulByXModAux(ZZX& h, const ZZX& a, const ZZX& f)
{
   long i, n, m;
   ZZ* hh;
   const ZZ *aa, *ff;

   ZZ t, z;


   n = deg(f);
   m = deg(a);

   if (m >= n || n == 0 || !IsOne(LeadCoeff(f)))
      Error("MulByXMod: bad args");

   if (m < 0) {
      clear(h);
      return;
   }

   if (m < n-1) {
      h.rep.SetLength(m+2);
      hh = h.rep.elts();
      aa = a.rep.elts();
      for (i = m+1; i >= 1; i--)
         hh[i] = aa[i-1];
      clear(hh[0]);
   }
   else {
      h.rep.SetLength(n);
      hh = h.rep.elts();
      aa = a.rep.elts();
      ff = f.rep.elts();
      negate(z, aa[n-1]);
      for (i = n-1; i >= 1; i--) {
         mul(t, z, ff[i]);
         add(hh[i], aa[i-1], t);
      }
      mul(hh[0], z, ff[0]);
      h.normalize();
   }
}

void MulByXMod(ZZX& h, const ZZX& a, const ZZX& f)
{
   if (&h == &f) {
      ZZX hh;
      MulByXModAux(hh, a, f);
      h = hh;
   }
   else
      MulByXModAux(h, a, f);
}

static
void EuclLength1(ZZ& l, const ZZX& a)
{
   long n = a.rep.length();
   long i;
 
   ZZ sum, t;

   clear(sum);
   for (i = 0; i < n; i++) {
      sqr(t, a.rep[i]);
      add(sum, sum, t);
   }

   abs(t, ConstTerm(a));
   mul(t, t, 2);
   add(t, t, 1);
   add(sum, sum, t);

   if (sum > 1) {
      SqrRoot(l, sum);
      add(l, l, 1);
   }
   else
      l = sum;
}


long CharPolyBound(const ZZX& a, const ZZX& f)
// This computes a bound on the size of the
// coefficients of the characterstic polynomial.
// It uses the characterization of the char poly as
// resultant_y(f(y), x-a(y)), and then interpolates this
// through complex primimitive (deg(f)+1)-roots of unity.

{
   if (IsZero(a) || IsZero(f))
      Error("CharPolyBound: bad args");

   ZZ t1, t2, t;
   EuclLength1(t1, a);
   EuclLength(t2, f);
   power(t1, t1, deg(f));
   power(t2, t2, deg(a));
   mul(t, t1, t2);
   return NumBits(t);
}


void SetCoeff(ZZX& x, long i, long a)
{
   if (a == 1) 
      SetCoeff(x, i);
   else {
      static ZZ aa;
      conv(aa, a);
      SetCoeff(x, i, aa);
   }
}


void CopyReverse(ZZX& x, const ZZX& a, long hi)

   // x[0..hi] = reverse(a[0..hi]), with zero fill
   // input may not alias output

{
   long i, j, n, m;

   n = hi+1;
   m = a.rep.length();

   x.rep.SetLength(n);

   const ZZ* ap = a.rep.elts();
   ZZ* xp = x.rep.elts();

   for (i = 0; i < n; i++) {
      j = hi-i;
      if (j < 0 || j >= m)
         clear(xp[i]);
      else
         xp[i] = ap[j];
   }

   x.normalize();
}

void reverse(ZZX& x, const ZZX& a, long hi)
{
   if (hi < 0) { clear(x); return; }
   if (NTL_OVERFLOW(hi, 1, 0))
      Error("overflow in reverse");

   if (&x == &a) {
      ZZX tmp;
      CopyReverse(tmp, a, hi);
      x = tmp;
   }
   else
      CopyReverse(x, a, hi);
}

void MulTrunc(ZZX& x, const ZZX& a, const ZZX& b, long n)
{
   ZZX t;
   mul(t, a, b);
   trunc(x, t, n);
}

void SqrTrunc(ZZX& x, const ZZX& a, long n)
{
   ZZX t;
   sqr(t, a);
   trunc(x, t, n);
}


void NewtonInvTrunc(ZZX& c, const ZZX& a, long e)
{
   ZZ x;

   if (ConstTerm(a) == 1)
      x = 1;
   else if (ConstTerm(a) == -1)
      x = -1;
   else
      Error("InvTrunc: non-invertible constant term");

   if (e == 1) {
      conv(c, x);
      return;
   }

   static vec_long E;
   E.SetLength(0);
   append(E, e);
   while (e > 1) {
      e = (e+1)/2;
      append(E, e);
   }

   long L = E.length();

   ZZX g, g0, g1, g2;


   g.rep.SetMaxLength(E[0]);
   g0.rep.SetMaxLength(E[0]);
   g1.rep.SetMaxLength((3*E[0]+1)/2);
   g2.rep.SetMaxLength(E[0]);

   conv(g, x);

   long i;

   for (i = L-1; i > 0; i--) {
      // lift from E[i] to E[i-1]

      long k = E[i];
      long l = E[i-1]-E[i];

      trunc(g0, a, k+l);

      mul(g1, g0, g);
      RightShift(g1, g1, k);
      trunc(g1, g1, l);

      mul(g2, g1, g);
      trunc(g2, g2, l);
      LeftShift(g2, g2, k);

      sub(g, g, g2);
   }

   c = g;
}


void InvTrunc(ZZX& c, const ZZX& a, long e)
{
   if (e < 0) Error("InvTrunc: bad args");

   if (e == 0) {
      clear(c);
      return;
   }

   if (NTL_OVERFLOW(e, 1, 0))
      Error("overflow in InvTrunc");

   NewtonInvTrunc(c, a, e);
}

NTL_END_IMPL
