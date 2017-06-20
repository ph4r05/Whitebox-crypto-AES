
#include <NTL/ZZ_pX.h>


// The mul & sqr routines use routines from ZZX, 
// which is faster for small degree polynomials.
// Define this macro to revert to old strategy.


#ifndef NTL_WIZARD_HACK

#include <NTL/ZZX.h>

#endif

#include <NTL/new.h>


#if (defined(NTL_GMP_LIP) || defined(NTL_GMP_HACK))
#define KARX 200
#else
#define KARX 80
#endif


NTL_START_IMPL




const ZZ_pX& ZZ_pX::zero()
{
   static ZZ_pX z;
   return z;
}


ZZ_pX& ZZ_pX::operator=(long a)
{
   conv(*this, a);
   return *this;
}


ZZ_pX& ZZ_pX::operator=(const ZZ_p& a)
{
   conv(*this, a);
   return *this;
}


istream& operator>>(istream& s, ZZ_pX& x)
{
   s >> x.rep;
   x.normalize();
   return s;
}

ostream& operator<<(ostream& s, const ZZ_pX& a)
{
   return s << a.rep;
}


void ZZ_pX::normalize()
{
   long n;
   const ZZ_p* p;

   n = rep.length();
   if (n == 0) return;
   p = rep.elts() + n;
   while (n > 0 && IsZero(*--p)) {
      n--;
   }
   rep.SetLength(n);
}


long IsZero(const ZZ_pX& a)
{
   return a.rep.length() == 0;
}


long IsOne(const ZZ_pX& a)
{
    return a.rep.length() == 1 && IsOne(a.rep[0]);
}

void GetCoeff(ZZ_p& x, const ZZ_pX& a, long i)
{
   if (i < 0 || i > deg(a))
      clear(x);
   else
      x = a.rep[i];
}

void SetCoeff(ZZ_pX& x, long i, const ZZ_p& a)
{
   long j, m;

   if (i < 0) 
      Error("SetCoeff: negative index");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m && IsZero(a)) return; 

   if (i > m) {
      /* careful: a may alias a coefficient of x */

      long alloc = x.rep.allocated();

      if (alloc > 0 && i >= alloc) {
         ZZ_pTemp aa_tmp;  ZZ_p& aa = aa_tmp.val();
         aa = a;
         x.rep.SetLength(i+1);
         x.rep[i] = aa;
      }
      else {
         x.rep.SetLength(i+1);
         x.rep[i] = a;
      }

      for (j = m+1; j < i; j++)
         clear(x.rep[j]);
   }
   else
      x.rep[i] = a;

   x.normalize();
}

void SetCoeff(ZZ_pX& x, long i, long a)
{
   if (a == 1) 
      SetCoeff(x, i);
   else {
      ZZ_pTemp TT;  ZZ_p& T = TT.val(); 
      conv(T, a);
      SetCoeff(x, i, T);
   }
}

void SetCoeff(ZZ_pX& x, long i)
{
   long j, m;

   if (i < 0) 
      Error("coefficient index out of range");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      x.rep.SetLength(i+1);
      for (j = m+1; j < i; j++)
         clear(x.rep[j]);
   }
   set(x.rep[i]);
   x.normalize();
}


void SetX(ZZ_pX& x)
{
   clear(x);
   SetCoeff(x, 1);
}


long IsX(const ZZ_pX& a)
{
   return deg(a) == 1 && IsOne(LeadCoeff(a)) && IsZero(ConstTerm(a));
}
      
      

const ZZ_p& coeff(const ZZ_pX& a, long i)
{
   if (i < 0 || i > deg(a))
      return ZZ_p::zero();
   else
      return a.rep[i];
}


const ZZ_p& LeadCoeff(const ZZ_pX& a)
{
   if (IsZero(a))
      return ZZ_p::zero();
   else
      return a.rep[deg(a)];
}

const ZZ_p& ConstTerm(const ZZ_pX& a)
{
   if (IsZero(a))
      return ZZ_p::zero();
   else
      return a.rep[0];
}



void conv(ZZ_pX& x, const ZZ_p& a)
{
   if (IsZero(a))
      x.rep.SetLength(0);
   else {
      x.rep.SetLength(1);
      x.rep[0] = a;

      // note: if a aliases x.rep[i], i > 0, this code
      //       will still work, since is is assumed that
      //       SetLength(1) will not relocate or destroy x.rep[i]
   }
}

void conv(ZZ_pX& x, long a)
{
   if (a == 0)
      clear(x);
   else if (a == 1)
      set(x);
   else {
      ZZ_pTemp TT; ZZ_p& T = TT.val();
      conv(T, a);
      conv(x, T);
   }
}

void conv(ZZ_pX& x, const ZZ& a)
{
   if (IsZero(a))
      clear(x);
   else {
      ZZ_pTemp TT; ZZ_p& T = TT.val();
      conv(T, a);
      conv(x, T);
   }
}

void conv(ZZ_pX& x, const vec_ZZ_p& a)
{
   x.rep = a;
   x.normalize();
}


void add(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b)
{
   long da = deg(a);
   long db = deg(b);
   long minab = min(da, db);
   long maxab = max(da, db);
   x.rep.SetLength(maxab+1);

   long i;
   const ZZ_p *ap, *bp; 
   ZZ_p* xp;

   for (i = minab+1, ap = a.rep.elts(), bp = b.rep.elts(), xp = x.rep.elts();
        i; i--, ap++, bp++, xp++)
      add(*xp, (*ap), (*bp));

   if (da > minab && &x != &a)
      for (i = da-minab; i; i--, xp++, ap++)
         *xp = *ap;
   else if (db > minab && &x != &b)
      for (i = db-minab; i; i--, xp++, bp++)
         *xp = *bp;
   else
      x.normalize();
}


void add(ZZ_pX& x, const ZZ_pX& a, const ZZ_p& b)
{
   long n = a.rep.length();
   if (n == 0) {
      conv(x, b);
   }
   else if (&x == &a) {
      add(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else if (x.rep.MaxLength() == 0) {
      x = a;
      add(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else {
      // ugly...b could alias a coeff of x

      ZZ_p *xp = x.rep.elts();
      add(xp[0], a.rep[0], b);
      x.rep.SetLength(n);
      xp = x.rep.elts();
      const ZZ_p *ap = a.rep.elts();
      long i;
      for (i = 1; i < n; i++)
         xp[i] = ap[i];
      x.normalize();
   }
}

void add(ZZ_pX& x, const ZZ_pX& a, long b)
{
   if (a.rep.length() == 0) {
      conv(x, b);
   }
   else {
      if (&x != &a) x = a;
      add(x.rep[0], x.rep[0], b);
      x.normalize();
   }
}

void sub(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b)
{
   long da = deg(a);
   long db = deg(b);
   long minab = min(da, db);
   long maxab = max(da, db);
   x.rep.SetLength(maxab+1);

   long i;
   const ZZ_p *ap, *bp; 
   ZZ_p* xp;

   for (i = minab+1, ap = a.rep.elts(), bp = b.rep.elts(), xp = x.rep.elts();
        i; i--, ap++, bp++, xp++)
      sub(*xp, (*ap), (*bp));

   if (da > minab && &x != &a)
      for (i = da-minab; i; i--, xp++, ap++)
         *xp = *ap;
   else if (db > minab)
      for (i = db-minab; i; i--, xp++, bp++)
         negate(*xp, *bp);
   else
      x.normalize();

}

void sub(ZZ_pX& x, const ZZ_pX& a, const ZZ_p& b)
{
   long n = a.rep.length();
   if (n == 0) {
      conv(x, b);
      negate(x, x);
   }
   else if (&x == &a) {
      sub(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else if (x.rep.MaxLength() == 0) {
      x = a;
      sub(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else {
      // ugly...b could alias a coeff of x

      ZZ_p *xp = x.rep.elts();
      sub(xp[0], a.rep[0], b);
      x.rep.SetLength(n);
      xp = x.rep.elts();
      const ZZ_p *ap = a.rep.elts();
      long i;
      for (i = 1; i < n; i++)
         xp[i] = ap[i];
      x.normalize();
   }
}

void sub(ZZ_pX& x, const ZZ_pX& a, long b)
{
   if (b == 0) {
      x = a;
      return;
   }

   if (a.rep.length() == 0) {
      x.rep.SetLength(1);
      x.rep[0] = b;
      negate(x.rep[0], x.rep[0]);
   }
   else {
      if (&x != &a) x = a;
      sub(x.rep[0], x.rep[0], b);
   }
   x.normalize();
}

void sub(ZZ_pX& x, const ZZ_p& a, const ZZ_pX& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val(); 
   T = a;

   negate(x, b);
   add(x, x, T);
}

void sub(ZZ_pX& x, long a, const ZZ_pX& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val(); 
   T = a;

   negate(x, b);
   add(x, x, T);
}

void negate(ZZ_pX& x, const ZZ_pX& a)
{
   long n = a.rep.length();
   x.rep.SetLength(n);

   const ZZ_p* ap = a.rep.elts();
   ZZ_p* xp = x.rep.elts();
   long i;

   for (i = n; i; i--, ap++, xp++)
      negate((*xp), (*ap));

}


#ifndef NTL_WIZARD_HACK

// These crossovers are tuned for a Pentium, but hopefully
// they should be OK on other machines as well.


const long SS_kbound = 40;
const double SS_rbound = 1.25;


void mul(ZZ_pX& c, const ZZ_pX& a, const ZZ_pX& b)
{
   if (IsZero(a) || IsZero(b)) {
      clear(c);
      return;
   }

   if (&a == &b) {
      sqr(c, a);
      return;
   }

   long k = ZZ_p::ModulusSize();
   long s = min(deg(a), deg(b)) + 1;

   if (s == 1 || (k == 1 && s < 40) || (k == 2 && s < 20) ||
                 (k == 3 && s < 12) || (k <= 5 && s < 8) ||
                 (k <= 12 && s < 4) )  {
      PlainMul(c, a, b);
   }
   else if (s < KARX) {
      ZZX A, B, C;
      conv(A, a);
      conv(B, b);
      KarMul(C, A, B);
      conv(c, C);
   }
   else {
      long mbits;
      mbits = NumBits(ZZ_p::modulus());
      if (k >= SS_kbound && 
          SSRatio(deg(a), mbits, deg(b), mbits) < SS_rbound) {
         ZZX A, B, C;
         conv(A, a);
         conv(B, b);
         SSMul(C, A, B);
         conv(c, C);
      }
      else {
         FFTMul(c, a, b);
      }
   }
}

void sqr(ZZ_pX& c, const ZZ_pX& a)
{
   if (IsZero(a)) {
      clear(c);
      return;
   }

   long k = ZZ_p::ModulusSize();
   long s = deg(a) + 1;

   if (s == 1 || (k == 1 && s < 50) || (k == 2 && s < 25) ||
                 (k == 3 && s < 25) || (k <= 6 && s < 12) ||
                 (k <= 8 && s < 8)  || (k == 9 && s < 6)  ||
                 (k <= 30 && s < 4) ) {

      PlainSqr(c, a);
   }
   else if (s < 80) {
      ZZX C, A;
      conv(A, a);
      KarSqr(C, A);
      conv(c, C);
   }
   else {
      long mbits;
      mbits = NumBits(ZZ_p::modulus());
      if (k >= SS_kbound && 
          SSRatio(deg(a), mbits, deg(a), mbits) < SS_rbound) {
         ZZX A, C;
         conv(A, a);
         SSSqr(C, A);
         conv(c, C);
      }
      else {
         FFTSqr(c, a);
      }
   }
}

#else

void mul(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b)
{
   if (&a == &b) {
      sqr(x, a);
      return;
   }

   if (deg(a) > NTL_ZZ_pX_FFT_CROSSOVER && deg(b) > NTL_ZZ_pX_FFT_CROSSOVER)
      FFTMul(x, a, b);
   else
      PlainMul(x, a, b);
}

void sqr(ZZ_pX& x, const ZZ_pX& a)
{
   if (deg(a) > NTL_ZZ_pX_FFT_CROSSOVER)
      FFTSqr(x, a);
   else
      PlainSqr(x, a);
}


#endif


void PlainMul(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b)
{
   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   if (da == 0) {
      mul(x, b, a.rep[0]);
      return;
   }

   if (db == 0) {
      mul(x, a, b.rep[0]);
      return;
   }

   long d = da+db;



   const ZZ_p *ap, *bp;
   ZZ_p *xp;
   
   ZZ_pX la, lb;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();

   if (&x == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   static ZZ t, accum;

   for (i = 0; i <= d; i++) {
      jmin = max(0, i-db);
      jmax = min(da, i);
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(bp[i-j]));
	 add(accum, accum, t);
      }
      conv(xp[i], accum);
   }
   x.normalize();
}

void PlainSqr(ZZ_pX& x, const ZZ_pX& a)
{
   long da = deg(a);

   if (da < 0) {
      clear(x);
      return;
   }

   long d = 2*da;

   const ZZ_p *ap;
   ZZ_p *xp;

   ZZ_pX la;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();


   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   long m, m2;
   static ZZ t, accum;

   for (i = 0; i <= d; i++) {
      jmin = max(0, i-da);
      jmax = min(da, i);
      m = jmax - jmin + 1;
      m2 = m >> 1;
      jmax = jmin + m2 - 1;
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(ap[i-j]));
	 add(accum, accum, t);
      }
      add(accum, accum, accum);
      if (m & 1) {
	 sqr(t, rep(ap[jmax + 1]));
	 add(accum, accum, t);
      }

      conv(xp[i], accum);
   }

   x.normalize();
}

void PlainDivRem(ZZ_pX& q, ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ_p *bp;
   ZZ_p *qp;
   ZZ *xp;


   ZZ_p LCInv, t;
   static ZZ s;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZ_pX: division by zero");

   if (da < db) {
      r = a;
      clear(q);
      return;
   }

   ZZ_pX lb;

   if (&q == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      inv(LCInv, bp[db]);
   }

   ZZVec x(da + 1, ZZ_pInfo->ExtendedModulusSize);

   for (i = 0; i <= da; i++)
      x[i] = rep(a.rep[i]);

   xp = x.elts();

   dq = da - db;
   q.rep.SetLength(dq+1);
   qp = q.rep.elts();

   for (i = dq; i >= 0; i--) {
      conv(t, xp[i+db]);
      if (!LCIsOne)
	 mul(t, t, LCInv);
      qp[i] = t;
      negate(t, t);

      for (j = db-1; j >= 0; j--) {
	 mul(s, rep(t), rep(bp[j]));
	 add(xp[i+j], xp[i+j], s);
      }
   }

   r.rep.SetLength(db);
   for (i = 0; i < db; i++)
      conv(r.rep[i], xp[i]);
   r.normalize();
}


void PlainRem(ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b, ZZVec& x)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ_p *bp;
   ZZ *xp;


   ZZ_p LCInv, t;
   static ZZ s;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZ_pX: division by zero");

   if (da < db) {
      r = a;
      return;
   }

   bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      inv(LCInv, bp[db]);
   }

   for (i = 0; i <= da; i++)
      x[i] = rep(a.rep[i]);

   xp = x.elts();

   dq = da - db;

   for (i = dq; i >= 0; i--) {
      conv(t, xp[i+db]);
      if (!LCIsOne)
	 mul(t, t, LCInv);
      negate(t, t);

      for (j = db-1; j >= 0; j--) {
	 mul(s, rep(t), rep(bp[j]));
	 add(xp[i+j], xp[i+j], s);
      }
   }

   r.rep.SetLength(db);
   for (i = 0; i < db; i++)
      conv(r.rep[i], xp[i]);
   r.normalize();
}


void PlainDivRem(ZZ_pX& q, ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b, ZZVec& x)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ_p *bp;
   ZZ_p *qp;
   ZZ *xp;


   ZZ_p LCInv, t;
   static ZZ s;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZ_pX: division by zero");

   if (da < db) {
      r = a;
      clear(q);
      return;
   }

   ZZ_pX lb;

   if (&q == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      inv(LCInv, bp[db]);
   }

   for (i = 0; i <= da; i++)
      x[i] = rep(a.rep[i]);

   xp = x.elts();

   dq = da - db;
   q.rep.SetLength(dq+1);
   qp = q.rep.elts();

   for (i = dq; i >= 0; i--) {
      conv(t, xp[i+db]);
      if (!LCIsOne)
	 mul(t, t, LCInv);
      qp[i] = t;
      negate(t, t);

      for (j = db-1; j >= 0; j--) {
	 mul(s, rep(t), rep(bp[j]));
	 add(xp[i+j], xp[i+j], s);
      }
   }

   r.rep.SetLength(db);
   for (i = 0; i < db; i++)
      conv(r.rep[i], xp[i]);
   r.normalize();
}


void PlainDiv(ZZ_pX& q, const ZZ_pX& a, const ZZ_pX& b)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ_p *bp;
   ZZ_p *qp;
   ZZ *xp;


   ZZ_p LCInv, t;
   static ZZ s;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZ_pX: division by zero");

   if (da < db) {
      clear(q);
      return;
   }

   ZZ_pX lb;

   if (&q == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      inv(LCInv, bp[db]);
   }

   ZZVec x(da + 1 - db, ZZ_pInfo->ExtendedModulusSize);

   for (i = db; i <= da; i++)
      x[i-db] = rep(a.rep[i]);

   xp = x.elts();

   dq = da - db;
   q.rep.SetLength(dq+1);
   qp = q.rep.elts();

   for (i = dq; i >= 0; i--) {
      conv(t, xp[i]);
      if (!LCIsOne)
	 mul(t, t, LCInv);
      qp[i] = t;
      negate(t, t);

      long lastj = max(0, db-i);

      for (j = db-1; j >= lastj; j--) {
	 mul(s, rep(t), rep(bp[j]));
	 add(xp[i+j-db], xp[i+j-db], s);
      }
   }
}

void PlainRem(ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b)
{
   long da, db, dq, i, j, LCIsOne;
   const ZZ_p *bp;
   ZZ *xp;


   ZZ_p LCInv, t;
   static ZZ s;

   da = deg(a);
   db = deg(b);

   if (db < 0) Error("ZZ_pX: division by zero");

   if (da < db) {
      r = a;
      return;
   }

   bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      inv(LCInv, bp[db]);
   }

   ZZVec x(da + 1, ZZ_pInfo->ExtendedModulusSize);

   for (i = 0; i <= da; i++)
      x[i] = rep(a.rep[i]);

   xp = x.elts();

   dq = da - db;

   for (i = dq; i >= 0; i--) {
      conv(t, xp[i+db]);
      if (!LCIsOne)
	 mul(t, t, LCInv);
      negate(t, t);

      for (j = db-1; j >= 0; j--) {
	 mul(s, rep(t), rep(bp[j]));
	 add(xp[i+j], xp[i+j], s);
      }
   }

   r.rep.SetLength(db);
   for (i = 0; i < db; i++)
      conv(r.rep[i], xp[i]);
   r.normalize();
}

void mul(ZZ_pX& x, const ZZ_pX& a, const ZZ_p& b)
{
   if (IsZero(b)) {
      clear(x);
      return;
   }

   if (IsOne(b)) {
      x = a;
      return;
   }

   ZZ_pTemp TT; ZZ_p& t = TT.val();

   long i, da;

   const ZZ_p *ap;
   ZZ_p* xp;


   t = b;

   da = deg(a);
   x.rep.SetLength(da+1);
   ap = a.rep.elts();
   xp = x.rep.elts();

   for (i = 0; i <= da; i++) 
      mul(xp[i], ap[i], t);

   x.normalize();
}

void mul(ZZ_pX& x, const ZZ_pX& a, long b)
{
   ZZ_pTemp TT;  ZZ_p& T = TT.val();
   conv(T, b);
   mul(x, a, T);
}


void PlainGCD(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b)
{
   ZZ_p t;

   if (IsZero(b))
      x = a;
   else if (IsZero(a))
      x = b;
   else {
      long n = max(deg(a),deg(b)) + 1;
      ZZ_pX u(INIT_SIZE, n), v(INIT_SIZE, n);
      ZZVec tmp(n, ZZ_pInfo->ExtendedModulusSize);

      u = a;
      v = b;
      do {
         PlainRem(u, u, v, tmp);
         swap(u, v);
      } while (!IsZero(v));

      x = u;
   }

   if (IsZero(x)) return;
   if (IsOne(LeadCoeff(x))) return;

   /* make gcd monic */


   inv(t, LeadCoeff(x)); 
   mul(x, x, t); 
}



         

void PlainXGCD(ZZ_pX& d, ZZ_pX& s, ZZ_pX& t, const ZZ_pX& a, const ZZ_pX& b)
{
   ZZ_p z;


   if (IsZero(b)) {
      set(s);
      clear(t);
      d = a;
   }
   else if (IsZero(a)) {
      clear(s);
      set(t);
      d = b;
   }
   else {
      long e = max(deg(a), deg(b)) + 1;

      ZZ_pX temp(INIT_SIZE, e), u(INIT_SIZE, e), v(INIT_SIZE, e), 
            u0(INIT_SIZE, e), v0(INIT_SIZE, e), 
            u1(INIT_SIZE, e), v1(INIT_SIZE, e), 
            u2(INIT_SIZE, e), v2(INIT_SIZE, e), q(INIT_SIZE, e);


      set(u1); clear(v1);
      clear(u2); set(v2);
      u = a; v = b;

      do {
         DivRem(q, u, u, v);
         swap(u, v);
         u0 = u2;
         v0 = v2;
         mul(temp, q, u2);
         sub(u2, u1, temp);
         mul(temp, q, v2);
         sub(v2, v1, temp);
         u1 = u0;
         v1 = v0;
      } while (!IsZero(v));

      d = u;
      s = u1;
      t = v1;
   }

   if (IsZero(d)) return;
   if (IsOne(LeadCoeff(d))) return;

   /* make gcd monic */

   inv(z, LeadCoeff(d));
   mul(d, d, z);
   mul(s, s, z);
   mul(t, t, z);
}


void MulMod(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b, const ZZ_pX& f)
{
   if (deg(a) >= deg(f) || deg(b) >= deg(f) || deg(f) == 0) 
      Error("MulMod: bad args");

   ZZ_pX t;

   mul(t, a, b);
   rem(x, t, f);
}

void SqrMod(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& f)
{
   if (deg(a) >= deg(f) || deg(f) == 0) Error("SqrMod: bad args");

   ZZ_pX t;

   sqr(t, a);
   rem(x, t, f);
}


void InvMod(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& f)
{
   if (deg(a) >= deg(f) || deg(f) == 0) Error("InvMod: bad args");

   ZZ_pX d, t;

   XGCD(d, x, t, a, f);
   if (!IsOne(d))
      Error("ZZ_pX InvMod: can't compute multiplicative inverse");
}

long InvModStatus(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& f)
{
   if (deg(a) >= deg(f) || deg(f) == 0) Error("InvModStatus: bad args");
   ZZ_pX d, t;

   XGCD(d, x, t, a, f);
   if (!IsOne(d)) {
      x = d;
      return 1;
   }
   else
      return 0;
}




static
void MulByXModAux(ZZ_pX& h, const ZZ_pX& a, const ZZ_pX& f)
{
   long i, n, m;
   ZZ_p* hh;
   const ZZ_p *aa, *ff;

   ZZ_p t, z;

   n = deg(f);
   m = deg(a);

   if (m >= n || n == 0) Error("MulByXMod: bad args");

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
      if (!IsOne(ff[n]))
         div(z, z, ff[n]);
      for (i = n-1; i >= 1; i--) {
         mul(t, z, ff[i]);
         add(hh[i], aa[i-1], t);
      }
      mul(hh[0], z, ff[0]);
      h.normalize();
   }
}


void MulByXMod(ZZ_pX& h, const ZZ_pX& a, const ZZ_pX& f)
{
   if (&h == &f) {
      ZZ_pX hh;
      MulByXModAux(hh, a, f);
      h = hh;
   }
   else
      MulByXModAux(h, a, f);
}



void random(ZZ_pX& x, long n)
{
   long i;

   x.rep.SetLength(n);

   for (i = 0; i < n; i++)
      random(x.rep[i]); 

   x.normalize();
}


void FFTRep::SetSize(long NewK)
{

   if (NewK < -1 || NewK >= NTL_BITS_PER_LONG-1)
      Error("bad arg to FFTRep::SetSize()");

   if (NewK <= MaxK) {
      k = NewK;
      return;
   }

   ZZ_pInfo->check();

   if (MaxK == -1)
      NumPrimes = ZZ_pInfo->NumPrimes;
   else {
      if (NumPrimes != ZZ_pInfo->NumPrimes)
         Error("FFTRep: inconsistent use");
   }

   long i, n;

   if (MaxK == -1) {
      tbl = (long **) NTL_MALLOC(NumPrimes, sizeof(long *), 0);
      if (!tbl)
         Error("out of space in FFTRep::SetSize()");
   }
   else {
      for (i = 0; i < NumPrimes; i++) 
         free(tbl[i]);
   }

   n = 1L << NewK;

   for (i = 0; i < NumPrimes; i++) {
      if ( !(tbl[i] = (long *) NTL_MALLOC(n, sizeof(long), 0)) )
         Error("out of space in FFTRep::SetSize()");
   }

   k = MaxK = NewK;
}

FFTRep::FFTRep(const FFTRep& R)
{
   k = MaxK = R.k;
   tbl = 0;
   NumPrimes = 0;

   if (k < 0) return;

   NumPrimes = R.NumPrimes;

   long i, j, n;
 
   tbl = (long **) NTL_MALLOC(NumPrimes, sizeof(long *), 0);
   if (!tbl)
      Error("out of space in FFTRep");

   n = 1L << k;

   for (i = 0; i < NumPrimes; i++) {
      if ( !(tbl[i] = (long *) NTL_MALLOC(n, sizeof(long), 0)) )
         Error("out of space in FFTRep");

      for (j = 0; j < n; j++)
         tbl[i][j] = R.tbl[i][j];
   }
}

FFTRep& FFTRep::operator=(const FFTRep& R)
{
   if (this == &R) return *this;

   if (MaxK >= 0 && R.MaxK >= 0 && NumPrimes != R.NumPrimes)
      Error("FFTRep: inconsistent use");

   if (R.k < 0) {
      k = -1;
      return *this;
   }

   NumPrimes = R.NumPrimes;

   if (R.k > MaxK) {
      long i, n;

      if (MaxK == -1) {
         tbl = (long **) NTL_MALLOC(NumPrimes, sizeof(long *), 0);
         if (!tbl)
            Error("out of space in FFTRep");
      }
      else {
         for (i = 0; i < NumPrimes; i++) 
            free(tbl[i]);
      }
   
      n = 1L << R.k;
   
      for (i = 0; i < NumPrimes; i++) {
         if ( !(tbl[i] = (long *) NTL_MALLOC(n, sizeof(long), 0)) )
            Error("out of space in FFTRep");
      }

      k = MaxK = R.k;
   }
   else {
      k = R.k;
   }

   long i, j, n;

   n = 1L << k;

   for (i = 0; i < NumPrimes; i++)
      for (j = 0; j < n; j++)
         tbl[i][j] = R.tbl[i][j];

   return *this;
}

FFTRep::~FFTRep()
{
   if (MaxK == -1)
      return;

   for (long i = 0; i < NumPrimes; i++)
      free(tbl[i]);

   free(tbl);
}



void ZZ_pXModRep::SetSize(long NewN)
{
   ZZ_pInfo->check();

   NumPrimes = ZZ_pInfo->NumPrimes;

   if (NewN < 0)
      Error("bad arg to ZZ_pXModRep::SetSize()");

   if (NewN <= MaxN) {
      n = NewN;
      return;
   }

   long i;
 

   if (MaxN == 0) {
      tbl = (long **) NTL_MALLOC(ZZ_pInfo->NumPrimes, sizeof(long *), 0);
      if (!tbl)
         Error("out of space in ZZ_pXModRep::SetSize()");
   }
   else {
      for (i = 0; i < ZZ_pInfo->NumPrimes; i++) 
         free(tbl[i]);
   }

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      if ( !(tbl[i] = (long *) NTL_MALLOC(NewN, sizeof(long), 0)) )
         Error("out of space in ZZ_pXModRep::SetSize()");
   }

   n = MaxN = NewN;
}

ZZ_pXModRep::~ZZ_pXModRep()
{
   if (MaxN == 0)
      return;

   long i;
   for (i = 0; i < NumPrimes; i++)
      free(tbl[i]);

   free(tbl);
}


static vec_long ModularRepBuf;


void ToModularRep(vec_long& x, const ZZ_p& a)
{
   ZZ_pInfo->check();
   ZZ_p_rem_struct_eval(ZZ_pInfo->rem_struct, &x[0], rep(a));
}


// NOTE: earlier versions used Kahan summation...
// we no longer do this, as it is less portable than I thought.



void FromModularRep(ZZ_p& x, const vec_long& a)
{
   ZZ_pInfo->check();

   long n = ZZ_pInfo->NumPrimes;
   static ZZ q, s, t;
   long i;
   double y;

   if (ZZ_p_crt_struct_special(ZZ_pInfo->crt_struct)) {
      ZZ_p_crt_struct_eval(ZZ_pInfo->crt_struct, t, &a[0]);
      x.LoopHole() = t;
      return;
   }
      

   if (ZZ_pInfo->QuickCRT) {
      y = 0;
      for (i = 0; i < n; i++)
         y += ((double) a[i])*ZZ_pInfo->x[i];

      conv(q, (y + 0.5)); 
   } else {
      long Q, r;
      static ZZ qq;

      y = 0;

      clear(q);

      for (i = 0; i < n; i++) {
         r = MulDivRem(Q, a[i], ZZ_pInfo->u[i], FFTPrime[i], ZZ_pInfo->x[i]);
         add(q, q, Q);
         y += r*FFTPrimeInv[i];
      }

      conv(qq, (y + 0.5));
      add(q, q, qq);
   }

   ZZ_p_crt_struct_eval(ZZ_pInfo->crt_struct, t, &a[0]);

   mul(s, q, ZZ_pInfo->MinusMModP);
   add(t, t, s);

   conv(x, t);
}




void ToFFTRep(FFTRep& y, const ZZ_pX& x, long k, long lo, long hi)
// computes an n = 2^k point convolution.
// if deg(x) >= 2^k, then x is first reduced modulo X^n-1.
{
   ZZ_pInfo->check();

   long n, i, j, m, j1;
   vec_long& t = ModularRepBuf;
   ZZ_p accum;


   if (k > ZZ_pInfo->MaxRoot) 
      Error("Polynomial too big for FFT");

   if (lo < 0)
      Error("bad arg to ToFFTRep");

   t.SetLength(ZZ_pInfo->NumPrimes);

   hi = min(hi, deg(x));

   y.SetSize(k);

   n = 1L << k;

   m = max(hi-lo + 1, 0);

   const ZZ_p *xx = x.rep.elts();

   for (j = 0; j < n; j++) {
      if (j >= m) {
         for (i = 0; i < ZZ_pInfo->NumPrimes; i++)
            y.tbl[i][j] = 0;
      }
      else {
         accum = xx[j+lo];
         for (j1 = j + n; j1 < m; j1 += n)
            add(accum, accum, xx[j1+lo]);
         ToModularRep(t, accum);
         for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
            y.tbl[i][j] = t[i];
         }
      }
   }


   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTFwd(yp, yp, k, i);
   }
}



void RevToFFTRep(FFTRep& y, const vec_ZZ_p& x, 
                 long k, long lo, long hi, long offset)
// computes an n = 2^k point convolution of X^offset*x[lo..hi] mod X^n-1
// using "inverted" evaluation points.

{
   ZZ_pInfo->check();

   long n, i, j, m, j1;
   vec_long& t = ModularRepBuf;
   ZZ_p accum;

   if (k > ZZ_pInfo->MaxRoot) 
      Error("Polynomial too big for FFT");

   if (lo < 0)
      Error("bad arg to ToFFTRep");

   t.SetLength(ZZ_pInfo->NumPrimes);

   hi = min(hi, x.length()-1);

   y.SetSize(k);

   n = 1L << k;

   m = max(hi-lo + 1, 0);

   const ZZ_p *xx = x.elts();

   offset = offset & (n-1);

   for (j = 0; j < n; j++) {
      if (j >= m) {
         for (i = 0; i < ZZ_pInfo->NumPrimes; i++)
            y.tbl[i][offset] = 0;
      }
      else {
         accum = xx[j+lo];
         for (j1 = j + n; j1 < m; j1 += n)
            add(accum, accum, xx[j1+lo]);
         ToModularRep(t, accum);
         for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
            y.tbl[i][offset] = t[i];

         }
      }

      offset = (offset + 1) & (n-1);
   }


   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTRev(yp, yp, k, i);
      FFTMulTwoInv(yp, yp, k, i);
   }

}

void FromFFTRep(ZZ_pX& x, FFTRep& y, long lo, long hi)

   // converts from FFT-representation to coefficient representation
   // only the coefficients lo..hi are computed
   

{
   ZZ_pInfo->check();

   long k, n, i, j, l;

   vec_long& t = ModularRepBuf;

   t.SetLength(ZZ_pInfo->NumPrimes);

   k = y.k;
   n = (1L << k);


   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTRev(yp, yp, k, i);
      FFTMulTwoInv(yp, yp, k, i);
   }

   hi = min(hi, n-1);
   l = hi-lo+1;
   l = max(l, 0);
   x.rep.SetLength(l);

   for (j = 0; j < l; j++) {
      for (i = 0; i < ZZ_pInfo->NumPrimes; i++) 
         t[i] = y.tbl[i][j+lo]; 

      FromModularRep(x.rep[j], t);
   }

   x.normalize();
}

void RevFromFFTRep(vec_ZZ_p& x, FFTRep& y, long lo, long hi)

   // converts from FFT-representation to coefficient representation
   // using "inverted" evaluation points.
   // only the coefficients lo..hi are computed
   

{
   ZZ_pInfo->check();

   long k, n, i, j, l;

   vec_long& t = ModularRepBuf;

   k = y.k;
   n = (1L << k);

   t.SetLength(ZZ_pInfo->NumPrimes);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTFwd(yp, yp, k, i);
   }

   hi = min(hi, n-1);
   l = hi-lo+1;
   l = max(l, 0);
   x.SetLength(l);

   for (j = 0; j < l; j++) {
      for (i = 0; i < ZZ_pInfo->NumPrimes; i++) 
         t[i] = y.tbl[i][j+lo]; 

      FromModularRep(x[j], t);
   }
}

void NDFromFFTRep(ZZ_pX& x, const FFTRep& y, long lo, long hi, FFTRep& z)
{
   ZZ_pInfo->check();

   long k, n, i, j, l;

   vec_long& t = ModularRepBuf;

   t.SetLength(ZZ_pInfo->NumPrimes);
   k = y.k;
   n = (1L << k);

   z.SetSize(k);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *zp = &z.tbl[i][0];
      const long *yp = &y.tbl[i][0];

      FFTRev(zp, yp, k, i);
      FFTMulTwoInv(zp, zp, k, i);
   }

   hi = min(hi, n-1);
   l = hi-lo+1;
   l = max(l, 0);
   x.rep.SetLength(l);

   for (j = 0; j < l; j++) {
      for (i = 0; i < ZZ_pInfo->NumPrimes; i++) 
         t[i] = z.tbl[i][j+lo]; 

      FromModularRep(x.rep[j], t);
   }

   x.normalize();
}

void NDFromFFTRep(ZZ_pX& x, FFTRep& y, long lo, long hi)
{
   FFTRep z;
   NDFromFFTRep(x, y, lo, hi, z);
}

void FromFFTRep(ZZ_p* x, FFTRep& y, long lo, long hi)

   // converts from FFT-representation to coefficient representation
   // only the coefficients lo..hi are computed
   

{
   ZZ_pInfo->check();

   long k, n, i, j;

   vec_long& t = ModularRepBuf;

   k = y.k;
   n = (1L << k);

   t.SetLength(ZZ_pInfo->NumPrimes);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *yp = &y.tbl[i][0];
      FFTRev(yp, yp, k, i);
      FFTMulTwoInv(yp, yp, k, i);
   }

   for (j = lo; j <= hi; j++) {
      if (j >= n)
         clear(x[j-lo]);
      else {
         for (i = 0; i < ZZ_pInfo->NumPrimes; i++) 
            t[i] = y.tbl[i][j]; 

         FromModularRep(x[j-lo], t);
      }
   }
}


void mul(FFTRep& z, const FFTRep& x, const FFTRep& y)
{
   ZZ_pInfo->check();

   long k, n, i, j;

   if (x.k != y.k) Error("FFT rep mismatch");

   k = x.k;
   n = 1L << k;

   z.SetSize(k);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *zp = &z.tbl[i][0];
      const long *xp = &x.tbl[i][0];
      const long *yp = &y.tbl[i][0];
      long q = FFTPrime[i];
      double qinv = ((double) 1)/((double) q);

      for (j = 0; j < n; j++)
         zp[j] = MulMod(xp[j], yp[j], q, qinv);
   }

}

void sub(FFTRep& z, const FFTRep& x, const FFTRep& y)
{
   ZZ_pInfo->check();

   long k, n, i, j;

   if (x.k != y.k) Error("FFT rep mismatch");

   k = x.k;
   n = 1L << k;

   z.SetSize(k);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *zp = &z.tbl[i][0];
      const long *xp = &x.tbl[i][0];
      const long *yp = &y.tbl[i][0];
      long q = FFTPrime[i];

      for (j = 0; j < n; j++)
         zp[j] = SubMod(xp[j], yp[j], q);
   }
}

void add(FFTRep& z, const FFTRep& x, const FFTRep& y)
{
   ZZ_pInfo->check();

   long k, n, i, j;

   if (x.k != y.k) Error("FFT rep mismatch");

   k = x.k;
   n = 1L << k;

   z.SetSize(k);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long *zp = &z.tbl[i][0];
      const long *xp = &x.tbl[i][0];
      const long *yp = &y.tbl[i][0];
      long q = FFTPrime[i];

      for (j = 0; j < n; j++)
         zp[j] = AddMod(xp[j], yp[j], q);
   }
}


void reduce(FFTRep& x, const FFTRep& a, long k)
  // reduces a 2^l point FFT-rep to a 2^k point FFT-rep
  // input may alias output
{
   ZZ_pInfo->check();

   long i, j, l, n;
   long* xp;
   const long* ap;

   l = a.k;
   n = 1L << k;

   if (l < k) Error("reduce: bad operands");

   x.SetSize(k);

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      ap = &a.tbl[i][0];   
      xp = &x.tbl[i][0];
      for (j = 0; j < n; j++) 
         xp[j] = ap[j << (l-k)];
   }
}

void AddExpand(FFTRep& x, const FFTRep& a)
//  x = x + (an "expanded" version of a)
{
   ZZ_pInfo->check();

   long i, j, l, k, n;

   l = x.k;
   k = a.k;
   n = 1L << k;

   if (l < k) Error("AddExpand: bad args");

   for (i = 0; i < ZZ_pInfo->NumPrimes; i++) {
      long q = FFTPrime[i];
      const long *ap = &a.tbl[i][0];
      long *xp = &x.tbl[i][0];
      for (j = 0; j < n; j++) {
         long j1 = j << (l-k);
         xp[j1] = AddMod(xp[j1], ap[j], q);
      }
   }
}



void ToZZ_pXModRep(ZZ_pXModRep& y, const ZZ_pX& x, long lo, long hi)
{
   ZZ_pInfo->check();

   long n, i, j;
   vec_long& t = ModularRepBuf;

   t.SetLength(ZZ_pInfo->NumPrimes);

   if (lo < 0)
      Error("bad arg to ToZZ_pXModRep");
   hi = min(hi, deg(x));
   n = max(hi-lo+1, 0);

   y.SetSize(n);

   const ZZ_p *xx = x.rep.elts();

   for (j = 0; j < n; j++) {
      ToModularRep(t, xx[j+lo]);
      for (i = 0; i < ZZ_pInfo->NumPrimes; i++)
         y.tbl[i][j] = t[i];
   }
}


void ToFFTRep(FFTRep& x, const ZZ_pXModRep& a, long k, long lo, long hi)
{
   ZZ_pInfo->check();

   vec_long s;
   long n, m, i, j;

   if (k < 0 || lo < 0)
      Error("bad args to ToFFTRep");

   if (hi > a.n-1) hi = a.n-1;

   n = 1L << k;
   m = max(hi-lo+1, 0);

   if (m > n)
      Error("bad args to ToFFTRep");

   s.SetLength(n);
   long *sp = s.elts();

   x.SetSize(k);

   long NumPrimes = ZZ_pInfo->NumPrimes;

   for (i = 0; i < NumPrimes; i++) {
      long *xp = &x.tbl[i][0];
      long *ap = (m == 0 ? 0 : &a.tbl[i][0]);
      for (j = 0; j < m; j++)
         sp[j] = ap[lo+j];
      for (j = m; j < n; j++)
         sp[j] = 0;
      
      FFTFwd(xp, sp, k, i);
   }
}






void FFTMul(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b)
{
   long k, d;

   if (IsZero(a) || IsZero(b)) {
      clear(x);
      return;
   }

   d = deg(a) + deg(b);
   k = NextPowerOfTwo(d+1);

   FFTRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);

   ToFFTRep(R1, a, k);
   ToFFTRep(R2, b, k);
   mul(R1, R1, R2);
   FromFFTRep(x, R1, 0, d);
}

void FFTSqr(ZZ_pX& x, const ZZ_pX& a)
{
   long k, d;

   if (IsZero(a)) {
      clear(x);
      return;
   }

   d = 2*deg(a);
   k = NextPowerOfTwo(d+1);

   FFTRep R1(INIT_SIZE, k);

   ToFFTRep(R1, a, k);
   mul(R1, R1, R1);
   FromFFTRep(x, R1, 0, d);
}


void CopyReverse(ZZ_pX& x, const ZZ_pX& a, long lo, long hi)

   // x[0..hi-lo] = reverse(a[lo..hi]), with zero fill
   // input may not alias output

{
   long i, j, n, m;

   n = hi-lo+1;
   m = a.rep.length();

   x.rep.SetLength(n);

   const ZZ_p* ap = a.rep.elts();
   ZZ_p* xp = x.rep.elts();

   for (i = 0; i < n; i++) {
      j = hi-i;
      if (j < 0 || j >= m)
         clear(xp[i]);
      else
         xp[i] = ap[j];
   }

   x.normalize();
} 

void copy(ZZ_pX& x, const ZZ_pX& a, long lo, long hi)

   // x[0..hi-lo] = a[lo..hi], with zero fill
   // input may not alias output

{
   long i, j, n, m;

   n = hi-lo+1;
   m = a.rep.length();

   x.rep.SetLength(n);

   const ZZ_p* ap = a.rep.elts();
   ZZ_p* xp = x.rep.elts();

   for (i = 0; i < n; i++) {
      j = lo + i;
      if (j < 0 || j >= m)
         clear(xp[i]);
      else
         xp[i] = ap[j];
   }

   x.normalize();
} 


void rem21(ZZ_pX& x, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long i, da, ds, n, kk;

   da = deg(a);
   n = F.n;

   if (da > 2*n-2)
      Error("bad args to rem(ZZ_pX,ZZ_pX,ZZ_pXModulus)");


   if (da < n) {
      x = a;
      return;
   }

   if (!F.UseFFT || da - n <= NTL_ZZ_pX_FFT_CROSSOVER) {
      PlainRem(x, a, F.f);
      return;
   }

   FFTRep R1(INIT_SIZE, F.l);
   ZZ_pX P1(INIT_SIZE, n);

   ToFFTRep(R1, a, F.l, n, 2*(n-1));
   mul(R1, R1, F.HRep);
   FromFFTRep(P1, R1, n-2, 2*n-4);

   ToFFTRep(R1, P1, F.k);
   mul(R1, R1, F.FRep);
   FromFFTRep(P1, R1, 0, n-1);

   ds = deg(P1);

   kk = 1L << F.k;

   x.rep.SetLength(n);
   const ZZ_p* aa = a.rep.elts();
   const ZZ_p* ss = P1.rep.elts();
   ZZ_p* xx = x.rep.elts();

   for (i = 0; i < n; i++) {
      if (i <= ds)
         sub(xx[i], aa[i], ss[i]);
      else
         xx[i] = aa[i];

      if (i + kk <= da)
         add(xx[i], xx[i], aa[i+kk]);
   }

   x.normalize();
}

void DivRem21(ZZ_pX& q, ZZ_pX& x, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long i, da, ds, n, kk;

   da = deg(a);
   n = F.n;

   if (da > 2*n-2)
      Error("bad args to rem(ZZ_pX,ZZ_pX,ZZ_pXModulus)");


   if (da < n) {
      x = a;
      clear(q);
      return;
   }

   if (!F.UseFFT || da - n <= NTL_ZZ_pX_FFT_CROSSOVER) {
      PlainDivRem(q, x, a, F.f);
      return;
   }

   FFTRep R1(INIT_SIZE, F.l);
   ZZ_pX P1(INIT_SIZE, n), qq;

   ToFFTRep(R1, a, F.l, n, 2*(n-1));
   mul(R1, R1, F.HRep);
   FromFFTRep(P1, R1, n-2, 2*n-4);
   qq = P1;

   ToFFTRep(R1, P1, F.k);
   mul(R1, R1, F.FRep);
   FromFFTRep(P1, R1, 0, n-1);

   ds = deg(P1);

   kk = 1L << F.k;

   x.rep.SetLength(n);
   const ZZ_p* aa = a.rep.elts();
   const ZZ_p* ss = P1.rep.elts();
   ZZ_p* xx = x.rep.elts();

   for (i = 0; i < n; i++) {
      if (i <= ds)
         sub(xx[i], aa[i], ss[i]);
      else
         xx[i] = aa[i];

      if (i + kk <= da)
         add(xx[i], xx[i], aa[i+kk]);
   }

   x.normalize();
   q = qq;
}

void div21(ZZ_pX& x, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long da, n;

   da = deg(a);
   n = F.n;

   if (da > 2*n-2)
      Error("bad args to rem(ZZ_pX,ZZ_pX,ZZ_pXModulus)");


   if (da < n) {
      clear(x);
      return;
   }

   if (!F.UseFFT || da - n <= NTL_ZZ_pX_FFT_CROSSOVER) {
      PlainDiv(x, a, F.f);
      return;
   }

   FFTRep R1(INIT_SIZE, F.l);
   ZZ_pX P1(INIT_SIZE, n);

   ToFFTRep(R1, a, F.l, n, 2*(n-1));
   mul(R1, R1, F.HRep);
   FromFFTRep(x, R1, n-2, 2*n-4);
}


void rem(ZZ_pX& x, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long da = deg(a);
   long n = F.n;

   if (n < 0) Error("rem: unitialized modulus");

   if (da <= 2*n-2) {
      rem21(x, a, F);
      return;
   }
   else if (!F.UseFFT || da - n <= NTL_ZZ_pX_FFT_CROSSOVER) {
      PlainRem(x, a, F.f);
      return;
   }

   ZZ_pX buf(INIT_SIZE, 2*n-1);

   long a_len = da+1;

   while (a_len > 0) {
      long old_buf_len = buf.rep.length();
      long amt = min(2*n-1-old_buf_len, a_len);

      buf.rep.SetLength(old_buf_len+amt);

      long i;

      for (i = old_buf_len+amt-1; i >= amt; i--)
         buf.rep[i] = buf.rep[i-amt];

      for (i = amt-1; i >= 0; i--)
         buf.rep[i] = a.rep[a_len-amt+i];

      buf.normalize();

      rem21(buf, buf, F);

      a_len -= amt;
   }

   x = buf;
}

void DivRem(ZZ_pX& q, ZZ_pX& r, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long da = deg(a);
   long n = F.n;

   if (n < 0) Error("uninitialized modulus");

   if (da <= 2*n-2) {
      DivRem21(q, r, a, F);
      return;
   }
   else if (!F.UseFFT || da - n <= NTL_ZZ_pX_FFT_CROSSOVER) {
      PlainDivRem(q, r, a, F.f);
      return;
   }

   ZZ_pX buf(INIT_SIZE, 2*n-1);
   ZZ_pX qbuf(INIT_SIZE, n-1);

   ZZ_pX qq;
   qq.rep.SetLength(da-n+1);

   long a_len = da+1;
   long q_hi = da-n+1;

   while (a_len > 0) {
      long old_buf_len = buf.rep.length();
      long amt = min(2*n-1-old_buf_len, a_len);

      buf.rep.SetLength(old_buf_len+amt);

      long i;

      for (i = old_buf_len+amt-1; i >= amt; i--)
         buf.rep[i] = buf.rep[i-amt];

      for (i = amt-1; i >= 0; i--)
         buf.rep[i] = a.rep[a_len-amt+i];

      buf.normalize();

      DivRem21(qbuf, buf, buf, F);
      long dl = qbuf.rep.length();
      a_len = a_len - amt;
      for(i = 0; i < dl; i++)
         qq.rep[a_len+i] = qbuf.rep[i];
      for(i = dl+a_len; i < q_hi; i++)
         clear(qq.rep[i]);
      q_hi = a_len;
   }

   r = buf;

   qq.normalize();
   q = qq;
}

void div(ZZ_pX& q, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long da = deg(a);
   long n = F.n;

   if (n < 0) Error("uninitialized modulus");

   if (da <= 2*n-2) {
      div21(q, a, F);
      return;
   }
   else if (!F.UseFFT || da - n <= NTL_ZZ_pX_FFT_CROSSOVER) {
      PlainDiv(q, a, F.f);
      return;
   }

   ZZ_pX buf(INIT_SIZE, 2*n-1);
   ZZ_pX qbuf(INIT_SIZE, n-1);

   ZZ_pX qq;
   qq.rep.SetLength(da-n+1);

   long a_len = da+1;
   long q_hi = da-n+1;

   while (a_len > 0) {
      long old_buf_len = buf.rep.length();
      long amt = min(2*n-1-old_buf_len, a_len);

      buf.rep.SetLength(old_buf_len+amt);

      long i;

      for (i = old_buf_len+amt-1; i >= amt; i--)
         buf.rep[i] = buf.rep[i-amt];

      for (i = amt-1; i >= 0; i--)
         buf.rep[i] = a.rep[a_len-amt+i];

      buf.normalize();

      a_len = a_len - amt;
      if (a_len > 0)
         DivRem21(qbuf, buf, buf, F);
      else
         div21(qbuf, buf, F);

      long dl = qbuf.rep.length();
      for(i = 0; i < dl; i++)
         qq.rep[a_len+i] = qbuf.rep[i];
      for(i = dl+a_len; i < q_hi; i++)
         clear(qq.rep[i]);
      q_hi = a_len;
   }

   qq.normalize();
   q = qq;
}



void MulMod(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b, const ZZ_pXModulus& F)
{
   long  da, db, d, n, k;

   da = deg(a);
   db = deg(b);
   n = F.n;

   if (n < 0) Error("MulMod: uninitialized modulus");

   if (da >= n || db >= n)
      Error("bad args to MulMod(ZZ_pX,ZZ_pX,ZZ_pX,ZZ_pXModulus)");

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   if (!F.UseFFT || da <= NTL_ZZ_pX_FFT_CROSSOVER || db <= NTL_ZZ_pX_FFT_CROSSOVER) {
      ZZ_pX P1;
      mul(P1, a, b);
      rem(x, P1, F);
      return;
   }

   d = da + db + 1;

   k = NextPowerOfTwo(d);
   k = max(k, F.k);

   FFTRep R1(INIT_SIZE, k), R2(INIT_SIZE, F.l);
   ZZ_pX P1(INIT_SIZE, n);

   ToFFTRep(R1, a, k);
   ToFFTRep(R2, b, k);
   mul(R1, R1, R2);
   NDFromFFTRep(P1, R1, n, d-1, R2); // save R1 for future use
   
   ToFFTRep(R2, P1, F.l);
   mul(R2, R2, F.HRep);
   FromFFTRep(P1, R2, n-2, 2*n-4);

   ToFFTRep(R2, P1, F.k);
   mul(R2, R2, F.FRep);
   reduce(R1, R1, F.k);
   sub(R1, R1, R2);
   FromFFTRep(x, R1, 0, n-1);
}

void SqrMod(ZZ_pX& x, const ZZ_pX& a, const ZZ_pXModulus& F)
{
   long  da, d, n, k;

   da = deg(a);
   n = F.n;

   if (n < 0) Error("SqrMod: uninitailized modulus");

   if (da >= n) 
      Error("bad args to SqrMod(ZZ_pX,ZZ_pX,ZZ_pXModulus)");

   if (!F.UseFFT || da <= NTL_ZZ_pX_FFT_CROSSOVER) {
      ZZ_pX P1;
      sqr(P1, a);
      rem(x, P1, F);
      return;
   }


   d = 2*da + 1;

   k = NextPowerOfTwo(d);
   k = max(k, F.k);

   FFTRep R1(INIT_SIZE, k), R2(INIT_SIZE, F.l);
   ZZ_pX P1(INIT_SIZE, n);

   ToFFTRep(R1, a, k);
   mul(R1, R1, R1);
   NDFromFFTRep(P1, R1, n, d-1, R2);  // save R1 for future use
   
   ToFFTRep(R2, P1, F.l);
   mul(R2, R2, F.HRep);
   FromFFTRep(P1, R2, n-2, 2*n-4);

   ToFFTRep(R2, P1, F.k);
   mul(R2, R2, F.FRep);
   reduce(R1, R1, F.k);
   sub(R1, R1, R2);
   FromFFTRep(x, R1, 0, n-1);
}

void PlainInvTrunc(ZZ_pX& x, const ZZ_pX& a, long m)

   /* x = (1/a) % X^m, input not output, constant term a is nonzero */

{
   long i, k, n, lb;
   static ZZ v, t;
   ZZ_p s;
   const ZZ_p* ap;
   ZZ_p* xp;
   

   n = deg(a);

   if (n < 0) Error("division by zero");

   inv(s, ConstTerm(a));

   if (n == 0) {
      conv(x, s);
      return;
   }

   ap = a.rep.elts();
   x.rep.SetLength(m);
   xp = x.rep.elts();

   xp[0] = s;

   long is_one = IsOne(s);

   for (k = 1; k < m; k++) {
      clear(v);
      lb = max(k-n, 0);
      for (i = lb; i <= k-1; i++) {
         mul(t, rep(xp[i]), rep(ap[k-i]));
         add(v, v, t);
      }
      conv(xp[k], v);
      negate(xp[k], xp[k]);
      if (!is_one) mul(xp[k], xp[k], s);
   }

   x.normalize();
}


void trunc(ZZ_pX& x, const ZZ_pX& a, long m)

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
      ZZ_p* xp;
      const ZZ_p* ap;

      n = min(a.rep.length(), m);
      x.rep.SetLength(n);

      xp = x.rep.elts();
      ap = a.rep.elts();

      for (i = 0; i < n; i++) xp[i] = ap[i];

      x.normalize();
   }
}

void CyclicReduce(ZZ_pX& x, const ZZ_pX& a, long m)

// computes x = a mod X^m-1

{
   long n = deg(a);
   long i, j;
   ZZ_p accum;

   if (n < m) {
      x = a;
      return;
   }

   if (&x != &a)
      x.rep.SetLength(m);

   for (i = 0; i < m; i++) {
      accum = a.rep[i];
      for (j = i + m; j <= n; j += m)
         add(accum, accum, a.rep[j]);
      x.rep[i] = accum;
   }

   if (&x == &a)
      x.rep.SetLength(m);

   x.normalize();
}



void InvTrunc(ZZ_pX& x, const ZZ_pX& a, long m)
{
   if (m < 0) Error("InvTrunc: bad args");

   if (m == 0) {
      clear(x);
      return;
   }

   if (NTL_OVERFLOW(m, 1, 0))
      Error("overflow in InvTrunc");

   if (&x == &a) {
      ZZ_pX la;
      la = a;
      if (m > NTL_ZZ_pX_NEWTON_CROSSOVER && deg(a) > 0)
         NewtonInvTrunc(x, la, m);
      else
         PlainInvTrunc(x, la, m);
   }
   else {
      if (m > NTL_ZZ_pX_NEWTON_CROSSOVER && deg(a) > 0)
         NewtonInvTrunc(x, a, m);
      else
         PlainInvTrunc(x, a, m);
   }
}
   


void build(ZZ_pXModulus& x, const ZZ_pX& f)
{
   x.f = f;
   x.n = deg(f);

   x.tracevec.SetLength(0);

   if (x.n <= 0)
      Error("build: deg(f) must be at least 1");

   if (x.n <= NTL_ZZ_pX_FFT_CROSSOVER + 1) {
      x.UseFFT = 0;
      return;
   }

   x.UseFFT = 1;

   x.k = NextPowerOfTwo(x.n);
   x.l = NextPowerOfTwo(2*x.n - 3);
   ToFFTRep(x.FRep, f, x.k);

   ZZ_pX P1(INIT_SIZE, x.n+1), P2(INIT_SIZE, x.n);

   CopyReverse(P1, f, 0, x.n);
   InvTrunc(P2, P1, x.n-1);

   CopyReverse(P1, P2, 0, x.n-2);
   ToFFTRep(x.HRep, P1, x.l);
}

ZZ_pXModulus::ZZ_pXModulus(const ZZ_pX& ff)
{
   build(*this, ff);
}

ZZ_pXMultiplier::ZZ_pXMultiplier(const ZZ_pX& b, const ZZ_pXModulus& F)
{
   build(*this, b, F); 
}

void build(ZZ_pXMultiplier& x, const ZZ_pX& b, 
                         const ZZ_pXModulus& F)
{
   long db;
   long n = F.n;

   if (n < 0) Error("build ZZ_pXMultiplier: uninitialized modulus");

   x.b = b;
   db = deg(b);

   if (db >= n) Error("build ZZ_pXMultiplier: deg(b) >= deg(f)");

   if (!F.UseFFT || db <= NTL_ZZ_pX_FFT_CROSSOVER) {
      x.UseFFT = 0;
      return;
   }

   x.UseFFT = 1;

   FFTRep R1(INIT_SIZE, F.l);
   ZZ_pX P1(INIT_SIZE, n);
   

   ToFFTRep(R1, b, F.l);
   reduce(x.B2, R1, F.k);
   mul(R1, R1, F.HRep);
   FromFFTRep(P1, R1, n-1, 2*n-3); 
   ToFFTRep(x.B1, P1, F.l);
}


void MulMod(ZZ_pX& x, const ZZ_pX& a, const ZZ_pXMultiplier& B,
                                      const ZZ_pXModulus& F)
{

   long n = F.n;
   long da;

   da = deg(a);

   if (da >= n)
      Error(" bad args to MulMod(ZZ_pX,ZZ_pX,ZZ_pXMultiplier,ZZ_pXModulus)");

   if (da < 0) {
      clear(x);
      return;
   }

   if (!B.UseFFT || !F.UseFFT || da <= NTL_ZZ_pX_FFT_CROSSOVER) {
      ZZ_pX P1;
      mul(P1, a, B.b);
      rem(x, P1, F);
      return;
   }

   ZZ_pX P1(INIT_SIZE, n), P2(INIT_SIZE, n);
   FFTRep R1(INIT_SIZE, F.l), R2(INIT_SIZE, F.l);

   ToFFTRep(R1, a, F.l);
   mul(R2, R1, B.B1);
   FromFFTRep(P1, R2, n-1, 2*n-3);

   reduce(R1, R1, F.k);
   mul(R1, R1, B.B2);
   ToFFTRep(R2, P1, F.k);
   mul(R2, R2, F.FRep);
   sub(R1, R1, R2);

   FromFFTRep(x, R1, 0, n-1);
}
   

void PowerXMod(ZZ_pX& hh, const ZZ& e, const ZZ_pXModulus& F)
{
   if (F.n < 0) Error("PowerXMod: uninitialized modulus");

   if (IsZero(e)) {
      set(hh);
      return;
   }

   long n = NumBits(e);
   long i;

   ZZ_pX h;

   h.SetMaxLength(F.n);
   set(h);

   for (i = n - 1; i >= 0; i--) {
      SqrMod(h, h, F);
      if (bit(e, i))
         MulByXMod(h, h, F);
   }

   if (e < 0) InvMod(h, h, F);

   hh = h;
}


void PowerXPlusAMod(ZZ_pX& hh, const ZZ_p& a, const ZZ& e, const ZZ_pXModulus& F)
{
   if (F.n < 0) Error("PowerXPlusAMod: uninitialized modulus");

   if (IsZero(e)) {
      set(hh);
      return;
   }

   ZZ_pX t1(INIT_SIZE, F.n), t2(INIT_SIZE, F.n);
   long n = NumBits(e);
   long i;

   ZZ_pX h;

   h.SetMaxLength(F.n);
   set(h);

   for (i = n - 1; i >= 0; i--) {
      SqrMod(h, h, F);
      if (bit(e, i)) {
         MulByXMod(t1, h, F);
         mul(t2, h, a);
         add(h, t1, t2);
      }
   }

   if (e < 0) InvMod(h, h, F);

   hh = h;
}


void PowerMod(ZZ_pX& h, const ZZ_pX& g, const ZZ& e, const ZZ_pXModulus& F)
{
   if (deg(g) >= F.n)
      Error("PowerMod: bad args");

   if (IsZero(e)) {
      set(h);
      return;
   }

   ZZ_pXMultiplier G;

   ZZ_pX res;

   long n = NumBits(e);
   long i;

   build(G, g, F);

   res.SetMaxLength(F.n);
   set(res);

   for (i = n - 1; i >= 0; i--) {
      SqrMod(res, res, F);
      if (bit(e, i))
         MulMod(res, res, G, F);
   }

   if (e < 0) InvMod(res, res, F);

   h = res;
}


void NewtonInvTrunc(ZZ_pX& x, const ZZ_pX& a, long m)
{
   x.SetMaxLength(m);
   long i, t, k;

   long log2_newton = NextPowerOfTwo(NTL_ZZ_pX_NEWTON_CROSSOVER)-1;
   PlainInvTrunc(x, a, 1L << log2_newton);

   t = NextPowerOfTwo(m);

   FFTRep R1(INIT_SIZE, t), R2(INIT_SIZE, t);
   ZZ_pX P1(INIT_SIZE, m/2);

   long a_len = min(m, a.rep.length());

   ZZ_pXModRep a_rep;
   ToZZ_pXModRep(a_rep, a, 0, a_len-1);

   k = 1L << log2_newton; 
   t = log2_newton;

   while (k < m) {
      long l = min(2*k, m);

      ToFFTRep(R1, x, t+1);
      ToFFTRep(R2, a_rep, t+1, 0, l-1); 
      mul(R2, R2, R1);
      FromFFTRep(P1, R2, k, l-1);
      
      ToFFTRep(R2, P1, t+1);
      mul(R2, R2, R1);
      FromFFTRep(P1, R2, 0, l-k-1);

      x.rep.SetLength(l);
      long y_len = P1.rep.length();
      for (i = k; i < l; i++) {
         if (i-k >= y_len)
            clear(x.rep[i]);
         else
            negate(x.rep[i], P1.rep[i-k]);
      }
      x.normalize();

      t++;
      k = l;
   }
}



void FFTDivRem(ZZ_pX& q, ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b)
{
   long n = deg(b);
   long m = deg(a);
   long k, l;

   if (m < n) {
      clear(q);
      r = a;
      return;
   }

   if (m >= 3*n) {
      ZZ_pXModulus B;
      build(B, b);
      DivRem(q, r, a, B);
      return;
   }

   ZZ_pX P1, P2, P3;

   CopyReverse(P3, b, 0, n);
   InvTrunc(P2, P3, m-n+1);
   CopyReverse(P1, P2, 0, m-n);

   k = NextPowerOfTwo(2*(m-n)+1);
   long k1 = NextPowerOfTwo(n);
   long mx = max(k1, k);

   FFTRep R1(INIT_SIZE, mx), R2(INIT_SIZE, mx);

   ToFFTRep(R1, P1, k);
   ToFFTRep(R2, a, k, n, m);
   mul(R1, R1, R2);
   FromFFTRep(P3, R1, m-n, 2*(m-n));
   
   l = 1L << k1;

   
   ToFFTRep(R1, b, k1);
   ToFFTRep(R2, P3, k1);
   mul(R1, R1, R2);
   FromFFTRep(P1, R1, 0, n-1);
   CyclicReduce(P2, a, l);
   trunc(r, P2, n);
   sub(r, r, P1);
   q = P3;
}




void FFTDiv(ZZ_pX& q, const ZZ_pX& a, const ZZ_pX& b)
{

   long n = deg(b);
   long m = deg(a);
   long k;

   if (m < n) {
      clear(q);
      return;
   }

   if (m >= 3*n) {
      ZZ_pXModulus B;
      build(B, b);
      div(q, a, B);
      return;
   }

   ZZ_pX P1, P2, P3;

   CopyReverse(P3, b, 0, n);
   InvTrunc(P2, P3, m-n+1);
   CopyReverse(P1, P2, 0, m-n);

   k = NextPowerOfTwo(2*(m-n)+1);

   FFTRep R1(INIT_SIZE, k), R2(INIT_SIZE, k);

   ToFFTRep(R1, P1, k);
   ToFFTRep(R2, a, k, n, m);
   mul(R1, R1, R2);
   FromFFTRep(q, R1, m-n, 2*(m-n));
}



void FFTRem(ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b)
{
   long n = deg(b);
   long m = deg(a);
   long k, l;

   if (m < n) {
      r = a;
      return;
   }

   if (m >= 3*n) {
      ZZ_pXModulus B;
      build(B, b);
      rem(r, a, B);
      return;
   }

   ZZ_pX P1, P2, P3;

   CopyReverse(P3, b, 0, n);
   InvTrunc(P2, P3, m-n+1);
   CopyReverse(P1, P2, 0, m-n);

   k = NextPowerOfTwo(2*(m-n)+1);
   long k1 = NextPowerOfTwo(n);
   long mx = max(k, k1);

   FFTRep R1(INIT_SIZE, mx), R2(INIT_SIZE, mx);

   ToFFTRep(R1, P1, k);
   ToFFTRep(R2, a, k, n, m);
   mul(R1, R1, R2);
   FromFFTRep(P3, R1, m-n, 2*(m-n));
   
   l = 1L << k1;
   
   ToFFTRep(R1, b, k1);
   ToFFTRep(R2, P3, k1);
   mul(R1, R1, R2);
   FromFFTRep(P3, R1, 0, n-1);
   CyclicReduce(P2, a, l);
   trunc(r, P2, n);
   sub(r, r, P3);
}


void DivRem(ZZ_pX& q, ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b)
{
   if (deg(b) > NTL_ZZ_pX_DIV_CROSSOVER && deg(a) - deg(b) > NTL_ZZ_pX_DIV_CROSSOVER)
      FFTDivRem(q, r, a, b);
   else
      PlainDivRem(q, r, a, b);
}

void div(ZZ_pX& q, const ZZ_pX& a, const ZZ_pX& b)
{
   if (deg(b) > NTL_ZZ_pX_DIV_CROSSOVER && deg(a) - deg(b) > NTL_ZZ_pX_DIV_CROSSOVER)
      FFTDiv(q, a, b);
   else
      PlainDiv(q, a, b);
}

void div(ZZ_pX& q, const ZZ_pX& a, const ZZ_p& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();

   inv(T, b);
   mul(q, a, T);
}

void div(ZZ_pX& q, const ZZ_pX& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();

   T = b;
   inv(T, T);
   mul(q, a, T);
}



void rem(ZZ_pX& r, const ZZ_pX& a, const ZZ_pX& b)
{
   if (deg(b) > NTL_ZZ_pX_DIV_CROSSOVER && deg(a) - deg(b) > NTL_ZZ_pX_DIV_CROSSOVER)
      FFTRem(r, a, b);
   else
      PlainRem(r, a, b);
}


long operator==(const ZZ_pX& a, long b)
{
   if (b == 0)
      return IsZero(a);

   if (b == 1)
      return IsOne(a);

   long da = deg(a);

   if (da > 0)
      return 0;

   ZZ_pTemp TT; ZZ_p& bb = TT.val();
   bb = b;

   if (da < 0)
      return IsZero(bb);

   return a.rep[0] == bb;
}

long operator==(const ZZ_pX& a, const ZZ_p& b)
{
   if (IsZero(b))
      return IsZero(a);

   long da = deg(a);

   if (da != 0)
      return 0;

   return a.rep[0] == b;
}

void power(ZZ_pX& x, const ZZ_pX& a, long e)
{
   if (e < 0) {
      Error("power: negative exponent");
   }

   if (e == 0) {
      x = 1;
      return;
   }

   if (a == 0 || a == 1) {
      x = a;
      return;
   }

   long da = deg(a);

   if (da == 0) {
      x = power(ConstTerm(a), e);
      return;
   }

   if (da > (NTL_MAX_LONG-1)/e)
      Error("overflow in power");

   ZZ_pX res;
   res.SetMaxLength(da*e + 1);
   res = 1;
   
   long k = NumBits(e);
   long i;

   for (i = k - 1; i >= 0; i--) {
      sqr(res, res);
      if (bit(e, i))
         mul(res, res, a);
   }

   x = res;
}

void reverse(ZZ_pX& x, const ZZ_pX& a, long hi)
{
   if (hi < 0) { clear(x); return; }
   if (NTL_OVERFLOW(hi, 1, 0))
      Error("overflow in reverse");

   if (&x == &a) {
      ZZ_pX tmp;
      CopyReverse(tmp, a, 0, hi);
      x = tmp;
   }
   else
      CopyReverse(x, a, 0, hi);
}

NTL_END_IMPL
