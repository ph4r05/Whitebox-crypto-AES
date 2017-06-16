

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>


#include <NTL/new.h>



NTL_START_IMPL




const ZZ& ZZ::zero()
{
   static ZZ z;
   return z;
}


const ZZ& ZZ_expo(long e)
{
   static ZZ expo_helper;
   conv(expo_helper, e);
   return expo_helper;
}




void AddMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   static ZZ B;
   conv(B, b);
   AddMod(x, a, B, n);
}


void SubMod(ZZ& x, const ZZ& a, long b, const ZZ& n)
{
   static ZZ B;
   conv(B, b);
   SubMod(x, a, B, n);
}

void SubMod(ZZ& x, long a, const ZZ& b, const ZZ& n)
{
   static ZZ A;
   conv(A, a);
   SubMod(x, A, b, n);
}



// ****** input and output

static long iodigits = 0;
static long ioradix = 0;

// iodigits is the greatest integer such that 10^{iodigits} < NTL_WSP_BOUND
// ioradix = 10^{iodigits}

static void InitZZIO()
{
   long x;

   x = (NTL_WSP_BOUND-1)/10;
   iodigits = 0;
   ioradix = 1;

   while (x) {
      x = x / 10;
      iodigits++;
      ioradix = ioradix * 10;
   }

   if (iodigits <= 0) Error("problem with I/O");
}

istream& operator>>(istream& s, ZZ& x)
{
   long c;
   long cval;
   long sign;
   long ndigits;
   long acc;
   static ZZ a;

   if (!s) Error("bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   SkipWhiteSpace(s);
   c = s.peek();

   if (c == '-') {
      sign = -1;
      s.get();
      c = s.peek();
   }
   else
      sign = 1;

   cval = CharToIntVal(c);

   if (cval < 0 || cval > 9) Error("bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (cval >= 0 && cval <= 9) {
      acc = acc*10 + cval;
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      s.get();
      c = s.peek();
      cval = CharToIntVal(c);
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;

   return s;
}


// The class _ZZ_local_stack should be defined in an empty namespace,
// but since I don't want to rely on namespaces, we just give it a funny 
// name to avoid accidental name clashes.

struct _ZZ_local_stack {
   long top;
   long alloc;
   long *elts;

   _ZZ_local_stack() { top = -1; alloc = 0; elts = 0; }
   ~_ZZ_local_stack() { }

   long pop() { return elts[top--]; }
   long empty() { return (top == -1); }
   void push(long x);
};

void _ZZ_local_stack::push(long x)
{
   if (alloc == 0) {
      alloc = 100;
      elts = (long *) NTL_MALLOC(alloc, sizeof(long), 0);
   }

   top++;

   if (top + 1 > alloc) {
      alloc = 2*alloc;
      elts = (long *) NTL_REALLOC(elts, alloc, sizeof(long), 0);
   }

   if (!elts) {
      Error("out of space in ZZ output");
   }

   elts[top] = x;
}


static
void PrintDigits(ostream& s, long d, long justify)
{
   static char *buf = 0;

   if (!buf) {
      buf = (char *) NTL_MALLOC(iodigits, 1, 0);
      if (!buf) Error("out of memory");
   }

   long i = 0;

   while (d) {
      buf[i] = IntValToChar(d % 10);
      d = d / 10;
      i++;
   }

   if (justify) {
      long j = iodigits - i;
      while (j > 0) {
         s << "0";
         j--;
      }
   }

   while (i > 0) {
      i--;
      s << buf[i];
   }
}
      

   

ostream& operator<<(ostream& s, const ZZ& a)
{
   static ZZ b;
   static _ZZ_local_stack S;
   long r;
   long k;

   if (!iodigits) InitZZIO();

   b = a;

   k = sign(b);

   if (k == 0) {
      s << "0";
      return s;
   }

   if (k < 0) {
      s << "-";
      negate(b, b);
   }

   do {
      r = DivRem(b, b, ioradix);
      S.push(r);
   } while (!IsZero(b));

   r = S.pop();
   PrintDigits(s, r, 0);

   while (!S.empty()) {
      r = S.pop();
      PrintDigits(s, r, 1);
   }
      
   return s;
}



long GCD(long a, long b)
{
   long u, v, t, x;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) Error("GCD: integer overflow");
      a = -a;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) Error("GCD: integer overflow");
      b = -b;
   }


   if (b==0)
      x = a;
   else {
      u = a;
      v = b;
      do {
         t = u % v;
         u = v; 
         v = t;
      } while (v != 0);

      x = u;
   }

   return x;
}

         

void XGCD(long& d, long& s, long& t, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) Error("XGCD: integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) Error("XGCD: integer overflow");
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d = u;
   s = u1;
   t = v1;
}
   

long InvMod(long a, long n)
{
   long d, s, t;

   XGCD(d, s, t, a, n);
   if (d != 1) Error("InvMod: inverse undefined");
   if (s < 0)
      return s + n;
   else
      return s;
}


long PowerMod(long a, long ee, long n)
{
   long x, y;

   unsigned long e;

   if (ee < 0)
      e = - ((unsigned long) ee);
   else
      e = ee;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = MulMod(x, y, n);
      y = MulMod(y, y, n);
      e = e >> 1;
   }

   if (ee < 0) x = InvMod(x, n);

   return x;
}

long ProbPrime(long n, long NumTests)
{
   long m, x, y, z;
   long i, j, k;

   if (n <= 1) return 0;


   if (n == 2) return 1;
   if (n % 2 == 0) return 0;

   if (n == 3) return 1;
   if (n % 3 == 0) return 0;

   if (n == 5) return 1;
   if (n % 5 == 0) return 0;

   if (n == 7) return 1;
   if (n % 7 == 0) return 0;

   if (n >= NTL_SP_BOUND) {
      return ProbPrime(to_ZZ(n), NumTests);
   }

   m = n - 1;
   k = 0;
   while((m & 1) == 0) {
      m = m >> 1;
      k++;
   }

   // n - 1 == 2^k * m, m odd

   for (i = 0; i < NumTests; i++) {
      do {
         x = RandomBnd(n);
      } while (x == 0);
      // x == 0 is not a useful candidtae for a witness!


      if (x == 0) continue;
      z = PowerMod(x, m, n);
      if (z == 1) continue;
   
      j = 0;
      do {
         y = z;
         z = MulMod(y, y, n);
         j++;
      } while (j != k && z != 1);

      if (z != 1 || y !=  n-1) return 0;
   }

   return 1;
}


long MillerWitness(const ZZ& n, const ZZ& x)
{
   ZZ m, y, z;
   long j, k;

   if (x == 0) return 0;

   add(m, n, -1);
   k = MakeOdd(m);
   // n - 1 == 2^k * m, m odd

   PowerMod(z, x, m, n);
   if (z == 1) return 0;

   j = 0;
   do {
      y = z;
      SqrMod(z, y, n);
      j++;
   } while (j != k && z != 1);

   if (z != 1) return 1;
   add(y, y, 1);
   if (y != n) return 1;
   return 0;
}


// ComputePrimeBound computes a reasonable bound for trial
// division in the Miller-Rabin test.
// It is computed a bit on the "low" side, since being a bit
// low doesn't hurt much, but being too high can hurt a lot.

static
long ComputePrimeBound(long bn)
{
   long wn = (bn+NTL_ZZ_NBITS-1)/NTL_ZZ_NBITS;

   long fn;

   if (wn <= 36)
      fn = wn/4 + 1;
   else
      fn = long(1.67*sqrt(double(wn)));

   long prime_bnd;

   if (NumBits(bn) + NumBits(fn) > NTL_SP_NBITS)
      prime_bnd = NTL_SP_BOUND;
   else
      prime_bnd = bn*fn;

   return prime_bnd;
}


long ProbPrime(const ZZ& n, long NumTrials)
{
   if (n <= 1) return 0;

   if (n.SinglePrecision()) {
      return ProbPrime(to_long(n), NumTrials);
   }


   long prime_bnd = ComputePrimeBound(NumBits(n));


   PrimeSeq s;
   long p;

   p = s.next();
   while (p && p < prime_bnd) {
      if (rem(n, p) == 0)
         return 0;

      p = s.next();
   }

   ZZ W;
   W = 2;

   // first try W == 2....the exponentiation
   // algorithm runs slightly faster in this case

   if (MillerWitness(n, W))
      return 0;


   long i;

   for (i = 0; i < NumTrials; i++) {
      do {
         RandomBnd(W, n);
      } while (W == 0);
      // W == 0 is not a useful candidate for a witness!

      if (MillerWitness(n, W)) 
         return 0;
   }

   return 1;
}


void RandomPrime(ZZ& n, long l, long NumTrials)
{
   if (l <= 1)
      Error("RandomPrime: l out of range");

   if (l == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }

   do {
      RandomLen(n, l);
      if (!IsOdd(n)) add(n, n, 1);
   } while (!ProbPrime(n, NumTrials));
}

void NextPrime(ZZ& n, const ZZ& m, long NumTrials)
{
   ZZ x;

   if (m <= 2) {
      n = 2;
      return;
   }

   x = m;

   while (!ProbPrime(x, NumTrials))
      add(x, x, 1);

   n = x;
}

long NextPrime(long m, long NumTrials)
{
   long x;

   if (m <= 2) 
      return 2;

   x = m;

   while (x < NTL_SP_BOUND && !ProbPrime(x, NumTrials))
      x++;

   if (x >= NTL_SP_BOUND)
      Error("NextPrime: no more primes");

   return x;
}



long NextPowerOfTwo(long m)
{
   long k; 
   unsigned long n, um;

   if (m < 0) return 0;

   um = m;
   n = 1;
   k = 0;

   while (n < um) {
      n = n << 1;
      k++;
   }

   if (k >= NTL_BITS_PER_LONG-1)
      Error("NextPowerOfTwo: overflow");

   return k;
}



long NumBits(long a)
{
   unsigned long aa;
   if (a < 0) 
      aa = - ((unsigned long) a);
   else
      aa = a;

   long k = 0;
   while (aa) {
      k++;
      aa = aa >> 1;
   }

   return k;
}


long bit(long a, long k)
{
   unsigned long aa;
   if (a < 0)
      aa = - ((unsigned long) a);
   else
      aa = a;

   if (k < 0 || k >= NTL_BITS_PER_LONG) 
      return 0;
   else
      return long((aa >> k) & 1);
}



long divide(ZZ& q, const ZZ& a, const ZZ& b)
{
   static ZZ qq, r;

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

   DivRem(qq, r, a, b);
   if (!IsZero(r)) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, const ZZ& b)
{
   static ZZ r;

   if (IsZero(b)) return IsZero(a);
   if (IsOne(b)) return 1;

   rem(r, a, b);
   return IsZero(r);
}

long divide(ZZ& q, const ZZ& a, long b)
{
   static ZZ qq;

   if (!b) {
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

   long r = DivRem(qq, a, b);
   if (r) return 0;
   q = qq;
   return 1;
}

long divide(const ZZ& a, long b)
{
   if (!b) return IsZero(a);
   if (b == 1) {
      return 1;
   }

   long r = rem(a,  b);
   return (r == 0);
}
   


long RandomPrime_long(long l, long NumTrials)
{
   if (l <= 1 || l >= NTL_BITS_PER_LONG)
      Error("RandomPrime: length out of range");

   long n;
   do {
      n = RandomLen_long(l);
   } while (!ProbPrime(n, NumTrials));

   return n;
}


PrimeSeq::PrimeSeq()
{
   movesieve = 0;
   movesieve_mem = 0;
   pshift = -1;
   pindex = -1;
   exhausted = 0;
}

PrimeSeq::~PrimeSeq()
{
   if (movesieve_mem)
      free(movesieve_mem);
}

long PrimeSeq::next()
{
   if (exhausted) {
      return 0;
   }

   if (pshift < 0) {
      shift(0);
      return 2;
   }

   for (;;) {
      char *p = movesieve;
      long i = pindex;

      while ((++i) < NTL_PRIME_BND) {
         if (p[i]) {
            pindex = i;
            return pshift + 2 * i + 3;
         }
      }

      long newshift = pshift + 2*NTL_PRIME_BND;

      if (newshift > 2 * NTL_PRIME_BND * (2 * NTL_PRIME_BND + 1)) {
         /* end of the road */
         exhausted = 1;
         return 0;
      }

      shift(newshift);
   }
}

static char *lowsieve = 0;

void PrimeSeq::shift(long newshift)
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibound;
   char *p;

   if (!lowsieve)
      start();

   pindex = -1;
   exhausted = 0;

   if (newshift < 0) {
      pshift = -1;
      return;
   }

   if (newshift == pshift) return;

   pshift = newshift;

   if (pshift == 0) {
      movesieve = lowsieve;
   } 
   else {
      if (!movesieve_mem) {
         movesieve_mem = (char *) NTL_MALLOC(NTL_PRIME_BND, 1, 0);
         if (!movesieve_mem) 
            Error("out of memory in PrimeSeq");
      }

      p = movesieve = movesieve_mem;
      for (i = 0; i < NTL_PRIME_BND; i++)
         p[i] = 1;

      jstep = 3;
      ibound = pshift + 2 * NTL_PRIME_BND + 1;
      for (i = 0; jstep * jstep <= ibound; i++) {
         if (lowsieve[i]) {
            if (!((jstart = (pshift + 2) / jstep + 1) & 1))
               jstart++;
            if (jstart <= jstep)
               jstart = jstep;
            jstart = (jstart * jstep - pshift - 3) / 2;
            for (j = jstart; j < NTL_PRIME_BND; j += jstep)
               p[j] = 0;
         }
         jstep += 2;
      }
   }
}


void PrimeSeq::start()
{
   long i;
   long j;
   long jstep;
   long jstart;
   long ibnd;
   char *p;

   p = lowsieve = (char *) NTL_MALLOC(NTL_PRIME_BND, 1, 0);
   if (!p)
      Error("out of memory in PrimeSeq");

   for (i = 0; i < NTL_PRIME_BND; i++)
      p[i] = 1;
      
   jstep = 1;
   jstart = -1;
   ibnd = (SqrRoot(2 * NTL_PRIME_BND + 1) - 3) / 2;
   for (i = 0; i <= ibnd; i++) {
      jstart += 2 * ((jstep += 2) - 1);
      if (p[i])
         for (j = jstart; j < NTL_PRIME_BND; j += jstep)
            p[j] = 0;
   }
}

void PrimeSeq::reset(long b)
{
   if (b > (2*NTL_PRIME_BND+1)*(2*NTL_PRIME_BND+1)) {
      exhausted = 1;
      return;
   }

   if (b <= 2) {
      shift(-1);
      return;
   }

   if ((b & 1) == 0) b++;

   shift(((b-3) / (2*NTL_PRIME_BND))* (2*NTL_PRIME_BND));
   pindex = (b - pshift - 3)/2 - 1;
}
 
long Jacobi(const ZZ& aa, const ZZ& nn)
{
   ZZ a, n;
   long t, k;
   long d;

   a = aa;
   n = nn;
   t = 1;

   while (a != 0) {
      k = MakeOdd(a);
      d = trunc_long(n, 3);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if (trunc_long(a, 2) == 3 && (d & 3) == 3) t = -t;
      swap(a, n);
      rem(a, a, n);
   }

   if (n == 1)
      return t;
   else
      return 0;
}


void SqrRootMod(ZZ& x, const ZZ& aa, const ZZ& nn)
{
   if (aa == 0 || aa == 1) {
      x = aa;
      return;
   }

   // at this point, we must have nn >= 5

   if (trunc_long(nn, 2) == 3) {  // special case, n = 3 (mod 4)
      ZZ n, a, e, z;

      n = nn;
      a  = aa;

      add(e, n, 1);
      RightShift(e, e, 2);

      PowerMod(z, a, e, n);
      x = z;

      return;
   }

   ZZ n, m;
   int h, nlen;

   n = nn;
   nlen = NumBits(n);

   sub(m, n, 1);
   h = MakeOdd(m);  // h >= 2


   if (nlen > 50 && h < SqrRoot(nlen)) {
      long i, j;
      ZZ a, b, a_inv, c, r, m1, d;

      a = aa;
      InvMod(a_inv, a, n);

      if (h == 2) 
         b = 2;
      else {
         do {
            RandomBnd(b, n);
         } while (Jacobi(b, n) != -1);
      }


      PowerMod(c, b, m, n);
      
      add(m1, m, 1);
      RightShift(m1, m1, 1);
      PowerMod(r, a, m1, n);

      for (i = h-2; i >= 0; i--) {
         SqrMod(d, r, n);
         MulMod(d, d, a_inv, n);
         for (j = 0; j < i; j++)
            SqrMod(d, d, n);
         if (!IsOne(d))
            MulMod(r, r, c, n);
         SqrMod(c, c, n);
      } 

      x = r;
      return;
   } 





   long i, k;
   ZZ ma, t, u, v, e;
   ZZ t1, t2, t3, t4;

   n = nn;
   NegateMod(ma, aa, n);

   // find t such that t^2 - 4*a is not a square

   MulMod(t1, ma, 4, n);
   do {
      RandomBnd(t, n);
      SqrMod(t2, t, n);
      AddMod(t2, t2, t1, n);
   } while (Jacobi(t2, n) != -1);

   // compute u*X + v = X^{(n+1)/2} mod f, where f = X^2 - t*X + a

   add(e, n, 1);
   RightShift(e, e, 1);

   u = 0;
   v = 1;

   k = NumBits(e);

   for (i = k - 1; i >= 0; i--) {
      add(t2, u, v);
      sqr(t3, t2);  // t3 = (u+v)^2
      sqr(t1, u);
      sqr(t2, v);
      sub(t3, t3, t1);
      sub(t3, t3, t2); // t1 = u^2, t2 = v^2, t3 = 2*u*v
      rem(t1, t1, n);
      mul(t4, t1, t);
      add(t4, t4, t3);
      rem(u, t4, n);

      mul(t4, t1, ma);
      add(t4, t4, t2);
      rem(v, t4, n);
      
      if (bit(e, i)) {
         MulMod(t1, u, t, n);
         AddMod(t1, t1, v, n);
         MulMod(v, u, ma, n);
         u = t1;
      }

   }

   x = v;
}



// Chinese Remaindering.
//
// This version in new to v3.7, and is significantly
// simpler and faster than the previous version.
//
// This function takes as input g, a, G, p,
// such that a > 0, 0 <= G < p, and gcd(a, p) = 1.
// It computes a' = a*p and g' such that 
//   * g' = g (mod a);
//   * g' = G (mod p);
//   * -a'/2 < g' <= a'/2.
// It then sets g := g' and a := a', and returns 1 iff g has changed.
//
// Under normal use, the input value g satisfies -a/2 < g <= a/2;
// however, this was not documented or enforced in earlier versions,
// so to maintain backward compatability, no restrictions are placed
// on g.  This routine runs faster, though, if -a/2 < g <= a/2,
// and the first thing the routine does is to make this condition
// hold.
//
// Also, under normal use, both a and p are odd;  however, the routine
// will still work even if this is not so.
//
// The routine is based on the following simple fact.
//
// Let -a/2 < g <= a/2, and let h satisfy
//   * g + a h = G (mod p);
//   * -p/2 < h <= p/2.
// Further, if p = 2*h and g > 0, set
//   g' := g - a h;
// otherwise, set
//   g' := g + a h.
// Then g' so defined satisfies the above requirements.
//
// It is trivial to see that g's satisfies the congruence conditions.
// The only thing is to check that the "balancing" condition
// -a'/2 < g' <= a'/2 also holds.


long CRT(ZZ& gg, ZZ& a, long G, long p)
{
   if (p >= NTL_SP_BOUND) {
      ZZ GG, pp;
      conv(GG, G);
      conv(pp, p);
      return CRT(gg, a, GG, pp);
   }

   long modified = 0;

   static ZZ g;

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   long p1;
   p1 = p >> 1;

   long a_inv;
   a_inv = rem(a, p);
   a_inv = InvMod(a_inv, p);

   long h;
   h = rem(g, p);
   h = SubMod(G, h, p);
   h = MulMod(h, a_inv, p);
   if (h > p1)
      h = h - p;

   if (h != 0) {
      modified = 1;

      if (!(p & 1) && g > 0 && (h == p1))
         MulSubFrom(g, a, h);
      else
         MulAddTo(g, a, h);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}

long CRT(ZZ& gg, ZZ& a, const ZZ& G, const ZZ& p)
{
   long modified = 0;

   ZZ g;

   if (!CRTInRange(gg, a)) {
      modified = 1;
      ZZ a1;
      rem(g, gg, a);
      RightShift(a1, a, 1);
      if (g > a1) sub(g, g, a);
   }
   else
      g = gg;


   ZZ p1;
   RightShift(p1, p, 1);

   ZZ a_inv;
   rem(a_inv, a, p);
   InvMod(a_inv, a_inv, p);

   ZZ h;
   rem(h, g, p);
   SubMod(h, G, h, p);
   MulMod(h, h, a_inv, p);
   if (h > p1)
      sub(h, h, p);

   if (h != 0) {
      modified = 1;
      ZZ ah;
      mul(ah, a, h);

      if (!IsOdd(p) && g > 0 &&  (h == p1))
         sub(g, g, ah);
      else
         add(g, g, ah);
   }

   mul(a, a, p);
   gg = g;

   return modified;
}



void sub(ZZ& x, const ZZ& a, long b)
{
   static ZZ B;
   conv(B, b);
   sub(x, a, B);
}

void sub(ZZ& x, long a, const ZZ& b)
{
   static ZZ A;
   conv(A, a);
   sub(x, A, b);
}


void power2(ZZ& x, long e)
{
   if (e < 0) Error("power2: negative exponent");
   set(x);
   LeftShift(x, x, e);
}

   
void conv(ZZ& x, const char *s)
{
   long c;
   long cval;
   long sign;
   long ndigits;
   long acc;
   long i = 0;

   static ZZ a;

   if (!s) Error("bad ZZ input");

   if (!iodigits) InitZZIO();

   a = 0;

   c = s[i];
   while (IsWhiteSpace(c)) {
      i++;
      c = s[i];
   }

   if (c == '-') {
      sign = -1;
      i++;
      c = s[i];
   }
   else
      sign = 1;

   cval = CharToIntVal(c);
   if (cval < 0 || cval > 9) Error("bad ZZ input");

   ndigits = 0;
   acc = 0;
   while (cval >= 0 && cval <= 9) {
      acc = acc*10 + cval;
      ndigits++;

      if (ndigits == iodigits) {
         mul(a, a, ioradix);
         add(a, a, acc);
         ndigits = 0;
         acc = 0;
      }

      i++;
      c = s[i];
      cval = CharToIntVal(c);
   }

   if (ndigits != 0) {
      long mpy = 1;
      while (ndigits > 0) {
         mpy = mpy * 10;
         ndigits--;
      }

      mul(a, a, mpy);
      add(a, a, acc);
   }

   if (sign == -1)
      negate(a, a);

   x = a;
}



void bit_and(ZZ& x, const ZZ& a, long b)
{
   static ZZ B;
   conv(B, b);
   bit_and(x, a, B);
}

void bit_or(ZZ& x, const ZZ& a, long b)
{
   static ZZ B;
   conv(B, b);
   bit_or(x, a, B);
}

void bit_xor(ZZ& x, const ZZ& a, long b)
{
   static ZZ B;
   conv(B, b);
   bit_xor(x, a, B);
}


long power_long(long a, long e)
{
   if (e < 0) Error("power_long: negative exponent");

   if (e == 0) return 1;

   if (a == 1) return 1;
   if (a == -1) {
      if (e & 1)
         return -1;
      else
         return 1;
   }

   // no overflow check --- result is computed correctly
   // modulo word size

   unsigned long res = 1;
   unsigned long aa = a;
   long i;

   for (i = 0; i < e; i++)
      res *= aa;

   return to_long(res);
}

//  RANDOM NUMBER GENERATION

// Idea for this PRNG.  Iteratively hash seed using md5 
// to get 256 bytes to initialize arc4.
// Then use arc4 to get a pseudo-random byte stream.

// I've taken care that the pseudo-random numbers generated by
// the routines RandomBnd, RandomBits, and RandomLen 
// are completely platform independent.

// I make use of the md5 compression function,
// which I've modified to work on 64-bit machines


/*
 *  BEGIN RSA's md5 stuff
 *
 */

/*
 **********************************************************************
 ** md5.c                                                            **
 ** RSA Data Security, Inc. MD5 Message Digest Algorithm             **
 ** Created: 2/17/90 RLR                                             **
 ** Revised: 1/91 SRD,AJ,BSK,JT Reference C Version                  **
 **********************************************************************
 */

/*
 **********************************************************************
 ** Copyright (C) 1990, RSA Data Security, Inc. All rights reserved. **
 **                                                                  **
 ** License to copy and use this software is granted provided that   **
 ** it is identified as the "RSA Data Security, Inc. MD5 Message     **
 ** Digest Algorithm" in all material mentioning or referencing this **
 ** software or this function.                                       **
 **                                                                  **
 ** License is also granted to make and use derivative works         **
 ** provided that such works are identified as "derived from the RSA **
 ** Data Security, Inc. MD5 Message Digest Algorithm" in all         **
 ** material mentioning or referencing the derived work.             **
 **                                                                  **
 ** RSA Data Security, Inc. makes no representations concerning      **
 ** either the merchantability of this software or the suitability   **
 ** of this software for any particular purpose.  It is provided "as **
 ** is" without express or implied warranty of any kind.             **
 **                                                                  **
 ** These notices must be retained in any copies of any part of this **
 ** documentation and/or software.                                   **
 **********************************************************************
 */


#if (NTL_BITS_PER_LONG <= 32)
#define TRUNC32(x) (x)
#else
#define TRUNC32(x) ((x) & ((1UL << 32)-1UL))
#endif

/* F, G and H are basic MD5 functions: selection, majority, parity */
#define F(x, y, z) (((x) & (y)) | ((~x) & (z)))
#define G(x, y, z) (((x) & (z)) | ((y) & (~z)))
#define H(x, y, z) ((x) ^ (y) ^ (z))
#define I(x, y, z) (TRUNC32((y) ^ ((x) | (~z)))) 

/* ROTATE_LEFT rotates x left n bits */
#define ROTATE_LEFT(x, n) (TRUNC32(((x) << (n)) | ((x) >> (32-(n)))))

/* FF, GG, HH, and II transformations for rounds 1, 2, 3, and 4 */
/* Rotation is separate from addition to prevent recomputation */
#define FF(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + F((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }
#define GG(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + G((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }
#define HH(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + H((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }
#define II(a, b, c, d, x, s, ac) \
  {(a) = TRUNC32((a) + I((b), (c), (d)) + (x) + (ac)); \
   (a) = ROTATE_LEFT((a), (s)); \
   (a) = TRUNC32((a) + (b)); \
  }



static
void MD5_default_IV(unsigned long *buf)
{
   buf[0] = 0x67452301UL;
   buf[1] = 0xefcdab89UL;
   buf[2] = 0x98badcfeUL;
   buf[3] = 0x10325476UL;
}



/* Basic MD5 step. Transform buf based on in.
 */

static
void MD5_compress(unsigned long *buf, unsigned long *in)
{
  unsigned long a = buf[0], b = buf[1], c = buf[2], d = buf[3];

  /* Round 1 */
#define S11 7
#define S12 12
#define S13 17
#define S14 22
  FF ( a, b, c, d, in[ 0], S11, 3614090360UL); /* 1 */
  FF ( d, a, b, c, in[ 1], S12, 3905402710UL); /* 2 */
  FF ( c, d, a, b, in[ 2], S13,  606105819UL); /* 3 */
  FF ( b, c, d, a, in[ 3], S14, 3250441966UL); /* 4 */
  FF ( a, b, c, d, in[ 4], S11, 4118548399UL); /* 5 */
  FF ( d, a, b, c, in[ 5], S12, 1200080426UL); /* 6 */
  FF ( c, d, a, b, in[ 6], S13, 2821735955UL); /* 7 */
  FF ( b, c, d, a, in[ 7], S14, 4249261313UL); /* 8 */
  FF ( a, b, c, d, in[ 8], S11, 1770035416UL); /* 9 */
  FF ( d, a, b, c, in[ 9], S12, 2336552879UL); /* 10 */
  FF ( c, d, a, b, in[10], S13, 4294925233UL); /* 11 */
  FF ( b, c, d, a, in[11], S14, 2304563134UL); /* 12 */
  FF ( a, b, c, d, in[12], S11, 1804603682UL); /* 13 */
  FF ( d, a, b, c, in[13], S12, 4254626195UL); /* 14 */
  FF ( c, d, a, b, in[14], S13, 2792965006UL); /* 15 */
  FF ( b, c, d, a, in[15], S14, 1236535329UL); /* 16 */

  /* Round 2 */
#define S21 5
#define S22 9
#define S23 14
#define S24 20
  GG ( a, b, c, d, in[ 1], S21, 4129170786UL); /* 17 */
  GG ( d, a, b, c, in[ 6], S22, 3225465664UL); /* 18 */
  GG ( c, d, a, b, in[11], S23,  643717713UL); /* 19 */
  GG ( b, c, d, a, in[ 0], S24, 3921069994UL); /* 20 */
  GG ( a, b, c, d, in[ 5], S21, 3593408605UL); /* 21 */
  GG ( d, a, b, c, in[10], S22,   38016083UL); /* 22 */
  GG ( c, d, a, b, in[15], S23, 3634488961UL); /* 23 */
  GG ( b, c, d, a, in[ 4], S24, 3889429448UL); /* 24 */
  GG ( a, b, c, d, in[ 9], S21,  568446438UL); /* 25 */
  GG ( d, a, b, c, in[14], S22, 3275163606UL); /* 26 */
  GG ( c, d, a, b, in[ 3], S23, 4107603335UL); /* 27 */
  GG ( b, c, d, a, in[ 8], S24, 1163531501UL); /* 28 */
  GG ( a, b, c, d, in[13], S21, 2850285829UL); /* 29 */
  GG ( d, a, b, c, in[ 2], S22, 4243563512UL); /* 30 */
  GG ( c, d, a, b, in[ 7], S23, 1735328473UL); /* 31 */
  GG ( b, c, d, a, in[12], S24, 2368359562UL); /* 32 */

  /* Round 3 */
#define S31 4
#define S32 11
#define S33 16
#define S34 23
  HH ( a, b, c, d, in[ 5], S31, 4294588738UL); /* 33 */
  HH ( d, a, b, c, in[ 8], S32, 2272392833UL); /* 34 */
  HH ( c, d, a, b, in[11], S33, 1839030562UL); /* 35 */
  HH ( b, c, d, a, in[14], S34, 4259657740UL); /* 36 */
  HH ( a, b, c, d, in[ 1], S31, 2763975236UL); /* 37 */
  HH ( d, a, b, c, in[ 4], S32, 1272893353UL); /* 38 */
  HH ( c, d, a, b, in[ 7], S33, 4139469664UL); /* 39 */
  HH ( b, c, d, a, in[10], S34, 3200236656UL); /* 40 */
  HH ( a, b, c, d, in[13], S31,  681279174UL); /* 41 */
  HH ( d, a, b, c, in[ 0], S32, 3936430074UL); /* 42 */
  HH ( c, d, a, b, in[ 3], S33, 3572445317UL); /* 43 */
  HH ( b, c, d, a, in[ 6], S34,   76029189UL); /* 44 */
  HH ( a, b, c, d, in[ 9], S31, 3654602809UL); /* 45 */
  HH ( d, a, b, c, in[12], S32, 3873151461UL); /* 46 */
  HH ( c, d, a, b, in[15], S33,  530742520UL); /* 47 */
  HH ( b, c, d, a, in[ 2], S34, 3299628645UL); /* 48 */

  /* Round 4 */
#define S41 6
#define S42 10
#define S43 15
#define S44 21
  II ( a, b, c, d, in[ 0], S41, 4096336452UL); /* 49 */
  II ( d, a, b, c, in[ 7], S42, 1126891415UL); /* 50 */
  II ( c, d, a, b, in[14], S43, 2878612391UL); /* 51 */
  II ( b, c, d, a, in[ 5], S44, 4237533241UL); /* 52 */
  II ( a, b, c, d, in[12], S41, 1700485571UL); /* 53 */
  II ( d, a, b, c, in[ 3], S42, 2399980690UL); /* 54 */
  II ( c, d, a, b, in[10], S43, 4293915773UL); /* 55 */
  II ( b, c, d, a, in[ 1], S44, 2240044497UL); /* 56 */
  II ( a, b, c, d, in[ 8], S41, 1873313359UL); /* 57 */
  II ( d, a, b, c, in[15], S42, 4264355552UL); /* 58 */
  II ( c, d, a, b, in[ 6], S43, 2734768916UL); /* 59 */
  II ( b, c, d, a, in[13], S44, 1309151649UL); /* 60 */
  II ( a, b, c, d, in[ 4], S41, 4149444226UL); /* 61 */
  II ( d, a, b, c, in[11], S42, 3174756917UL); /* 62 */
  II ( c, d, a, b, in[ 2], S43,  718787259UL); /* 63 */
  II ( b, c, d, a, in[ 9], S44, 3951481745UL); /* 64 */

  buf[0] = TRUNC32(buf[0] + a);
  buf[1] = TRUNC32(buf[1] + b);
  buf[2] = TRUNC32(buf[2] + c);
  buf[3] = TRUNC32(buf[3] + d);
}


/*
 *  END RSA's md5 stuff
 *
 */


static
void words_from_bytes(unsigned long *txtl, unsigned char *txtc, long n)
{
   long i;
   unsigned long v;

   for (i = 0; i < n; i++) {
      v = txtc[4*i];
      v += ((unsigned long) (txtc[4*i+1])) << 8;
      v += ((unsigned long) (txtc[4*i+2])) << 16;
      v += ((unsigned long) (txtc[4*i+3])) << 24;
      txtl[i] = v;
   }
}

static 
void bytes_from_words(unsigned char *txtc, unsigned long *txtl, long n)
{
   long i;
   unsigned long v;

   for (i = 0; i < n; i++) {
      v = txtl[i];
      txtc[4*i] = v & 255;
      v = v >> 8;
      txtc[4*i+1] = v & 255;
      v = v >> 8;
      txtc[4*i+2] = v & 255;
      v = v >> 8;
      txtc[4*i+3] = v & 255;
   }
}


static
void MD5_compress1(unsigned long *buf, unsigned char *in, long n)
{
   unsigned long txtl[16];
   unsigned char txtc[64]; 
   long i, j, k;

   if (n < 0) n = 0;

   i = 0;
   while (i < n) {
      k = n-i;
      if (k > 64) k = 64;
      for (j = 0; j < k; j++)
         txtc[j] = in[i+j];
      for (; j < 64; j++)
         txtc[j] = 0;
      words_from_bytes(txtl, txtc, 16);
      MD5_compress(buf, txtl);
      i += k;
   }
}


// the "cipherpunk" version of arc4 

struct _ZZ_arc4_key
{      
    unsigned char state[256];       
    unsigned char x;        
    unsigned char y;
};


static inline
void swap_byte(unsigned char *a, unsigned char *b)
{
    unsigned char swapByte; 
    
    swapByte = *a; 
    *a = *b;      
    *b = swapByte;
}

static
void prepare_key(unsigned char *key_data_ptr, 
                 long key_data_len, _ZZ_arc4_key *key)
{
    unsigned char index1;
    unsigned char index2;
    unsigned char* state;
    long counter;     
    
    state = &key->state[0];         
    for(counter = 0; counter < 256; counter++)              
       state[counter] = counter;               
    key->x = 0;     
    key->y = 0;     
    index1 = 0;     
    index2 = 0;             
    for(counter = 0; counter < 256; counter++)      
    {               
         index2 = (key_data_ptr[index1] + state[counter] + index2) & 255;                
         swap_byte(&state[counter], &state[index2]);            

         index1 = (index1 + 1) % key_data_len;  
    }       
}



static
void arc4(unsigned char *buffer_ptr, long buffer_len, _ZZ_arc4_key *key)
{ 
    unsigned char x;
    unsigned char y;
    unsigned char* state;
    unsigned char xorIndex;
    long counter;              
    
    x = key->x;     
    y = key->y;     
    
    state = &key->state[0];         
    for(counter = 0; counter < buffer_len; counter ++)      
    {               
         x = (x + 1) & 255;
         y = (state[x] + y) & 255;
         swap_byte(&state[x], &state[y]);                        
              
         xorIndex = (state[x] + state[y]) & 255;
              
         buffer_ptr[counter] = state[xorIndex];         
     }               
     key->x = x;     
     key->y = y;
}

// global state information for PRNG

static long ran_initialized = 0;
static _ZZ_arc4_key ran_key;

static unsigned long default_md5_tab[16] = {
744663023UL, 1011602954UL, 3163087192UL, 3383838527UL, 
3305324122UL, 3197458079UL, 2266495600UL, 2760303563UL, 
346234297UL, 1919920720UL, 1896169861UL, 2192176675UL, 
2027150322UL, 2090160759UL, 2134858730UL, 1131796244UL
};



static
void build_arc4_tab(unsigned char *seed_bytes, const ZZ& s)
{
   long nb = NumBytes(s);
   
   unsigned char *txt;

   typedef unsigned char u_char;
   txt = NTL_NEW_OP u_char[nb + 68];
   if (!txt) Error("out of memory");

   BytesFromZZ(txt + 4, s, nb);

   bytes_from_words(txt + nb + 4, default_md5_tab, 16);

   unsigned long buf[4];

   unsigned long i;
   for (i = 0; i < 16; i++) {
      MD5_default_IV(buf);
      bytes_from_words(txt, &i, 1);

      MD5_compress1(buf, txt, nb + 68);

      bytes_from_words(seed_bytes + 16*i, buf, 4);
   }

   delete [] txt;
}


void SetSeed(const ZZ& s)
{
   unsigned char seed_bytes[256];

   build_arc4_tab(seed_bytes, s);
   prepare_key(seed_bytes, 256, &ran_key);

   // we discard the first 1024 bytes of the arc4 stream, as this is
   // recommended practice.

   arc4(seed_bytes, 256, &ran_key);
   arc4(seed_bytes, 256, &ran_key);
   arc4(seed_bytes, 256, &ran_key);
   arc4(seed_bytes, 256, &ran_key);

   ran_initialized = 1;
}

static 
void ran_bytes(unsigned char *bytes, long n)
{
   if (!ran_initialized) SetSeed(ZZ::zero());
   arc4(bytes, n, &ran_key);
}


unsigned long RandomWord()
{
   unsigned char buf[NTL_BITS_PER_LONG/8];
   long i;
   unsigned long res;

   ran_bytes(buf, NTL_BITS_PER_LONG/8);

   res = 0;
   for (i = NTL_BITS_PER_LONG/8 - 1; i >= 0; i--) {
      res = res << 8;
      res = res | buf[i];
   }

   return res;
}

long RandomBits_long(long l)
{
   if (l <= 0) return 0;
   if (l >= NTL_BITS_PER_LONG) 
      Error("RandomBits: length too big");

   unsigned char buf[NTL_BITS_PER_LONG/8];
   unsigned long res;
   long i;

   long nb = (l+7)/8;
   ran_bytes(buf, nb);

   res = 0;
   for (i = nb - 1; i >= 0; i--) {
      res = res << 8;
      res = res | buf[i];
   }

   return long(res & ((1UL << l)-1UL)); 
}

unsigned long RandomBits_ulong(long l)
{
   if (l <= 0) return 0;
   if (l > NTL_BITS_PER_LONG) 
      Error("RandomBits: length too big");

   unsigned char buf[NTL_BITS_PER_LONG/8];
   unsigned long res;
   long i;

   long nb = (l+7)/8;
   ran_bytes(buf, nb);

   res = 0;
   for (i = nb - 1; i >= 0; i--) {
      res = res << 8;
      res = res | buf[i];
   }

   if (l < NTL_BITS_PER_LONG)
      res = res & ((1UL << l)-1UL);

   return res;
}

long RandomLen_long(long l)
{
   if (l <= 0) return 0;
   if (l == 1) return 1;
   if (l >= NTL_BITS_PER_LONG) 
      Error("RandomLen: length too big");

   return RandomBits_long(l-1) + (1L << (l-1)); 
}


void RandomBits(ZZ& x, long l)
{
   if (l <= 0) {
      x = 0;
      return;
   }

   if (NTL_OVERFLOW(l, 1, 0))
      Error("RandomBits: length too big");

   long nb = (l+7)/8;

   static unsigned char *buf = 0;
   static long buf_len = 0;

   if (nb > buf_len) {
      if (buf) delete [] buf;
      buf_len = ((nb + 1023)/1024)*1024; // allocate in 1024-byte lots
      typedef unsigned char u_char;
      buf = NTL_NEW_OP u_char[buf_len];
      if (!buf) Error("out of memory");
   }

   ran_bytes(buf, nb);

   static ZZ res;

   ZZFromBytes(res, buf, nb);
   trunc(res, res, l);

   x = res;
}

void RandomLen(ZZ& x, long l)
{
   if (l <= 0) {
      x = 0;
      return;
   }

   if (l == 1) {
      x = 1;
      return;
   }

   if (NTL_OVERFLOW(l, 1, 0))
      Error("RandomLen: length too big");

   // pre-allocate space to avoid two allocations
   long nw = (l + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS;
   x.SetSize(nw);

   RandomBits(x, l-1);
   SetBit(x, l-1);
}


const long RandomBndExcess = 8;


void RandomBnd(ZZ& x, const ZZ& bnd)
{
   if (bnd <= 1) {
      x = 0;
      return;
   }

   long k = NumBits(bnd);

   if (weight(bnd) == 1) {
      RandomBits(x, k-1);
      return;
   }

   long l = k + RandomBndExcess;

   static ZZ t, r, t1;

   do {
      RandomBits(t, l);
      rem(r, t, bnd);
      sub(t1, bnd, r);
      add(t, t, t1);
   } while (NumBits(t) > l);

   x = r;
}

long RandomBnd(long bnd)
{
   if (bnd <= 1) return 0;

   long k = NumBits(bnd);

   if (((bnd - 1) & bnd) == 0) 
      return RandomBits_long(k-1);

   long l = k + RandomBndExcess;

   if (l > NTL_BITS_PER_LONG-2) {
      static ZZ Bnd, res;

      Bnd = bnd;
      RandomBnd(res, Bnd);
      return to_long(res);
   }

   long t, r;

   do {
      t = RandomBits_long(l);
      r = t % bnd;
   } while (t + bnd - r > (1L << l)); 

   return r;
}




// More prime generation stuff...

static
double Log2(double x)
{
   static double log2 = log(2.0);
   return log(x)/log2;
}

// Define p(k,t) to be the conditional probability that a random, odd, k-bit 
// number is composite, given that it passes t iterations of the 
// Miller-Rabin test.
// This routine returns 0 or 1, and if it returns 1 then
// p(k,t) <= 2^{-n}.
// This basically encodes the estimates of Damgard, Landrock, and Pomerance;
// it uses floating point arithmetic, but is coded in such a way
// that its results should be correct, assuming that the log function
// is computed with reasonable precision.
// 
// It is assumed that k >= 3 and t >= 1; if this does not hold,
// then 0 is returned.

static
long ErrBoundTest(long kk, long tt, long nn)

{
   const double fudge = (1.0 + 1024.0/NTL_FDOUBLE_PRECISION);
   const double log2_3 = Log2(3.0);
   const double log2_7 = Log2(7.0);
   const double log2_20 = Log2(20.0);

   double k = kk;
   double t = tt;
   double n = nn;

   if (k < 3 || t < 1) return 0;
   if (n < 1) return 1;

   // the following test is largely academic
   if (9*t > NTL_FDOUBLE_PRECISION) Error("ErrBoundTest: t too big");

   double log2_k = Log2(k);

   if ((n + log2_k)*fudge <= 2*t)
      return 1;

   if ((2*log2_k + 4.0 + n)*fudge <= 2*sqrt(k))
      return 2;

   if ((t == 2 && k >= 88) || (3 <= t && 9*t <= k && k >= 21)) {
      if ((1.5*log2_k + t + 4.0 + n)*fudge <= 0.5*Log2(t) + 2*(sqrt(t*k)))
         return 3;
   }

   if (k <= 9*t && 4*t <= k && k >= 21) {
      if ( ((log2_3 + log2_7 + log2_k + n)*fudge <= log2_20 + 5*t)  &&
           ((log2_3 + (15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t) &&
           ((2*log2_3 + 2 + log2_k + n)*fudge <= k/4 + 3*t) )
         return 4; 
   }

   if (4*t >= k && k >= 21) {
      if (((15.0/4.0)*log2_k + n)*fudge <= log2_7 + k/2 + 2*t)
         return 5;
   }

   return 0;
}


void GenPrime(ZZ& n, long k, long err)
{
   if (k <= 1) Error("GenPrime: bad length");

   if (k > (1L << 20)) Error("GenPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }


   long t;

   t = 1;
   while (!ErrBoundTest(k, t, err))
      t++;

   RandomPrime(n, k, t);
}


long GenPrime_long(long k, long err)
{
   if (k <= 1) Error("GenPrime: bad length");

   if (k >= NTL_BITS_PER_LONG) Error("GenPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         return 3;
      else
         return 2;
   }

   long t;

   t = 1;
   while (!ErrBoundTest(k, t, err))
      t++;

   return RandomPrime_long(k, t);
}


void GenGermainPrime(ZZ& n, long k, long err)
{
   if (k <= 1) Error("GenGermainPrime: bad length");

   if (k > (1L << 20)) Error("GenGermainPrime: length too large");

   if (err < 1) err = 1;
   if (err > 512) err = 512;

   if (k == 2) {
      if (RandomBnd(2))
         n = 3;
      else
         n = 2;

      return;
   }


   long prime_bnd = ComputePrimeBound(k);

   if (NumBits(prime_bnd) >= k/2)
      prime_bnd = (1L << (k/2-1));


   ZZ two;
   two = 2;

   ZZ n1;

   
   PrimeSeq s;

   ZZ iter;
   iter = 0;


   for (;;) {
      iter++;

      RandomLen(n, k);
      if (!IsOdd(n)) add(n, n, 1);

      s.reset(3);
      long p;

      long sieve_passed = 1;

      p = s.next();
      while (p && p < prime_bnd) {
         long r = rem(n, p);

         if (r == 0) {
            sieve_passed = 0;
            break;
         }

         // test if 2*r + 1 = 0 (mod p)
         if (r == p-r-1) {
            sieve_passed = 0;
            break;
         }

         p = s.next();
      }

      if (!sieve_passed) continue;


      if (MillerWitness(n, two)) continue;

      // n1 = 2*n+1
      mul(n1, n, 2);
      add(n1, n1, 1);


      if (MillerWitness(n1, two)) continue;

      // now do t M-R iterations...just to make sure
 
      // First compute the appropriate number of M-R iterations, t
      // The following computes t such that 
      //       p(k,t)*8/k <= 2^{-err}/(5*iter^{1.25})
      // which suffices to get an overall error probability of 2^{-err}.
      // Note that this method has the advantage of not requiring 
      // any assumptions on the density of Germain primes.

      long err1 = max(1, err + 7 + (5*NumBits(iter) + 3)/4 - NumBits(k));
      long t;
      t = 1;
      while (!ErrBoundTest(k, t, err1))
         t++;

      ZZ W;
      long MR_passed = 1;

      long i;
      for (i = 1; i <= t; i++) {
         do {
            RandomBnd(W, n);
         } while (W == 0);
         // W == 0 is not a useful candidate witness!

         if (MillerWitness(n, W)) {
            MR_passed = 0;
            break;
         }
      }

      if (MR_passed) break;
   }
}

long GenGermainPrime_long(long k, long err)
{
   if (k >= NTL_BITS_PER_LONG-1)
      Error("GenGermainPrime_long: length too long");

   ZZ n;
   GenGermainPrime(n, k, err);
   return to_long(n);
}


NTL_END_IMPL
