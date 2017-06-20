
#include <NTL/xdouble.h>
#include <NTL/RR.h>


#include <NTL/new.h>

NTL_START_IMPL



long xdouble::oprec = 10;

void xdouble::SetOutputPrecision(long p)
{
   if (p < 1) p = 1;

   if (NTL_OVERFLOW(p, 1, 0)) 
      Error("xdouble: output precision too big");

   oprec = p;
}

void xdouble::normalize() 
{
   if (x == 0) 
      e = 0;
   else if (x > 0) {
      while (x < NTL_XD_HBOUND_INV) { x *= NTL_XD_BOUND; e--; }
      while (x > NTL_XD_HBOUND) { x *= NTL_XD_BOUND_INV; e++; }
   }
   else {
      while (x > -NTL_XD_HBOUND_INV) { x *= NTL_XD_BOUND; e--; }
      while (x < -NTL_XD_HBOUND) { x *= NTL_XD_BOUND_INV; e++; }
   }

   if (e >= NTL_OVFBND)
      Error("xdouble: overflow");

   if (e <= -NTL_OVFBND)
      Error("xdouble: underflow");
}
   


xdouble to_xdouble(double a)
{
   if (a == 0 || a == 1 || (a > 0 && a >= NTL_XD_HBOUND_INV && a <= NTL_XD_HBOUND)
       || (a < 0 && a <= -NTL_XD_HBOUND_INV && a >= -NTL_XD_HBOUND)) {
      
      return xdouble(a, 0); 

   }

   if (!IsFinite(&a))
      Error("double to xdouble conversion: non finite value");

   xdouble z = xdouble(a, 0);
   z.normalize();
   return z;
}


void conv(double& xx, const xdouble& a)
{
   double x;
   long e;

   x = a.x;
   e = a.e;

   while (e > 0) { x *= NTL_XD_BOUND; e--; }
   while (e < 0) { x *= NTL_XD_BOUND_INV; e++; }

   xx = x;
}




xdouble operator+(const xdouble& a, const xdouble& b)
{
   xdouble z;

   if (a.x == 0) 
      return b;

   if (b.x == 0)
     return a;
      

   if (a.e == b.e) {
      z.x = a.x + b.x;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else if (a.e > b.e) {
      if (a.e > b.e+1)
         return a;

      z.x = a.x + b.x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else {
      if (b.e > a.e+1)
         return b;

      z.x = a.x*NTL_XD_BOUND_INV + b.x;
      z.e = b.e;
      z.normalize();
      return z;
   }
}


xdouble operator-(const xdouble& a, const xdouble& b)
{
   xdouble z;

   if (a.x == 0)
      return -b;

   if (b.x == 0)
      return a;

   if (a.e == b.e) {
      z.x = a.x - b.x;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else if (a.e > b.e) {
      if (a.e > b.e+1)
         return a;

      z.x = a.x - b.x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return z;
   }
   else {
      if (b.e > a.e+1)
         return -b;

      z.x = a.x*NTL_XD_BOUND_INV - b.x;
      z.e = b.e;
      z.normalize();
      return z;
   }
}

xdouble operator-(const xdouble& a)
{
   xdouble z;
   z.x = -a.x;
   z.e = a.e;
   return z;
}

xdouble operator*(const xdouble& a, const xdouble& b)
{
   xdouble z;

   z.e = a.e + b.e;
   z.x = a.x * b.x;
   z.normalize();
   return z;
}

xdouble operator/(const xdouble& a, const xdouble& b)
{
   xdouble z;

   if (b.x == 0) Error("xdouble division by 0");

   z.e = a.e - b.e;
   z.x = a.x / b.x;
   z.normalize();
   return z;
}



long compare(const xdouble& a, const xdouble& b)
{
   xdouble z = a - b;

   if (z.x < 0)
      return -1;
   else if (z.x == 0)
      return 0;
   else
      return 1;
}

long sign(const xdouble& z)
{
   if (z.x < 0)
      return -1;
   else if (z.x == 0)
      return 0;
   else
      return 1;
}
   


xdouble trunc(const xdouble& a)
{
   if (a.x >= 0)
      return floor(a);
   else
      return ceil(a);
}


xdouble floor(const xdouble& aa)
{
   xdouble z;

   xdouble a = aa;
   ForceToMem(&a.x);

   if (a.e == 0) {
      z.x = floor(a.x);
      z.e = 0;
      z.normalize();
      return z;
   }
   else if (a.e > 0) {
      return a;
   }
   else {
      if (a.x < 0)
         return to_xdouble(-1);
      else
         return to_xdouble(0);
   }
}

xdouble ceil(const xdouble& aa)
{
   xdouble z;

   xdouble a = aa;
   ForceToMem(&a.x);

   if (a.e == 0) {
      z.x = ceil(a.x);
      z.e = 0;
      z.normalize();
      return z;
   }
   else if (a.e > 0) {
      return a;
   }
   else {
      if (a.x < 0)
         return to_xdouble(0);
      else
         return to_xdouble(1);
   }
}

xdouble to_xdouble(const ZZ& a)
{
   long old_p = RR::precision();
   RR::SetPrecision(NTL_DOUBLE_PRECISION);
   
   static RR t;
   conv(t, a);

   double x;
   conv(x, t.mantissa());

   xdouble y, z, res;

   conv(y, x);
   power2(z, t.exponent());

   res = y*z;

   RR::SetPrecision(old_p);

   return res;
}

void conv(ZZ& x, const xdouble& a)
{
   xdouble b = floor(a);
   long old_p = RR::precision();
   RR::SetPrecision(NTL_DOUBLE_PRECISION);
   static RR t;
   conv(t, b);
   conv(x, t);
   RR::SetPrecision(old_p);
}


xdouble fabs(const xdouble& a)
{
   xdouble z;

   z.e = a.e;
   z.x = fabs(a.x);
   return z;
}

xdouble sqrt(const xdouble& a)
{
   if (a == 0)
      return to_xdouble(0);

   if (a < 0)
      Error("xdouble: sqrt of negative number");

   xdouble t;

   if (a.e & 1) {
      t.e = (a.e - 1)/2;
      t.x = sqrt(a.x * NTL_XD_BOUND);
   }
   else {
      t.e = a.e/2;
      t.x = sqrt(a.x);
   }

   t.normalize();

   return t;
}
      

void power(xdouble& z, const xdouble& a, const ZZ& e)
{
   xdouble b, res;

   b = a;

   res = 1;
   long n = NumBits(e);
   long i;

   for (i = n-1; i >= 0; i--) {
      res = res * res;
      if (bit(e, i))
         res = res * b;
   }

   if (sign(e) < 0) 
      z = 1/res;
   else
      z = res;
}




void power(xdouble& z, const xdouble& a, long e)
{
   static ZZ E;
   E = e;
   power(z, a, E);
}
   

   


void power2(xdouble& z, long e)
{
   long hb = NTL_XD_HBOUND_LOG;
   long b = 2*hb;

   long q, r;

   q = e/b;
   r = e%b;

   while (r >= hb) {
      r -= b;
      q++;
   }

   while (r < -hb) {
      r += b;
      q--;
   }

   if (q >= NTL_OVFBND)
      Error("xdouble: overflow");

   if (q <= -NTL_OVFBND)
      Error("xdouble: underflow");

   double x = _ntl_ldexp(1.0, r);

   z.x = x;
   z.e = q;
}


void MulAdd(xdouble& z, const xdouble& a, const xdouble& b, const xdouble& c)
// z = a + b*c
{
   double x;
   long e;

   e = b.e + c.e;
   x = b.x * c.x;

   if (x == 0) { 
      z = a;
      return;
   }

   if (a.x == 0) {
      z.e = e;
      z.x = x;
      z.normalize();
      return;
   }
      

   if (a.e == e) {
      z.x = a.x + x;
      z.e = e;
      z.normalize();
      return;
   }
   else if (a.e > e) {
      if (a.e > e+1) {
         z = a;
         return;
      }

      z.x = a.x + x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return;
   }
   else {
      if (e > a.e+1) {
         z.x = x;
         z.e = e;
         z.normalize();
         return;
      }

      z.x = a.x*NTL_XD_BOUND_INV + x;
      z.e = e;
      z.normalize();
      return;
   }
}

void MulSub(xdouble& z, const xdouble& a, const xdouble& b, const xdouble& c)
// z = a - b*c
{
   double x;
   long e;

   e = b.e + c.e;
   x = b.x * c.x;

   if (x == 0) { 
      z = a;
      return;
   }

   if (a.x == 0) {
      z.e = e;
      z.x = -x;
      z.normalize();
      return;
   }
      

   if (a.e == e) {
      z.x = a.x - x;
      z.e = e;
      z.normalize();
      return;
   }
   else if (a.e > e) {
      if (a.e > e+1) {
         z = a;
         return;
      }

      z.x = a.x - x*NTL_XD_BOUND_INV;
      z.e = a.e;
      z.normalize();
      return;
   }
   else {
      if (e > a.e+1) {
         z.x = -x;
         z.e = e;
         z.normalize();
         return;
      }

      z.x = a.x*NTL_XD_BOUND_INV - x;
      z.e = e;
      z.normalize();
      return;
   }
}

double log(const xdouble& a)
{
   static double LogBound = log(NTL_XD_BOUND);
   if (a.x <= 0) {
      Error("log(xdouble): argument must be positive");
   }

   return log(a.x) + a.e*LogBound;
}

xdouble xexp(double x)
{
   const double LogBound = log(NTL_XD_BOUND);

   double y = x/LogBound;
   double iy = floor(y+0.5);

   if (iy >= NTL_OVFBND)
      Error("xdouble: overflow");

   if (iy <= -NTL_OVFBND)
      Error("xdouble: underflow");


   double fy = y - iy;

   xdouble res;
   res.e = long(iy);
   res.x = exp(fy*LogBound);
   res.normalize();
   return res;
}

/**************  input / output routines **************/


void ComputeLn2(RR&);
void ComputeLn10(RR&);

long ComputeMax10Power()
{
   long old_p = RR::precision();
   RR::SetPrecision(NTL_BITS_PER_LONG);

   RR ln2, ln10;
   ComputeLn2(ln2);
   ComputeLn10(ln10);

   long k = to_long( to_RR(NTL_OVFBND/2) * ln2 / ln10 );

   RR::SetPrecision(old_p);
   return k;
}


xdouble PowerOf10(const ZZ& e)
{
   static long init = 0;
   static xdouble v10k;
   static long k;

   if (!init) {
      long old_p = RR::precision();
      k = ComputeMax10Power();
      RR::SetPrecision(NTL_DOUBLE_PRECISION);
      v10k = to_xdouble(power(to_RR(10), k)); 
      RR::SetPrecision(old_p);
      init = 1;
   }

   ZZ e1;
   long neg;

   if (e < 0) {
      e1 = -e;
      neg = 1;
   }
   else {
      e1 = e;
      neg = 0;
   }

   long r;
   ZZ q;

   r = DivRem(q, e1, k);

   long old_p = RR::precision();
   RR::SetPrecision(NTL_DOUBLE_PRECISION);
   xdouble x1 = to_xdouble(power(to_RR(10), r));
   RR::SetPrecision(old_p);

   xdouble x2 = power(v10k, q);
   xdouble x3 = x1*x2;

   if (neg) x3 = 1/x3;

   return x3;
}




ostream& operator<<(ostream& s, const xdouble& a)
{
   if (a == 0) {
      s << "0";
      return s;
   }

   long old_p = RR::precision();
   long temp_p = long(log(fabs(log(fabs(a))) + 1.0)/log(2.0)) + 10; 

   RR::SetPrecision(temp_p);

   RR ln2, ln10, log_2_10;
   ComputeLn2(ln2);
   ComputeLn10(ln10);
   log_2_10 = ln10/ln2;
   ZZ log_10_a = to_ZZ(
  (to_RR(a.e)*to_RR(2*NTL_XD_HBOUND_LOG) + log(fabs(a.x))/log(2.0))/log_2_10);

   RR::SetPrecision(old_p);

   xdouble b;
   long neg;

   if (a < 0) {
      b = -a;
      neg = 1;
   }
   else {
      b = a;
      neg = 0;
   }

   ZZ k = xdouble::OutputPrecision() - log_10_a;

   xdouble c, d;

   c = PowerOf10(to_ZZ(xdouble::OutputPrecision()));
   d = PowerOf10(log_10_a);

   b = b / d;
   b = b * c;

   while (b < c) {
      b = b * 10.0;
      k++;
   }

   while (b >= c) {
      b = b / 10.0;
      k--;
   }

   b = b + 0.5;
   k = -k;

   ZZ B;
   conv(B, b);

   long bp_len = xdouble::OutputPrecision()+10;

   char *bp = NTL_NEW_OP char[bp_len];

   if (!bp) Error("xdouble output: out of memory");

   long len, i;

   len = 0;
   do {
      if (len >= bp_len) Error("xdouble output: buffer overflow");
      bp[len] = IntValToChar(DivRem(B, B, 10));
      len++;
   } while (B > 0);

   for (i = 0; i < len/2; i++) {
      char tmp;
      tmp = bp[i];
      bp[i] = bp[len-1-i];
      bp[len-1-i] = tmp;
   }

   i = len-1;
   while (bp[i] == '0') i--;

   k += (len-1-i);
   len = i+1;

   bp[len] = '\0';

   if (k > 3 || k < -len - 3) {
      // use scientific notation

      if (neg) s << "-";
      s << "0." << bp << "e" << (k + len);
   }
   else {
      long kk = to_long(k);

      if (kk >= 0) {
         if (neg) s << "-";
         s << bp;
         for (i = 0; i < kk; i++) 
            s << "0";
      }
      else if (kk <= -len) {
         if (neg) s << "-";
         s << "0.";
         for (i = 0; i < -len-kk; i++)
            s << "0";
         s << bp;
      }
      else {
         if (neg) s << "-";
         for (i = 0; i < len+kk; i++)
            s << bp[i];
   
         s << ".";
   
         for (i = len+kk; i < len; i++)
            s << bp[i];
      }
   }

   delete [] bp;
   return s;
}

istream& operator>>(istream& s, xdouble& x)
{
   long c;
   long cval;
   long sign;
   ZZ a, b;

   if (!s) Error("bad xdouble input");

   c = s.peek();
   while (IsWhiteSpace(c)) {
      s.get();
      c = s.peek();
   }

   if (c == '-') {
      sign = -1;
      s.get();
      c = s.peek();
   }
   else
      sign = 1;

   long got1 = 0;
   long got_dot = 0;
   long got2 = 0;

   a = 0;
   b = 1;

   cval = CharToIntVal(c);

   if (cval >= 0 && cval <= 9) {
      got1 = 1;

      while (cval >= 0 && cval <= 9) {
         mul(a, a, 10);
         add(a, a, cval);
         s.get();
         c = s.peek();
         cval = CharToIntVal(c);
      }
   }

   if (c == '.') {
      got_dot = 1;

      s.get();
      c = s.peek();
      cval = CharToIntVal(c);

      if (cval >= 0 && cval <= 9) {
         got2 = 1;
   
         while (cval >= 0 && cval <= 9) {
            mul(a, a, 10);
            add(a, a, cval);
            mul(b, b, 10);
            s.get();
            c = s.peek();
            cval = CharToIntVal(c);
         }
      }
   }

   if (got_dot && !got1 && !got2)  Error("bad xdouble input");

   ZZ e;

   long got_e = 0;
   long e_sign;

   if (c == 'e' || c == 'E') {
      got_e = 1;

      s.get();
      c = s.peek();

      if (c == '-') {
         e_sign = -1;
         s.get();
         c = s.peek();
      }
      else if (c == '+') {
         e_sign = 1;
         s.get();
         c = s.peek();
      }
      else
         e_sign = 1;

      cval = CharToIntVal(c);

      if (cval < 0 || cval > 9) Error("bad xdouble input");

      e = 0;
      while (cval >= 0 && cval <= 9) {
         mul(e, e, 10);
         add(e, e, cval);
         s.get();
         c = s.peek();
         cval = CharToIntVal(c);
      }
   }

   if (!got1 && !got2 && !got_e) Error("bad xdouble input");

   xdouble t1, t2, v;

   if (got1 || got2) {
      conv(t1, a);
      conv(t2, b);
      v = t1/t2;
   }
   else
      v = 1;

   if (sign < 0)
      v = -v;

   if (got_e) {
      if (e_sign < 0) negate(e, e);
      t1 = PowerOf10(e);
      v = v * t1;
   }

   x = v;
   return s;
}


xdouble to_xdouble(const char *s)
{
   long c;
   long cval;
   long sign;
   ZZ a, b;
   long i=0;

   if (!s) Error("bad xdouble input");

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

   long got1 = 0;
   long got_dot = 0;
   long got2 = 0;

   a = 0;
   b = 1;

   cval = CharToIntVal(c);

   if (cval >= 0 && cval <= 9) {
      got1 = 1;

      while (cval >= 0 && cval <= 9) {
         mul(a, a, 10);
         add(a, a, cval);
         i++;
         c = s[i];
         cval = CharToIntVal(c);
      }
   }

   if (c == '.') {
      got_dot = 1;

      i++;
      c = s[i];
      cval = CharToIntVal(c);

      if (cval >= 0 && cval <= 9) {
         got2 = 1;
   
         while (cval >= 0 && cval <= 9) {
            mul(a, a, 10);
            add(a, a, cval);
            mul(b, b, 10);
            i++;
            c = s[i];
            cval = CharToIntVal(c);
         }
      }
   }

   if (got_dot && !got1 && !got2)  Error("bad xdouble input");

   ZZ e;

   long got_e = 0;
   long e_sign;

   if (c == 'e' || c == 'E') {
      got_e = 1;

      i++;
      c = s[i];

      if (c == '-') {
         e_sign = -1;
         i++;
         c = s[i];
      }
      else if (c == '+') {
         e_sign = 1;
         i++;
         c = s[i];
      }
      else
         e_sign = 1;

      cval = CharToIntVal(c);

      if (cval < 0 || cval > 9) Error("bad xdouble input");

      e = 0;
      while (cval >= 0 && cval <= 9) {
         mul(e, e, 10);
         add(e, e, cval);
         i++;
         c = s[i];
         cval = CharToIntVal(c);
      }
   }

   if (!got1 && !got2 && !got_e) Error("bad xdouble input");

   xdouble t1, t2, v;

   if (got1 || got2) {
      conv(t1, a);
      conv(t2, b);
      v = t1/t2;
   }
   else
      v = 1;

   if (sign < 0)
      v = -v;

   if (got_e) {
      if (e_sign < 0) negate(e, e);
      t1 = PowerOf10(e);
      v = v * t1;
   }

   return v;
}

NTL_END_IMPL
