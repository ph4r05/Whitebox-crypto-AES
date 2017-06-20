/*
Copyright (C) 1997, 1998, 1999, 2000 Victor Shoup

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*****************************************************

The quad_float package is derived from the doubledouble package of
Keith Briggs.  However, the version employed in NTL has been extensively 
modified.  Below, I attach the copyright notice from the original
doubledouble package, which is currently available at 

   http://www.labs.bt.com/people/briggsk2

*****************************************************

Copyright (C) 1997 Keith Martin Briggs

Except where otherwise indicated,
this program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <NTL/quad_float.h>
#include <NTL/RR.h>

#include <float.h>

#include <NTL/new.h>

NTL_START_IMPL

#if (defined(__GNUC__) && NTL_EXT_DOUBLE && (defined(__i386__) || defined(__i486__) || defined(__i586__)))

#if (!defined(NTL_X86_FIX) && !defined(NTL_NO_X86_FIX))

#define NTL_X86_FIX

#endif

#endif


#if (NTL_EXT_DOUBLE && !defined(NTL_X86_FIX))

#define DOUBLE volatile double

#else

#define DOUBLE double

#endif


#ifdef NTL_X86_FIX


#define START_FIX \
  volatile unsigned short __old_cw, __new_cw; \
  asm volatile ("fnstcw %0":"=m" (__old_cw)); \
  __new_cw = (__old_cw & ~0x300) | 0x200; \
  asm volatile ("fldcw %0": :"m" (__new_cw));


#define END_FIX  asm volatile ("fldcw %0": :"m" (__old_cw));

#else

#define START_FIX
#define END_FIX

#endif


static 
void normalize(quad_float& z, const double& xhi, const double& xlo)
{
START_FIX
   DOUBLE u, v;

   u = xhi + xlo; 
   v = xhi - u;    
   v = v + xlo;    

   z.hi = u;
   z.lo = v;
END_FIX
}



#if (NTL_BITS_PER_LONG >= NTL_DOUBLE_PRECISION)
quad_float to_quad_float(long n)
{
START_FIX
   DOUBLE xhi, xlo;
   DOUBLE u, v;

   xhi = double(n);

   // Because we are assuming 2's compliment integer
   // arithmetic, the following prevents long(xhi) from overflowing.

   if (n > 0)
      xlo = double(n+long(-xhi));
   else
      xlo = double(n-long(xhi));

   // renormalize...just to be safe
  
   u = xhi + xlo;
   v = xhi - u;
   v = v + xlo;
END_FIX
   return quad_float(u, v);
}

quad_float to_quad_float(unsigned long n)
{
START_FIX
   DOUBLE xhi, xlo, t;
   DOUBLE u, v;

   const double bnd = double(1L << (NTL_BITS_PER_LONG-2))*4.0;

   xhi = double(n);
   
   if (xhi >= bnd)
      t = xhi - bnd;
   else
      t = xhi;

   // we use the "to_long" function here to be as portable as possible.
   long llo = to_long(n - (unsigned long)(t));
   xlo = double(llo);

   // renormalize...just to be safe

   u = xhi + xlo;
   v = xhi - u;
   v = v + xlo;
END_FIX
   return quad_float(u, v);
}
#endif


long quad_float::oprec = 10;

void quad_float::SetOutputPrecision(long p)
{
   if (p < 1) p = 1;

   if (NTL_OVERFLOW(p, 1, 0)) 
      Error("quad_float: output precision too big");

   oprec = p;
}


quad_float operator +(const quad_float& x, const quad_float& y ) {
START_FIX
        DOUBLE    H, h, T, t, S, s, e, f;
        DOUBLE    t1;

        S = x.hi + y.hi;
        T = x.lo + y.lo;
        e = S - x.hi;
        f = T - x.lo;

        t1 = S-e;
        t1 = x.hi-t1;
        s = y.hi-e;
        s = s + t1;
        
        t1 = T-f;
        t1 = x.lo-t1;
        t = y.lo-f;
        t = t + t1;


        s = s + T;
        H = S + s;
        h = S - H;
        h = h + s;

        h = h + t;
        e = H + h; 
        f = H - e;
        f = f + h;
END_FIX
        return quad_float(e, f);
}

quad_float& operator +=(quad_float& x, const quad_float& y ) {
START_FIX
        DOUBLE    H, h, T, t, S, s, e, f;
        DOUBLE    t1;

        S = x.hi + y.hi;
        T = x.lo + y.lo;
        e = S - x.hi;
        f = T - x.lo;

        t1 = S-e;
        t1 = x.hi-t1;
        s = y.hi-e;
        s = s + t1;
        
        t1 = T-f;
        t1 = x.lo-t1;
        t = y.lo-f;
        t = t + t1;


        s = s + T;
        H = S + s;
        h = S - H;
        h = h + s;

        h = h + t;
        e = H + h; 
        f = H - e;
        f = f + h;

        x.hi = e;
        x.lo = f;
END_FIX
        return x;
}

quad_float operator -(const quad_float& x, const quad_float& y ) {
START_FIX
        DOUBLE    H, h, T, t, S, s, e, f;
        DOUBLE    t1, yhi, ylo;

        yhi = -y.hi;
        ylo = -y.lo;

        S = x.hi + yhi;
        T = x.lo + ylo;
        e = S - x.hi;
        f = T - x.lo;

        t1 = S-e;
        t1 = x.hi-t1;
        s = yhi-e;
        s = s + t1;
        
        t1 = T-f;
        t1 = x.lo-t1;
        t = ylo-f;
        t = t + t1;


        s = s + T;
        H = S + s;
        h = S - H;
        h = h + s;

        h = h + t;
        e = H + h; 
        f = H - e;
        f = f + h;

END_FIX
        return quad_float(e, f);
}

quad_float& operator -=(quad_float& x, const quad_float& y ) {
START_FIX
	DOUBLE    H, h, T, t, S, s, e, f;
        DOUBLE    t1, yhi, ylo;

        yhi = -y.hi;
        ylo = -y.lo;

        S = x.hi + yhi;
        T = x.lo + ylo;
        e = S - x.hi;
        f = T - x.lo;

        t1 = S-e;
        t1 = x.hi-t1;
        s = yhi-e;
        s = s + t1;
        
        t1 = T-f;
        t1 = x.lo-t1;
        t = ylo-f;
        t = t + t1;


        s = s + T;
        H = S + s;
        h = S - H;
        h = h + s;

        h = h + t;
        e = H + h; 
        f = H - e;
        f = f + h;

        x.hi = e;
        x.lo = f;
END_FIX
        return x;
}

quad_float operator -(const quad_float& x)
{
START_FIX
   DOUBLE xhi, xlo, u, v;

   xhi = -x.hi;
   xlo = -x.lo;

   // it is a good idea to renormalize here, just in case
   // the rounding rule depends on sign, and thus we will
   // maintain the "normal form" for quad_float's.
  
   u = xhi + xlo;
   v = xhi - u;
   v = v + xlo;

END_FIX
   return quad_float(u, v);
}


quad_float operator *(const quad_float& x,const quad_float& y ) {
START_FIX
  DOUBLE hx, tx, hy, ty, C, c;
  DOUBLE t1, t2;

  C = NTL_QUAD_FLOAT_SPLIT*x.hi;
  hx = C-x.hi;
  c = NTL_QUAD_FLOAT_SPLIT*y.hi;
  hx = C-hx;
  tx = x.hi-hx;
  hy = c-y.hi;
  C = x.hi*y.hi;
  hy = c-hy;
  ty = y.hi-hy;

  // c = ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(x.hi*y.lo+x.lo*y.hi);
  
  t1 = hx*hy;
  t1 = t1-C;
  t2 = hx*ty;
  t1 = t1+t2;
  t2 = tx*hy;
  t1 = t1+t2;
  t2 = tx*ty;
  c = t1+t2;
  t1 = x.hi*y.lo;
  t2 = x.lo*y.hi;
  t1 = t1+t2;
  c = c + t1;


  hx = C+c;
  tx = C-hx;
  tx = tx+c;

END_FIX
  return quad_float(hx, tx);
}

quad_float& operator *=(quad_float& x,const quad_float& y ) {
START_FIX
  DOUBLE hx, tx, hy, ty, C, c;
  DOUBLE t1, t2;

  C = NTL_QUAD_FLOAT_SPLIT*x.hi;
  hx = C-x.hi;
  c = NTL_QUAD_FLOAT_SPLIT*y.hi;
  hx = C-hx;
  tx = x.hi-hx;
  hy = c-y.hi;
  C = x.hi*y.hi;
  hy = c-hy;
  ty = y.hi-hy;

  // c = ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(x.hi*y.lo+x.lo*y.hi);
  
  t1 = hx*hy;
  t1 = t1-C;
  t2 = hx*ty;
  t1 = t1+t2;
  t2 = tx*hy;
  t1 = t1+t2;
  t2 = tx*ty;
  c = t1+t2;
  t1 = x.hi*y.lo;
  t2 = x.lo*y.hi;
  t1 = t1+t2;
  c = c + t1;


  hx = C+c;
  tx = C-hx;
  tx = tx+c;

  x.hi = hx;
  x.lo = tx;
END_FIX
  return x;
}

quad_float operator /(const quad_float& x, const quad_float& y ) {
START_FIX
  DOUBLE hc, tc, hy, ty, C, c, U, u;
  DOUBLE t1;

  C = x.hi/y.hi;
  c = NTL_QUAD_FLOAT_SPLIT*C;
  hc = c-C;
  u = NTL_QUAD_FLOAT_SPLIT*y.hi;
  hc = c-hc;
  tc = C-hc;
  hy = u-y.hi;
  U = C * y.hi;
  hy = u-hy;
  ty = y.hi-hy;

  // u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;

  u = hc*hy;
  u = u-U;
  t1 = hc*ty;
  u = u+t1;
  t1 = tc*hy;
  u = u+t1;
  t1 = tc*ty;
  u = u+t1;

  // c = ((((x.hi-U)-u)+x.lo)-C*y.lo)/y.hi;

  c = x.hi-U;
  c = c-u;
  c = c+x.lo;
  t1 = C*y.lo;
  c = c - t1;
  c = c/y.hi;
  
  hy = C+c;
  ty = C-hy;
  ty = ty+c;

END_FIX
  return quad_float(hy, ty);
}

quad_float& operator /=(quad_float& x, const quad_float& y ) {
START_FIX
  DOUBLE hc, tc, hy, ty, C, c, U, u;
  DOUBLE t1;

  C = x.hi/y.hi;
  c = NTL_QUAD_FLOAT_SPLIT*C;
  hc = c-C;
  u = NTL_QUAD_FLOAT_SPLIT*y.hi;
  hc = c-hc;
  tc = C-hc;
  hy = u-y.hi;
  U = C * y.hi;
  hy = u-hy;
  ty = y.hi-hy;

  // u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;

  u = hc*hy;
  u = u-U;
  t1 = hc*ty;
  u = u+t1;
  t1 = tc*hy;
  u = u+t1;
  t1 = tc*ty;
  u = u+t1;

  // c = ((((x.hi-U)-u)+x.lo)-C*y.lo)/y.hi;

  c = x.hi-U;
  c = c-u;
  c = c+x.lo;
  t1 = C*y.lo;
  c = c - t1;
  c = c/y.hi;
  
  hy = C+c;
  ty = C-hy;
  ty = ty+c;

  x.hi = hy;
  x.lo = ty;
END_FIX
  return x;
}


quad_float sqrt(const quad_float& y) {
  if (y.hi < 0.0) 
    Error("Quad: attempto to take square root of negative number");
  if (y.hi == 0.0) return quad_float(0.0,0.0);

  double c;
  c = sqrt(y.hi);
  ForceToMem(&c);  // This is fairly paranoid, but it doesn't cost too much.

START_FIX

  DOUBLE p,q,hx,tx,u,uu,cc;
  DOUBLE t1;

  p = NTL_QUAD_FLOAT_SPLIT*c; 
  hx = (c-p); 
  hx = hx+p; 
  tx = c-hx;
  p = hx*hx;
  q = hx*tx;
  q = q+q;

  u = p+q;
  uu = p-u;
  uu = uu+q;
  t1 = tx*tx;
  uu = uu+t1;


  cc = y.hi-u;
  cc = cc-uu;
  cc = cc+y.lo;
  t1 = c+c;
  cc = cc/t1;

  hx = c+cc;
  tx = c-hx;
  tx = tx+cc;
END_FIX
  return quad_float(hx, tx);
}



void power(quad_float& z, const quad_float& a, long e)
{
   quad_float res, u;
   unsigned long k;

   if (e < 0)
      k = -((unsigned long) e);
   else
      k = e;

   res = 1.0;
   u = a;

   while (k) {
      if (k & 1)
         res = res * u;

      k = k >> 1;
      if (k)
         u = u * u;
   }

   if (e < 0)
      z = 1.0/res;
   else
      z = res;
}


void power2(quad_float& z, long e)
{
   z.hi = _ntl_ldexp(1.0, e);
   z.lo = 0;
}


long to_long(const quad_float& x)
{
   double fhi, flo;

   fhi = floor(x.hi);

   if (fhi == x.hi) 
      flo = floor(x.lo);
   else
      flo = 0;

   // the following code helps to prevent unnecessary integer overflow,
   // and guarantees that to_long(to_quad_float(a)) == a, for all long a,
   // provided long's are not too wide.

   if (fhi > 0)
      return long(flo) - long(-fhi);
   else
      return long(fhi) + long(flo);
}



// This version of ZZ to quad_float coversion relies on the
// precise rounding rules implemented by the ZZ to double conversion.


void conv(quad_float& z, const ZZ& a)
{
   double xhi, xlo;

   conv(xhi, a);

   if (!IsFinite(&xhi)) {
      z.hi = xhi;
      z.lo = 0;
      return;
   }

   static ZZ t;

   conv(t, xhi);
   sub(t, a, t);

   conv(xlo, t);

   normalize(z, xhi, xlo);

   // The following is just paranoia.
   if (fabs(z.hi) < NTL_FDOUBLE_PRECISION && z.lo != 0)
      Error("internal error: ZZ to quad_float conversion");
} 

void conv(ZZ& z, const quad_float& x)
{ 
   static ZZ t1, t2, t3;

   double fhi, flo;

   fhi = floor(x.hi);

   if (fhi == x.hi) {
      flo = floor(x.lo);

      conv(t1, fhi);
      conv(t2, flo);

      add(z, t1, t2);
   }
   else
      conv(z, fhi);
}



ostream& operator<<(ostream& s, const quad_float& a)
{
   quad_float aa = a;

   if (!IsFinite(&aa)) {
      s << "NaN";
      return s;
   }

   long old_p = RR::precision();
   long old_op = RR::OutputPrecision();

   RR::SetPrecision(long(3.33*quad_float::oprec) + 10);
   RR::SetOutputPrecision(quad_float::oprec);

   static RR t;

   conv(t, a);
   s << t;

   RR::SetPrecision(old_p);
   RR::SetOutputPrecision(old_op);

   return s;
}

istream& operator>>(istream& s, quad_float& x)
{
   long old_p = RR::precision();
   RR::SetPrecision(4*NTL_DOUBLE_PRECISION);

   static RR t;
   s >> t;
   conv(x, t);

   RR::SetPrecision(old_p);
   return s;
}

void random(quad_float& x)
{
   long old_p = RR::precision();
   RR::SetPrecision(4*NTL_DOUBLE_PRECISION);

   static RR t;
   random(t);
   conv(x, t);

   RR::SetPrecision(old_p);
}

quad_float random_quad_float()
{
   quad_float x;
   random(x);
   return x;
}
      
long IsFinite(quad_float *x)
{
   return IsFinite(&x->hi) && IsFinite(&x->lo);
}


long PrecisionOK()
{
START_FIX
   long k;
   DOUBLE l1 = (double)1;
   DOUBLE lh = 1/(double)2;
   DOUBLE epsilon;
   DOUBLE fudge, oldfudge;

   epsilon = l1;
   fudge = l1+l1;

   k = 0;

   do {
      k++;
      epsilon = epsilon * lh;
      oldfudge = fudge;
      fudge = l1 + epsilon;
   } while (fudge > l1 && fudge < oldfudge);

END_FIX
   return k == NTL_DOUBLE_PRECISION;
}

quad_float floor(const quad_float& x)
{
   double fhi = floor(x.hi);

   if (fhi != x.hi)
      return quad_float(fhi, 0.0);
   else {
      double flo = floor(x.lo);
      quad_float z;
      normalize(z, fhi, flo);
      return z;
   }
}


quad_float ceil(const quad_float& x) { 
  return -floor(-x);
}

quad_float trunc(const quad_float& x) { 
  if (x>=0.0) return floor(x); else return -floor(-x);
}



long compare(const quad_float& x, const quad_float& y)
{
   if (x.hi > y.hi) 
      return 1;
   else if (x.hi < y.hi)
      return -1;
   else if (x.lo > y.lo)
      return 1;
   else if (x.lo < y.lo) 
      return -1;
   else
      return 0;
}


quad_float fabs(const quad_float& x) 
{ if (x.hi>=0.0) return x; else return -x; }

quad_float to_quad_float(const char *s)
{
   quad_float x;
   long old_p = RR::precision();
   RR::SetPrecision(4*NTL_DOUBLE_PRECISION);

   static RR t;
   conv(t, s);
   conv(x, t);

   RR::SetPrecision(old_p);
   return x;
}


quad_float ldexp(const quad_float& x, long exp) { // x*2^exp
   double xhi, xlo;
   quad_float z;

   xhi = _ntl_ldexp(x.hi, exp);
   xlo = _ntl_ldexp(x.lo, exp);

   normalize(z, xhi, xlo);
   return z;
}


quad_float exp(const quad_float& x) { // New version 97 Aug 05
/*
!  Calculate a quadruple-precision exponential
!  Method:
!   x    x.log2(e)    nint[x.log2(e)] + frac[x.log2(e)]
!  e  = 2          = 2
!
!                     iy    fy
!                  = 2   . 2
!  Then
!   fy    y.loge(2)
!  2   = e
!
!  Now y.loge(2) will be less than 0.3466 in absolute value.
!  This is halved and a Pade aproximation is used to approximate e^x over
!  the region (-0.1733, +0.1733).   This approximation is then squared.
*/
  if (x.hi<DBL_MIN_10_EXP*2.302585092994045684017991) 
    return to_quad_float(0.0);
  if (x.hi>DBL_MAX_10_EXP*2.302585092994045684017991) {
    Error("exp(quad_float): overflow");
  }

  // changed this from "const" to "static" in v5.3, since "const"
  // causes the initialization to be performed with *every* invocation.
  static quad_float Log2 = 
    to_quad_float("0.6931471805599453094172321214581765680755");

  quad_float y,temp,ysq,sum1,sum2;
  long iy;
  y=x/Log2;
  temp = floor(y+0.5);
  iy = to_long(temp);
  y=(y-temp)*Log2;
  y=ldexp(y,-1L);
  ysq=y*y;
  sum1=y*((((ysq+3960.0)*ysq+2162160.0)*ysq+302702400.0)*ysq+8821612800.0);
  sum2=(((90.0*ysq+110880.0)*ysq+30270240.0)*ysq+2075673600.0)*ysq+17643225600.0;
/*
!                     sum2 + sum1         2.sum1
! Now approximation = ----------- = 1 + ----------- = 1 + 2.temp
!                     sum2 - sum1       sum2 - sum1
!
! Then (1 + 2.temp)^2 = 4.temp.(1 + temp) + 1
*/
  temp=sum1/(sum2-sum1);
  y=temp*(temp+1);
  y=ldexp(y,2L);
  return ldexp(y+1,iy);
}

quad_float log(const quad_float& t) { // Newton method. See Bailey, MPFUN
  if (t.hi <= 0.0) {
    Error("log(quad_float): argument must be positive");
  }
  double s1 = log(t.hi);
  ForceToMem(&s1);  // Again, this is fairly paranoid.
  quad_float s;
  s = s1;
  quad_float e;
  e=exp(s);
  return s+(t-e)/e;  // Newton step
}

long operator> (const quad_float& x, const quad_float& y) {
   return (x.hi> y.hi) || (x.hi==y.hi && x.lo> y.lo); }
long operator>=(const quad_float& x, const quad_float& y) {
   return (x.hi>y.hi) || (x.hi==y.hi && x.lo>=y.lo); }
long operator< (const quad_float& x, const quad_float& y) {
   return (x.hi< y.hi) || (x.hi==y.hi && x.lo< y.lo); }
long operator<=(const quad_float& x, const quad_float& y) {
   return (x.hi<y.hi) || (x.hi==y.hi && x.lo<=y.lo); }
long operator==(const quad_float& x, const quad_float& y)
   { return x.hi==y.hi && x.lo==y.lo; }
long operator!=(const quad_float& x, const quad_float& y)
   { return x.hi!=y.hi || x.lo!=y.lo; }


NTL_END_IMPL
