
#include <NTL/lip.h>

#include <NTL/ctools.h>


#include <stdlib.h>
#include <stdio.h>
#include <math.h>






#define MustAlloc(c, len)  (!(c) || ((c)[-1] >> 1) < (len))
   /* Fast test to determine if allocation is necessary */


static 
void zhalt(char *c)
{
   fprintf(stderr,"fatal error:\n   %s\nexit...\n",c);
   fflush(stderr);
   _ntl_abort();
}



#define NTL_GMP_MUL_CROSS (4)
#define NTL_GMP_SQR_CROSS (4)
#define NTL_GMP_DIV_CROSS (4)
#define NTL_GMP_REM_CROSS (4)
#define NTL_GMP_QREM_CROSS (4)
#define NTL_GMP_SQRT_CROSS (2)
#define NTL_GMP_GCD_CROSS (1)
#define NTL_GMP_XGCD_CROSS (4)
#define NTL_GMP_INVMOD_CROSS (2)

#ifdef NTL_GMP_HACK

/* using GMP */

#ifdef NTL_SINGLE_MUL
#error "sorry...not implemented"
#endif

#if (NTL_NBITS != NTL_NBITS_MAX)
#error "sorry...not implemented"
#endif

int _ntl_gmp_hack = 1;


#include <gmp.h>


static mp_limb_t *gspace_data = 0;
static long gspace_size = 0;


static void alloc_gspace(long n)
{
   if (n <= gspace_size) return;
   
   if (n <= gspace_size*1.1)
      n = ((long) (gspace_size*1.1)) + 10;

   if (gspace_data)
      gspace_data = 
         (mp_limb_t *) NTL_REALLOC(gspace_data, n, sizeof(mp_limb_t), 0);   
   else
      gspace_data = (mp_limb_t *) NTL_MALLOC(n, sizeof(mp_limb_t), 0);

   if (!gspace_data) zhalt("alloc_gspace: out of memory");

   gspace_size = n;
}


#define NTL_GSPACE(n) \
{ long _n = (n);  if (_n > gspace_size) alloc_gspace(_n); } 


static mpz_t gt_1, gt_2, gt_3, gt_4, gt_5;
static long gt_init = 0;

static
void do_gt_init()
{
   gt_init = 1;
   mpz_init(gt_1);
   mpz_init(gt_2);
   mpz_init(gt_3);
   mpz_init(gt_4);
   mpz_init(gt_5);
}

#define GT_INIT { if (!gt_init) do_gt_init(); }

#include "lip_gmp_aux_impl.h"


/* Check that we really have it...otherwise, some compilers
 * only issue warnings for missing functions.
 */


#ifndef HAVE_LIP_GMP_AUX

#error "do not have lip_gmp_aux.h"

#endif

static
void lip_to_mpz(const long *x, mpz_t y)
{
   long sx;
   long neg;
   long gsx;

   if (!x || (x[0] == 1 && x[1] == 0)) {
      mpz_set_ui(y, 0);
      return;
   }

   sx = x[0];
   if (sx < 0) {
      neg = 1;
      sx = -sx;
   }
   else
      neg = 0;

   gsx = L_TO_G(sx);

   if (gsx > y->_mp_alloc) _mpz_realloc(y, gsx);
   lip_to_gmp(x+1, y->_mp_d, sx);
   while (y->_mp_d[gsx-1] == 0) gsx--;
   if (neg) gsx = -gsx;
   y->_mp_size = gsx;
}

static
void mpz_to_lip(long **x, mpz_t y)
{
   long gsy;
   long neg;
   long sx;
   long *xp;

   gsy = y->_mp_size;

   if (gsy == 0) {
      _ntl_zzero(x);
      return;
   }

   if (gsy < 0) {
      neg = 1;
      gsy = -gsy;
   }
   else
      neg = 0;

   sx = G_TO_L(gsy);

   xp = *x;
   if (MustAlloc(xp, sx)) { 
      _ntl_zsetlength(x, sx); 
      xp = *x;
   }

   gmp_to_lip1(xp+1, y->_mp_d, sx);

   while (xp[sx] == 0) sx--;
   if (neg) sx = -sx;
   xp[0] = sx; 
}


void _ntl_gmp_powermod(long **x, long *a, long *e, long *n)
{
   long neg = 0;

   GT_INIT


   lip_to_mpz(a, gt_2);
   lip_to_mpz(e, gt_3);
   lip_to_mpz(n, gt_4);

   if (mpz_sgn(gt_3) < 0) {
      mpz_neg(gt_3, gt_3);
      neg = 1;
   }

   mpz_powm(gt_1, gt_2, gt_3, gt_4);

   /* This fixes a bug in GMP 3.0.1 */
   if (mpz_cmp(gt_1, gt_4) >= 0) 
      mpz_sub(gt_1, gt_1, gt_4);

   if (neg) {
      if (!mpz_invert(gt_1, gt_1, gt_4))
         zhalt("Modular inverse not defined");
   }

   mpz_to_lip(x, gt_1);
}


#else

/* not using GMP */

int _ntl_gmp_hack = 0;


#endif



#define DENORMALIZE  NTL_FDOUBLE_PRECISION


#define MIN_SETL	(4)
   /* _ntl_zsetlength allocates a multiple of MIN_SETL digits */



#if (defined(NTL_SINGLE_MUL) && !defined(NTL_FAST_INT_MUL))
#define MulLo(rres,a,b) \
{ \
   double _x = ((double) a) * ((double) b) + (DENORMALIZE); \
   unsigned long _x_lo; \
   NTL_FetchLo(_x_lo, _x); \
   rres =  _x_lo; \
}
#else
#define MulLo(rres,a,b) rres = (a)*(b)
#endif






/*
 * definitions of zaddmulp, zxmulp, zaddmulpsq for the various
 * long integer arithmentic implementation options.
 */


#if (defined(NTL_SINGLE_MUL))

#define zaddmulp(a,b,d,t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   unsigned long  __lhi, __llo;\
   double __lx = ((double) ((_a) + (_t))) + ((double) _b)*((double) _d); \
  __lx += (DENORMALIZE); \
   NTL_FetchHiLo(__lhi,__llo,__lx);\
   __lhi = ((__lhi<<6)|(__llo>>26)) & 0x3FFFFFF; \
   __llo &= 0x3FFFFFF; \
   (a) = __llo;\
   (t) = __lhi;\
}

#define zxmulp(a,b,d,t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   unsigned long  __lhi, __llo;\
   double __lx = ((double) ((_t))) + ((double) _b)*((double) _d); \
  __lx += (DENORMALIZE); \
   NTL_FetchHiLo(__lhi,__llo,__lx);\
   __lhi = ((__lhi<<6)|(__llo>>26)) & 0x3FFFFFF; \
   __llo &= 0x3FFFFFF; \
   (a) = __llo;\
   (t) = __lhi;\
}

#define zaddmulpsq(a,b,t) \
{ \
   long _a = (a), _b = (b); \
   unsigned long  __lhi, __llo; \
   double __lx = ((double) (_a)) + ((double) _b)*((double) _b); \
   __lx += (DENORMALIZE); \
   NTL_FetchHiLo(__lhi,__llo,__lx);\
   __lhi = ((__lhi<<6)|(__llo>>26)) & 0x3FFFFFF; \
   __llo &= 0x3FFFFFF; \
   (a) =  __llo;\
   (t) =  __lhi;\
}

#elif (defined(NTL_LONG_LONG))


#if (!defined(NTL_CLEAN_INT))

/*
 * One might get slightly better code with this version.
 */

#define zaddmulp(a, b, d, t) { \
   NTL_LL_TYPE _pp = ((NTL_LL_TYPE) (b)) * ((NTL_LL_TYPE) (d)) + ((t)+(a)); \
   (a) = ((long)(_pp)) & NTL_RADIXM; \
   (t) = (long) (_pp >> NTL_NBITS); \
} 


#define zxmulp(a, b, d, t) { \
   NTL_LL_TYPE _pp = ((NTL_LL_TYPE) (b)) * ((NTL_LL_TYPE) (d)) + (t); \
   (a) = ((long)(_pp)) & NTL_RADIXM; \
   (t) = (long) (_pp >> NTL_NBITS); \
} 

#define zaddmulpsq(a,b,t) { \
   NTL_LL_TYPE _pp = ((NTL_LL_TYPE) (b)) * ((NTL_LL_TYPE) (b)) + (a); \
   (a) = ((long)(_pp)) & NTL_RADIXM; \
   (t) = (long) (_pp >> NTL_NBITS); \
}

#else

/*
 * This version conforms to the language standard when d is non-negative.
 * Some compilers may emit sub-optimal code, though.
 */



#define zaddmulp(a, b, d, t) { \
   NTL_LL_TYPE _pp = ((NTL_LL_TYPE) (b)) * ((NTL_LL_TYPE) (d)) + ((t)+(a)); \
   (a) = (long) (_pp & NTL_RADIXM); \
   (t) = (long) (_pp >> NTL_NBITS); \
} 


#define zxmulp(a, b, d, t) { \
   NTL_LL_TYPE _pp = ((NTL_LL_TYPE) (b)) * ((NTL_LL_TYPE) (d)) + (t); \
   (a) = (long) (_pp & NTL_RADIXM); \
   (t) = (long) (_pp >> NTL_NBITS); \
} 

#define zaddmulpsq(a,b,t) { \
   NTL_LL_TYPE _pp = ((NTL_LL_TYPE) (b)) * ((NTL_LL_TYPE) (b)) + (a); \
   (a) = (long) (_pp & NTL_RADIXM); \
   (t) = (long) (_pp >> NTL_NBITS); \
}


#endif


#elif (defined(NTL_AVOID_FLOAT))
   

#define zaddmulp(  a,  b,  d,  t) { \
        unsigned long _b1 = b & NTL_RADIXROOTM; \
        unsigned long _d1 = d & NTL_RADIXROOTM; \
        unsigned long _bd,_b1d1,_m,_aa= (a) + (t); \
	unsigned long _ld = (d>>NTL_NBITSH); \
	unsigned long _lb = (b>>NTL_NBITSH); \
 \
        _bd=_lb*_ld; \
        _b1d1=_b1*_d1; \
        _m=(_lb+_b1)*(_ld+_d1) - _bd - _b1d1; \
        _aa += ( _b1d1+ ((_m&NTL_RADIXROOTM)<<NTL_NBITSH)); \
        (t) = (_aa >> NTL_NBITS) + _bd + (_m>>NTL_NBITSH); \
        (a) = _aa & NTL_RADIXM; \
}





#define zxmulp(  a,  b,  d,  t) { \
        unsigned long _b1 = b & NTL_RADIXROOTM; \
        unsigned long _d1 = d & NTL_RADIXROOTM; \
        unsigned long _bd,_b1d1,_m,_aa= (t); \
	unsigned long _ld = (d>>NTL_NBITSH); \
	unsigned long _lb = (b>>NTL_NBITSH); \
 \
        _bd=_lb*_ld; \
        _b1d1=_b1*_d1; \
        _m=(_lb+_b1)*(_ld+_d1) - _bd - _b1d1; \
        _aa += ( _b1d1+ ((_m&NTL_RADIXROOTM)<<NTL_NBITSH)); \
        (t) = (_aa >> NTL_NBITS) + _bd + (_m>>NTL_NBITSH); \
        (a) = _aa & NTL_RADIXM; \
}


#define zaddmulpsq(_a, _b, _t) \
{ \
	long _lb = (_b); \
	long _b1 = (_b) & NTL_RADIXROOTM; \
	long _aa = (_a) + _b1 * _b1; \
 \
	_b1 = (_b1 * (_lb >>= NTL_NBITSH) << 1) + (_aa >> NTL_NBITSH); \
	_aa = (_aa & NTL_RADIXROOTM) + ((_b1 & NTL_RADIXROOTM) << NTL_NBITSH); \
	(_t) = _lb * _lb + (_b1 >> NTL_NBITSH) + (_aa >> NTL_NBITS); \
	(_a) = (_aa & NTL_RADIXM); \
}



#else

/* default long integer arithemtic */
/* various "software pipelining" routines are also defined */


/*
 * The macros CARRY_TYPE and CARRY_CONV are only used in the submul 
 * logic.
 */


#if (defined(NTL_CLEAN_INT))

#define CARRY_TYPE unsigned long
#define CARRY_CONV(x) (-((long)(-x)))

#else

#define CARRY_TYPE long
#define CARRY_CONV(x) (x)

#endif


#if (NTL_BITS_PER_LONG <= NTL_NBITS + 2)

#if (NTL_ARITH_RIGHT_SHIFT && !defined(NTL_CLEAN_INT))
/* value right-shifted is -1..1 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ); \
   _t2 = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   _t1 = (_t1 & NTL_RADIXM) + _a +_t; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}


#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* value shifted is -1..1 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b); \
   long _t1 = _b*_b; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ); \
   _t2 = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   _t1 = (_t1 & NTL_RADIXM) + _a; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/*
 * In the following definition of zam_init, the value _ds is computed so that
 * it is slightly bigger than s*NTL_RADIX_INV.  This has the consequence that
 * the value _hi is equal to floor(_b*_s/NTL_RADIX) or
 * floor(_b*_s/NTL_RADIX) + 1, assuming only that (1) conversion of "small"
 * integer to doubles is exact, (2) multiplication by powers of 2 is exact, and
 * (3) multiplication of two general doubles yields a result with relative
 * error 1/2^{NTL_DOUBLE_PRECISION-1}.  These assumptions are very
 * conservative, and in fact, the IEEE floating point standard would guarantee
 * this result *without* making _ds slightly bigger.
 */

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = ((_s << 1)+1)*(NTL_FRADIX_INV/2.0); \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* value shifted is 0..3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -1..+1 */
#define zsx_loop(a,t,nb) \
{ \
   long  _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}


/* value shifted is -2..+1 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is -1..+1 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is -2..+1 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}

#elif (!defined(NTL_CLEAN_INT))

/* right shift is not arithmetic */

/* value right-shifted is 0..2 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   _t2 = _t2 + ( ((unsigned long) (_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS ); \
   _t1 = (_t1 & NTL_RADIXM) + _a +_t; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}


#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* value shifted is 0..2 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b); \
   long _t1 = _b*_b; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ) - 1; \
   _t2 = _t2 + ( ((unsigned long) (_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS ); \
   _t1 = (_t1 & NTL_RADIXM) + _a; \
   (t) = _t2 + (((unsigned long)_t1) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = ((_s << 1)+1)*(NTL_FRADIX_INV/2.0); \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* value shifted is 0..3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..2 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi += 2; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is 0..2 */
#define zsx_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* value shifted is 0..3 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _hi += 2; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}

#else
/* clean int version */


/* value right-shifted is 0..2 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   unsigned long _t1 = ((unsigned long) _b)*((unsigned long) _d); \
   unsigned long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   _t2 = _t2 + ( (_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS ); \
   _t1 = (_t1 & NTL_RADIXM) + ((unsigned long) _a) + ((unsigned long) _t); \
   (t) = (long) (_t2 + (_t1 >> NTL_NBITS)); \
   (a) = (long) (_t1 & NTL_RADIXM); \
}


#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   unsigned long _t1 =  ((unsigned long) _b)*((unsigned long) _d) + ((unsigned long) _t); \
   unsigned long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = (long) (_t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_t1 & NTL_RADIXM); \
}

/* value shifted is 0..2 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b); \
   unsigned long _t1 = ((unsigned long) _b)*((unsigned long) _b); \
   unsigned long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ) - 1; \
   _t2 = _t2 + ( (_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS ); \
   _t1 = (_t1 & NTL_RADIXM) + ((unsigned long) _a); \
   (t) = (long) (_t2 + (_t1 >> NTL_NBITS)); \
   (a) = (long) (_t1 & NTL_RADIXM); \
}

#define zam_decl double _ds; long _s; unsigned long _hi, _lo;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = ((_s << 1)+1)*(NTL_FRADIX_INV/2.0); \
   _lo = ((unsigned long) _b)*((unsigned long) _s); \
   _hi = (long) (((double) _b)*_ds); \
}

/* value shifted is 0..3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   unsigned long _vv; \
   double _yy; \
   _vv = ((unsigned long) _nb)*((unsigned long)_s); \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + ((unsigned long) _a) + ((unsigned long) _t); \
   _hi--; \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..2 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   unsigned long _vv; \
   double _yy; \
   _vv = ((unsigned long) _nb)*((unsigned long) _s); \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + ((unsigned long) _t); \
   _hi--; \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a); unsigned long _t = (t); long _nb = (nb); \
   unsigned long _vv; \
   double _yy; \
   _vv = ((unsigned long) _nb)*((unsigned long) _s); \
   _yy = ((double) _nb)*_ds; \
   _hi += 2; \
   _lo = ((unsigned long) _a) + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = (long) (_lo & NTL_RADIXM); \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* value shifted is 0..3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + ((unsigned long) _a) + ((unsigned long) _t); \
   _hi--; \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
}

/* value shifted is 0..2 */
#define zsx_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + ((unsigned long) _t); \
   _hi--; \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
}

/* value shifted is 0..3 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a); unsigned long _t = (t); \
   _hi += 2; \
   _lo = ((unsigned long) _a) + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = (long) (_lo & NTL_RADIXM); \
}

#endif
/* end of arithmemtic-right-shift if-then else */
   
#else
/*  NTL_BITS_PER_LONG > NTL_NBITS + 2, and certain optimizations can be
    made.  Useful on 64-bit machines.  */

#if (NTL_ARITH_RIGHT_SHIFT && !defined(NTL_CLEAN_INT))
/* shift is -1..+3 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _a + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ); \
   (t) = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ); \
   (t) = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* shift is -1..+2 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b), _t = (t); \
   long _t1 = _b*_b + _a; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ); \
   (t) = _t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = _s*NTL_FRADIX_INV; \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* shift is -1..+3 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _a + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -1..+2 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -3..+1 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is -1..+3 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is -1..+2 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + _t; \
   (t) = _hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is -3..+1 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _a + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}

#elif (!defined(NTL_CLEAN_INT))
/* right shift is not arithmetic */

/* shift is 0..4 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _a + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   long _t1 =  _b*_d + _t; \
   long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

/* shift is 0..3 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b), _t = (t); \
   long _t1 = _b*_b + _a; \
   long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ) - 1; \
   (t) = _t2 + (((unsigned long)(_t1 - (_t2 << NTL_NBITS))) >> NTL_NBITS); \
   (a) = _t1 & NTL_RADIXM; \
}

#define zam_decl double _ds; long _hi, _lo, _s;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = _s*NTL_FRADIX_INV; \
   _lo = _b*_s; \
   _hi = (long) (((double) _b)*_ds); \
}

/* shift is 0..4 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi--; \
   _lo = _lo + _a + _t; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..3 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi--; \
   _lo = _lo + _t; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..4 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   long _vv; \
   double _yy; \
   _vv = _nb*_s; \
   _yy = ((double) _nb)*_ds; \
   _hi += 3; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..4 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + _a + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is 0..3 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + _t; \
   _hi--; \
   (t) = _hi + (((unsigned long)(_lo - (_hi<<NTL_NBITS))) >> NTL_NBITS); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is 0..4 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _hi += 3; \
   _lo = _a + _t - _lo; \
   (t) = (((unsigned long)(_lo + (_hi<<NTL_NBITS))) >> NTL_NBITS) - _hi; \
   (a) = _lo & NTL_RADIXM; \
}
#else

/* clean int version */

/* shift is 0..4 */
#define zaddmulp(a, b, d, t) \
{ \
   long _a = (a), _b = (b), _d = (d), _t = (t); \
   unsigned long _t1 =  ((unsigned long) _b)*((unsigned long) _d) + ((unsigned long) _a) + ((unsigned long) _t); \
   unsigned long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = (long) (_t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_t1 & NTL_RADIXM); \
}

#define zxmulp(a, b, d, t) \
{ \
   long _b = (b), _d = (d), _t = (t); \
   unsigned long _t1 =  ((unsigned long) _b)*((unsigned long) _d) + ((unsigned long) _t); \
   unsigned long _t2 = (long) ( ((double) _b)*(((double) _d)*NTL_FRADIX_INV) ) - 1; \
   (t) = (long) (_t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_t1 & NTL_RADIXM); \
}

/* shift is 0..3 */
#define zaddmulpsq(a, b, t) \
{ \
   long _a = (a), _b = (b), _t = (t); \
   unsigned long _t1 = ((unsigned long) _b)*((unsigned long) _b) + ((unsigned long) _a); \
   unsigned long _t2 = (long) ( ((double) _b)*(((double) _b)*NTL_FRADIX_INV) ) - 1; \
   (t) = (long) (_t2 + ((_t1 - (_t2 << NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_t1 & NTL_RADIXM); \
}

#define zam_decl double _ds; long _s; unsigned long _hi, _lo;

#define zam_init(b,s) \
{ \
   long _b = (b); \
   _s = (s); \
   _ds = _s*NTL_FRADIX_INV; \
   _lo = ((unsigned long) _b)*((unsigned long) _s); \
   _hi = (long) (((double) _b)*_ds); \
}

/* shift is 0..4 */
#define zam_loop(a,t,nb) \
{ \
   long _a = (a), _t = (t), _nb = (nb); \
   unsigned long _vv; \
   double _yy; \
   _vv = ((unsigned long) _nb)*((unsigned long) _s); \
   _yy = ((double) _nb)*_ds; \
   _hi--; \
   _lo = _lo + ((unsigned long) _a) + ((unsigned long) _t); \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..3 */
#define zsx_loop(a,t,nb) \
{ \
   long _t = (t), _nb = (nb); \
   unsigned long _vv; \
   double _yy; \
   _vv = ((unsigned long) _nb)*((unsigned long) _s); \
   _yy = ((double) _nb)*_ds; \
   _hi--; \
   _lo = _lo + ((unsigned long) _t); \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..4 */
#define zam_subloop(a,t,nb) \
{ \
   long _a = (a); unsigned long _t = (t); long _nb = (nb); \
   unsigned long _vv; \
   double _yy; \
   _vv = ((unsigned long) _nb)*((unsigned long) _s); \
   _yy = ((double) _nb)*_ds; \
   _hi += 3; \
   _lo = ((unsigned long) _a) + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = (long) (_lo & NTL_RADIXM); \
   _lo = _vv; \
   _hi = (long) _yy; \
}

/* shift is 0..4 */
#define zam_finish(a,t) \
{ \
   long _a = (a), _t = (t); \
   _lo = _lo + ((unsigned long) _a) + ((unsigned long) _t); \
   _hi--; \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = _lo & NTL_RADIXM; \
}

/* shift is 0..3 */
#define zsx_finish(a,t) \
{ \
   long _t = (t); \
   _lo = _lo + ((unsigned long) _t); \
   _hi--; \
   (t) = (long) (_hi + ((_lo - (_hi<<NTL_NBITS)) >> NTL_NBITS)); \
   (a) = (long) (_lo & NTL_RADIXM); \
}

/* shift is 0..4 */
#define zam_subfinish(a,t) \
{ \
   long _a = (a); unsigned long _t = (t); \
   _hi += 3; \
   _lo = ((unsigned long) _a) + _t - _lo; \
   (t) = ((_lo + (_hi<<NTL_NBITS)) >> NTL_NBITS) - _hi; \
   (a) = (long) (_lo & NTL_RADIXM); \
}

#endif
/* end of arithmetic-right-shift if-then-else */

#endif
/* end of "NTL_BITS_PER_LONG <= NTL_NBITS + 2" if-then-else */

#endif
/* end of long-integer-implementation if-then-else */






static
void zaddmulone(long *lama, long *lamb)
{ 
   long lami; 
   long lams = 0; 
 
   lams = 0; 
   for (lami = (*lamb++); lami > 0; lami--) { 
      lams += (*lama + *lamb++); 
      *lama++ = lams & NTL_RADIXM; 
      lams >>= NTL_NBITS; 
   } 
   *lama += lams; 
}

#if (NTL_ARITH_RIGHT_SHIFT && !defined(NTL_CLEAN_INT))

static
void zsubmulone(long *lama, long *lamb)
{ 
   long lami; 
   long lams = 0; 
 
   lams = 0; 
   for (lami = (*lamb++); lami > 0; lami--) { 
      lams += (*lama - *lamb++); 
      *lama++ = lams & NTL_RADIXM; 
      lams >>= NTL_NBITS; 
   } 
   *lama += lams; 
}


#else


static
void zsubmulone(long *lama, long *lamb)
{ 
   long lami; 
   long lams = 0; 
 
   lams = 0; 
   for (lami = (*lamb++); lami > 0; lami--) { 
      lams = *lama - *lamb++ - lams; 
      *lama++ = lams & NTL_RADIXM; 
      lams = (lams < 0);
   } 
   *lama -= lams; 
}

#endif


/*
 * definitions of zaddmul, zsxmul, zaddmulsq for the various 
 * long integer implementation options.
 */


#if (defined(NTL_SINGLE_MUL))

static
void zaddmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + (*ama) + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = lo + (*ama) + carry;
   *ama = lo & 0x3FFFFFF;
   carry = hi + (lo >> 26);
   ama++; 
   *ama += carry; 
}

static
void zsxmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = lo + carry;
   *ama = lo & 0x3FFFFFF;
   carry = hi + (lo >> 26);
   ama++; 
   *ama = carry; 
}

static
void zaddmulsq(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = ams;
   double dams = (double) (*amb++);
   double xx;
   double yy;
   unsigned long lo, lo_wd;
   unsigned long hi, hi_wd;

   xx =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams + DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + (*ama) + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
      xx = yy;
   } 
   if (i==1) {
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = lo + (*ama) + carry;
      *ama = lo & 0x3FFFFFF;
      carry = hi + (lo >> 26);
      ama++; 
   }
   *ama += carry; 
}

#elif (defined(NTL_AVOID_FLOAT) || defined(NTL_LONG_LONG))

static
void zaddmul(long lams, long *lama, long *lamb)
{
        long lami;
        long lamcarry = 0;

        for (lami = (*lamb++); lami > 0; lami--)
        {
                zaddmulp(*lama, *lamb, lams, lamcarry);
                lama++;
                lamb++;
        }
        *lama += lamcarry;
}

static
void zsxmul(long lams, long *lama, long *lamb)
{
        long lami;
        long lamcarry = 0;

        for (lami = (*lamb++); lami > 0; lami--)
        {
                zxmulp(*lama, *lamb, lams, lamcarry);
                lama++;
                lamb++;
        }
        *lama = lamcarry;
}

static
void zaddmulsq(long lsqi, long *lsqa, long *lsqb)
{
        long lsqs = *(lsqb);
        long lsqcarry = 0;

        lsqb++;
        for (; lsqi > 0; lsqi--)
        {
                zaddmulp(*lsqa, *lsqb, lsqs, lsqcarry);
                lsqa++;
                lsqb++;
        }
        *lsqa += lsqcarry;
}


#else
/* default long integer arithmetic */

static
void zaddmul(long lams, long *lama, long *lamb)
{ 
   long lami = (*lamb++)-1; 
   long lamcarry = 0; 
   zam_decl;

   zam_init(*lamb, lams);
   lamb++;

 
   for (; lami > 0; lami--) { 
      zam_loop(*lama, lamcarry, *lamb);
      lama++; 
      lamb++; 
   } 
   zam_finish(*lama, lamcarry);
   lama++;
   *lama += lamcarry; 
}



static
void zsxmul(long lams, long *lama, long *lamb)
{ 
   long lami = (*lamb++)-1; 
   long lamcarry = 0; 
   zam_decl;

   zam_init(*lamb, lams);
   lamb++;

 
   for (; lami > 0; lami--) { 
      zsx_loop(*lama, lamcarry, *lamb);
      lama++; 
      lamb++; 
   } 
   zsx_finish(*lama, lamcarry);
   lama++;
   *lama = lamcarry; 
}



static
void zaddmulsq(long lsqi, long *lsqa, long *lsqb)
{ 
   long lsqs; 
   long lsqcarry; 
   zam_decl

   if (lsqi <= 0) return;

   lsqs = *lsqb;
   lsqcarry = 0;

   lsqb++; 
   zam_init(*lsqb, lsqs);
   lsqb++;
   lsqi--;
   for (; lsqi > 0; lsqi--) { 
      zam_loop(*lsqa, lsqcarry, *lsqb);
      lsqa++; 
      lsqb++; 
   } 
   zam_finish(*lsqa, lsqcarry);
   lsqa++;
   *lsqa += lsqcarry; 
}


#endif







/*
 * definition of zsubmul for the various long integer implementation options.
 * Note that zsubmul is only called with a positive first argument.
 */



#if (defined(NTL_SINGLE_MUL))

#if (NTL_ARITH_RIGHT_SHIFT)

static
void zsubmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = (*ama) + carry - lo;
      *ama = lo & 0x3FFFFFF;
      carry = (((long)lo) >> 26) - hi;
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = (*ama) + carry - lo;
   *ama = lo & 0x3FFFFFF;
   carry = (((long)lo) >> 26) - hi;
   ama++; 
   *ama += carry; 
}

#else

static
void zsubmul(long ams, long *ama, long *amb) 
{ 
   long carry = 0; 
   long i = (*amb++);
   double dams = (double) ams;
   double xx;
   double yy;
   unsigned long lo_wd, lo;
   unsigned long hi_wd, hi;

   xx  =  ((double) (*amb++))*dams + DENORMALIZE;
   for (; i > 1; i--) { 
      yy =  ((double) (*amb++))*dams +DENORMALIZE;
      NTL_FetchHiLo(hi_wd, lo_wd, xx);
      lo = lo_wd & 0x3FFFFFF;
      hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
      lo = (*ama) + carry - lo;
      *ama = lo & 0x3FFFFFF;
      carry = ((lo + (1L << 27)) >> 26) - hi - 2;
      ama++; 
      xx = yy;
   } 

   NTL_FetchHiLo(hi_wd, lo_wd, xx);
   lo = lo_wd & 0x3FFFFFF;
   hi = ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
   lo = (*ama) + carry - lo;
   *ama = lo & 0x3FFFFFF;
   carry = ((lo + (1L << 27)) >> 26) - hi - 2;
   ama++; 
   *ama += carry; 
}

#endif



#elif (defined(NTL_AVOID_FLOAT) || (defined(NTL_LONG_LONG) && defined(NTL_CLEAN_INT)))

static void
zsubmul(
        long r,
        _ntl_verylong a,
        _ntl_verylong b
        )
{
        long rd = NTL_RADIX - r;
        long i;
        long carry = NTL_RADIX;

        for (i = (*b++); i > 0; i--)
        {
                zaddmulp(*a, *b, rd, carry);
                a++;
                carry += NTL_RADIXM - (*b++);
        }
        *a += carry - NTL_RADIX; /* unnormalized */
}

#elif (defined(NTL_LONG_LONG))

/*
 * NOTE: the implementation of zaddmulp for the NTL_LONG_LONG option
 * will work on most machines even when the single-precision 
 * multiplicand is negative;  however, the C language standard does
 * not guarantee correct behaviour in this case, which is why the above
 * implementation is used when NTL_CLEAN_INT is set. 
 */

static
void zsubmul(long lams, long *lama, long *lamb)
{
        long lami;
        long lamcarry = 0;

        lams = -lams;

        for (lami = (*lamb++); lami > 0; lami--)
        {
                zaddmulp(*lama, *lamb, lams, lamcarry);
                lama++;
                lamb++;
        }
        *lama += lamcarry;
}


#else
/* default long integer arithmetic */

static
void zsubmul(long lams, long *lama, long *lamb)
{ 
   long lami = (*lamb++)-1; 
   CARRY_TYPE lamcarry = 0; 
   zam_decl;

   zam_init(*lamb, lams);
   lamb++;

 
   for (; lami > 0; lami--) { 
      zam_subloop(*lama, lamcarry, *lamb);
      lama++; 
      lamb++; 
   } 
   zam_subfinish(*lama, lamcarry);
   lama++;
   *lama += CARRY_CONV(lamcarry); 
}

#endif





/*
 *
 * zdiv21 returns quot, numhigh so
 *
 * quot = (numhigh*NTL_RADIX + numlow)/denom;
 * numhigh  = (numhigh*NTL_RADIX + numlow)%denom;
 * Assumes 0 <= numhigh < denom < NTL_RADIX and 0 <= numlow < NTL_RADIX.
 */


#if (defined(NTL_CLEAN_INT))

/*
 * This "clean" version relies on the guaranteed semantics of
 * unsigned integer arithmetic.
 */

#define zdiv21(numhigh, numlow, denom, deninv, quot) \
{ \
   unsigned long udenom = denom; \
   unsigned long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) + \
                        (double) (numlow)) * (deninv)); \
   unsigned long lr21 = (((unsigned long) numhigh) << NTL_NBITS) + \
                        ((unsigned long) numlow)  - udenom*lq21 ; \
 \
   if (lr21 >> (NTL_BITS_PER_LONG-1)) { \
      lq21--; \
      lr21 += udenom; \
   } \
   else if (lr21 >= udenom) { \
      lr21 -= udenom; \
      lq21++; \
   } \
   quot = (long) lq21; \
   numhigh = (long) lr21; \
}


#else

/*
 * This "less clean" version relies on wrap-around semantics for
 * signed integer arithmetic.
 */

#define zdiv21(numhigh, numlow, denom, deninv, quot) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
          + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   if (lr21 < 0) { \
      lq21--; \
      lr21 += denom; \
   } \
   else if (lr21 >= denom) { \
      lr21 -= denom; \
      lq21++; \
   } \
   quot = lq21; \
   numhigh = lr21; \
}

#endif




/*
 * zrem21 behaves just like zdiv21, except the only the remainder is computed.
 */

#if (defined(NTL_CLEAN_INT) || (defined(NTL_AVOID_BRANCHING)  && !NTL_ARITH_RIGHT_SHIFT))
#define NTL_CLEAN_SPMM
#endif

#if (defined(NTL_CLEAN_SPMM)  && !defined(NTL_AVOID_BRANCHING))

#define zrem21(numhigh, numlow, denom, deninv) \
{ \
   unsigned long udenom = denom; \
   unsigned long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) + \
                        (double) (numlow)) * (deninv)); \
   unsigned long lr21 = (((unsigned long) numhigh) << NTL_NBITS) + \
                        ((unsigned long) numlow)  - udenom*lq21 ; \
 \
   if (lr21 >> (NTL_BITS_PER_LONG-1)) { \
      lr21 += udenom; \
   } \
   else if (lr21 >= udenom) { \
      lr21 -= udenom; \
   } \
   numhigh = (long) lr21; \
}

#elif (defined(NTL_CLEAN_SPMM)  && defined(NTL_AVOID_BRANCHING))


#define zrem21(numhigh, numlow, denom, deninv) \
{ \
   unsigned long udenom = denom; \
   unsigned long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) + \
                        (double) (numlow)) * (deninv)); \
   unsigned long lr21 = (((unsigned long) numhigh) << NTL_NBITS) + \
                        ((unsigned long) numlow)  - udenom*lq21 ; \
   lr21 += (-(lr21 >> (NTL_BITS_PER_LONG-1))) & udenom; \
   lr21 -= udenom; \
   lr21 += (-(lr21 >> (NTL_BITS_PER_LONG-1))) & udenom; \
   numhigh = (long) lr21; \
}


#elif (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING))

#define zrem21(numhigh, numlow, denom, deninv) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
          + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   lr21 += (lr21 >> (NTL_BITS_PER_LONG-1)) & denom; \
   lr21 -= denom; \
   lr21 += (lr21 >> (NTL_BITS_PER_LONG-1)) & denom; \
   numhigh = lr21; \
}

#else

#define zrem21(numhigh, numlow, denom, deninv) \
{ \
   long lr21; \
   long lq21 = (long) (((NTL_FRADIX * (double) (numhigh)) \
      + (double) (numlow)) * (deninv)); \
   long lp21; \
   MulLo(lp21, lq21, denom); \
   lr21 = (numhigh << NTL_NBITS) + numlow - lp21; \
   if (lr21 < 0) lr21 += denom; \
   else if (lr21 >= denom) lr21 -= denom; \
   numhigh = lr21; \
}

#endif


long _ntl_zmaxalloc(_ntl_verylong x)
{
   if (!x) 
      return 0; 
   else
      return (x[-1] >> 1) - 1;
}
      

void _ntl_zsetlength(_ntl_verylong *v, long len)
{
   _ntl_verylong x = *v;

   if (len < 0)
      zhalt("negative size allocation in _ntl_zsetlength");

   if (NTL_OVERFLOW(len, NTL_NBITS, 0))
      zhalt("size too big in _ntl_zsetlength");

#ifdef L_TO_G_CHECK_LEN
   /* this makes sure that numbers don't get too big for GMP */
   if (len >= (1L << (NTL_BITS_PER_INT-4)))
      zhalt("size too big for GMP");
#endif

   if (x) {
      long oldlen = x[-1];
      long fixed = oldlen & 1;

      oldlen = oldlen >> 1;

      if (fixed) {
         if (len > oldlen) 
            zhalt("internal error: can't grow this _ntl_verylong");
         else
            return;
      }

      if (len <= oldlen) return;

      len++;  /* always allocate at least one more than requested */

      oldlen = (long) (oldlen * 1.2); /* always increase by at least 20% */
      if (len < oldlen)
         len = oldlen;

      /* round up to multiple of MIN_SETL */
      len = ((len+(MIN_SETL-1))/MIN_SETL)*MIN_SETL; 

      /* test len again */
      if (NTL_OVERFLOW(len, NTL_NBITS, 0))
         zhalt("size too big in _ntl_zsetlength");

      x[-1] = len << 1;
      if (!(x = (_ntl_verylong)NTL_REALLOC(&(x[-1]), 
                  len, sizeof(long), 2*sizeof(long)))) {
         zhalt("reallocation failed in _ntl_zsetlength");
      }
   }
   else {
      len++; /* as above, always allocate one more than requested */
      len = ((len+(MIN_SETL-1))/MIN_SETL)*MIN_SETL; 

      /* test len again */
      if (NTL_OVERFLOW(len, NTL_NBITS, 0))
         zhalt("size too big in _ntl_zsetlength");


      if (!(x = (_ntl_verylong)NTL_MALLOC(len, 
                  sizeof(long), 2*sizeof(long)))) {
         zhalt("allocation failed in _ntl_zsetlength");
      }
      x[0] = len << 1;
      x[1] = 1;
      x[2] = 0;
   }

   *v = x+1;
}

void _ntl_zfree(_ntl_verylong *x)
{
   _ntl_verylong y;

   if (!(*x))
      return;

   if ((*x)[-1] & 1)
      zhalt("Internal error: can't free this _ntl_verylong");

   y = (*x - 1);
   free((void*)y);
   *x = 0;
}






long _ntl_zround_correction(_ntl_verylong a, long k, long residual)
{
   long direction;
   long p;
   long sgn;
   long bl;
   long wh;
   long i;

   if (a[0] > 0)
      sgn = 1;
   else
      sgn = -1;

   p = k - 1;
   bl = (p/NTL_NBITS);
   wh = 1L << (p - NTL_NBITS*bl);
   bl++;

   if (a[bl] & wh) {
      /* bit is 1...we have to see if lower bits are all 0
         in order to implement "round to even" */

      if (a[bl] & (wh - 1)) 
         direction = 1;
      else {
         i = bl - 1;
         while (i > 0 && a[i] == 0) i--;
         if (i > 0)
            direction = 1;
         else
            direction = 0;
      }

      /* use residual to break ties */

      if (direction == 0 && residual != 0) {
         if (residual == sgn)
            direction = 1;
         else 
            direction = -1;
      }

      if (direction == 0) {
         /* round to even */

         wh = wh << 1;
         if (wh == NTL_RADIX) {
            wh = 1;
            bl++;
         }

         if (a[bl] & wh)
            direction = 1;
         else
            direction = -1;
      }
   }
   else
      direction = -1;

   if (direction == 1)
      return sgn;

   return 0;
}



double _ntl_zdoub_aux(_ntl_verylong n)
{
   double res;
   long i;

   if (!n)
      return ((double) 0);
   if ((i = n[0]) < 0)
      i = -i;
   res = (double) (n[i--]);
   for (; i; i--)
      res = res * NTL_FRADIX + (double) (n[i]);
   if (n[0] > 0)
      return (res);
   return (-res);
}



double _ntl_zdoub(_ntl_verylong n)
{
   static _ntl_verylong tmp = 0;

   long s;
   long shamt;
   long correction;
   double x;

   s = _ntl_z2log(n);
   shamt = s - NTL_DOUBLE_PRECISION;

   if (shamt <= 0)
      return _ntl_zdoub_aux(n);

   _ntl_zrshift(n, shamt, &tmp);

   correction = _ntl_zround_correction(n, shamt, 0);

   if (correction) _ntl_zsadd(tmp, correction, &tmp);

   x = _ntl_zdoub_aux(tmp);

   x = _ntl_ldexp(x, shamt);

   return x;
}


double _ntl_zlog(_ntl_verylong n)
{
   static _ntl_verylong tmp = 0;

   static double log_2;
   static long init = 0;

   long s;
   long shamt;
   long correction;
   double x;

   if (!init) {
      log_2 = log(2.0);
      init = 1;
   }

   if (_ntl_zsign(n) <= 0)
      zhalt("log argument <= 0");

   s = _ntl_z2log(n);
   shamt = s - NTL_DOUBLE_PRECISION;

   if (shamt <= 0)
      return log(_ntl_zdoub_aux(n));

   _ntl_zrshift(n, shamt, &tmp);

   correction = _ntl_zround_correction(n, shamt, 0);

   if (correction) _ntl_zsadd(tmp, correction, &tmp);

   x = _ntl_zdoub_aux(tmp);

   return log(x) + shamt*log_2;
}



void _ntl_zdoubtoz(double a, _ntl_verylong *xx)
{
   _ntl_verylong x;
   long neg, i, t, sz;


   a = floor(a);

   if (!_ntl_IsFinite(&a))
      zhalt("_ntl_zdoubtoz: attempt to convert non-finite value");


   if (a < 0) {
      a = -a;
      neg = 1;
   }
   else
      neg = 0;

   if (a == 0) {
      _ntl_zzero(xx);
      return;
   }

   sz = 1;
   a = a*NTL_FRADIX_INV;

   while (a >= 1) {
      a = a*NTL_FRADIX_INV;
      sz++;
   }
         
   x = *xx;
   if (MustAlloc(x, sz)) {
      _ntl_zsetlength(&x, sz);
      *xx = x;
   }

   for (i = sz; i > 0; i--) {
      a = a*NTL_FRADIX;
      t = (long) a;
      x[i] = t;
      a = a - t;
   }

   x[0] = (neg ? -sz : sz);
}

void _ntl_zzero(_ntl_verylong *aa)
{
   if (!(*aa)) _ntl_zsetlength(aa, 1);
   (*aa)[0] =  1;
   (*aa)[1] =  0;
}

/* same as _ntl_zzero, except does not unnecessarily allocate space */

void _ntl_zzero1(_ntl_verylong *aa)
{
   if (!(*aa)) return;
   (*aa)[0] =  1;
   (*aa)[1] =  0;
}

void _ntl_zone(_ntl_verylong *aa)
{
   if (!(*aa)) _ntl_zsetlength(aa, 1);
   (*aa)[0] =  1;
   (*aa)[1] =  1;
}

void _ntl_zcopy(_ntl_verylong a, _ntl_verylong *bb)
{
   long i;
   _ntl_verylong b = *bb;

   if (!a) {
      _ntl_zzero(bb);
      return;
   }
   if (a != b) {
      if ((i = *a) < 0)
         i = (-i);
      if (MustAlloc(b, i)) {
         _ntl_zsetlength(&b, i);
         *bb = b;
      }
      for (; i >= 0; i--)
         *b++ = *a++;
   }
}

/* same as _ntl_zcopy, but does not unnecessarily allocate space */

void _ntl_zcopy1(_ntl_verylong a, _ntl_verylong *bb)
{
   long i;
   _ntl_verylong b = *bb;

   if (!a) {
      _ntl_zzero1(bb);
      return;
   }
   if (a != b) {
      if ((i = *a) < 0)
         i = (-i);
      if (MustAlloc(b, i)) {
         _ntl_zsetlength(&b, i);
         *bb = b;
      }
      for (; i >= 0; i--)
         *b++ = *a++;
   }
}

void _ntl_zintoz(long d, _ntl_verylong *aa)
{
   long i;
   long anegative;
   unsigned long d1, d2;
   _ntl_verylong a = *aa;

   anegative = 0;
   if (d < 0) {
      anegative = 1;
      d1 = - ((unsigned long) d);  /* careful: avoid overflow */
   }
   else
      d1 = d;

   i = 0;
   d2 = d1;
   do {
      d2 >>= NTL_NBITS;
      i++;
   }
   while (d2 > 0);

   if (MustAlloc(a, i)) {
      _ntl_zsetlength(&a, i);
      *aa = a;
   }

   i = 0;
   a[1] = 0;
   while (d1 > 0) {
      a[++i] = d1 & NTL_RADIXM;
      d1 >>= NTL_NBITS;
   }
   if (i > 0)
      a[0] = i;
   else
      a[0] = 1;

   if (anegative)
      a[0] = (-a[0]);
}


/* same as _ntl_zintoz, but does not unnecessarily allocate space */

void _ntl_zintoz1(long d, _ntl_verylong *aa)
{
   long i;
   long anegative;
   unsigned long d1, d2;
   _ntl_verylong a = *aa;

   if (!d && !a) return;

   anegative = 0;
   if (d < 0) {
      anegative = 1;
      d1 = - ((unsigned long) d);  /* careful: avoid overlow */
   }
   else
      d1 = d;

   i = 0;
   d2 = d1;
   do {
      d2 >>= NTL_NBITS;
      i++;
   }
   while (d2 > 0);

   if (MustAlloc(a, i)) {
      _ntl_zsetlength(&a, i);
      *aa = a;
   }

   i = 0;
   a[1] = 0;
   while (d1 > 0) {
      a[++i] = d1 & NTL_RADIXM;
      d1 >>= NTL_NBITS;
   }
   if (i > 0)
      a[0] = i;
   else
      a[0] = 1;

   if (anegative)
      a[0] = (-a[0]);
}


void _ntl_zuintoz(unsigned long d, _ntl_verylong *aa)
{
   long i;
   unsigned long d1, d2;
   _ntl_verylong a = *aa;

   d1 = d;
   i = 0;
   d2 = d1;
   do {
      d2 >>= NTL_NBITS;
      i++;
   }
   while (d2 > 0);

   if (MustAlloc(a, i)) {
      _ntl_zsetlength(&a, i);
      *aa = a;
   }

   i = 0;
   a[1] = 0;
   while (d1 > 0) {
      a[++i] = d1 & NTL_RADIXM;
      d1 >>= NTL_NBITS;
   }
   if (i > 0)
      a[0] = i;
   else
      a[0] = 1;
}


unsigned long _ntl_ztouint(_ntl_verylong a)
{
   unsigned long d;
   long sa;

   if (!a)
      return (0);

   if ((sa = *a) < 0)
      sa = -sa;

   d = (unsigned long) (*(a += sa));
   while (--sa) {
      d <<= NTL_NBITS;
      d += (unsigned long) (*(--a));
   }

   if ((*(--a)) < 0)
      return (-d);
   return (d);
}


long _ntl_ztoint(_ntl_verylong a)
{
   unsigned long res = _ntl_ztouint(a);
   return NTL_ULONG_TO_LONG(res);
}



long _ntl_zcompare(_ntl_verylong a, _ntl_verylong b)
{
   long sa;
   long sb;

   if (!a) {
      if (!b)
         return (0);
      if (b[0] < 0)
         return (1);
      if (b[0] > 1)
         return (-1);
      if (b[1])
         return (-1);
      return (0);
   }
   if (!b) {
      if (a[0] < 0)
         return (-1);
      if (a[0] > 1)
         return (1);
      if (a[1])
         return (1);
      return (0);
   }

   if ((sa = *a) > (sb = *b))
      return (1);
   if (sa < sb)
      return (-1);
   if (sa < 0)
      sa = (-sa);
   a += sa;
   b += sa;
   for (; sa; sa--) {
      long diff = *a - *b;

      if (diff > 0) {
         if (sb < 0)
            return (-1);
         return (1);
      }
      if (diff < 0) {
         if (sb < 0)
            return (1);
         return (-1);
      }

      a--;
      b--;
   }
   return (0);
}

void _ntl_znegate(_ntl_verylong *aa)
{
   _ntl_verylong a = *aa;

   if (!a)
      return;
   if (a[1] || a[0] != 1)
      a[0] = (-a[0]);
}

void _ntl_zsadd(_ntl_verylong a, long d, _ntl_verylong *b)
{
   static _ntl_verylong x = 0;

   _ntl_zintoz(d, &x);
   _ntl_zadd(a, x, b);
}


void
_ntl_zadd(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *cc)
{
   long sa;
   long sb;
   long anegative;
   _ntl_verylong c;
   long a_alias, b_alias;

   if (!a) {
      if (b)
         _ntl_zcopy(b, cc);
      else
         _ntl_zzero(cc);
      return;
   }

   if (!b) {
      _ntl_zcopy(a, cc);
      return;
   }

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   if ((anegative = ((sa = a[0]) < 0)) == ((sb = b[0]) < 0)) {
      /* signs a and b are the same */
      _ntl_verylong pc;
      long carry;
      long i;
      long maxab;

      if (anegative) {
         sa = -sa;
         sb = -sb;
      }

      if (sa < sb) {
         i = sa;
         maxab = sb;
      }
      else {
         i = sb;
         maxab = sa;
      }

      if (MustAlloc(c, maxab+1)) {
         _ntl_zsetlength(&c, maxab + 1);
         if (a_alias) a = c; 
         if (b_alias) b = c; 
         *cc = c;
      }

      pc = c;
      carry = 0;

      do {
         long t = (*(++a)) + (*(++b)) + carry;
         carry = t >> NTL_NBITS;
         *(++pc) = t & NTL_RADIXM;
         i--;
      } while (i);

      i = sa-sb;
      if (!i) goto i_exit;

      if (i < 0) {
         i = -i;
         a = b;
      }

      if (!carry) goto carry_exit;

      for (;;) {
         long t = (*(++a)) + 1;
         carry = t >> NTL_NBITS;
         *(++pc) = t & NTL_RADIXM;
         i--;
         if (!i) goto i_exit;
         if (!carry) goto carry_exit;
      }

      i_exit:
      if (carry) {
         *(++pc) = 1;
         maxab++;
      }
      *c = anegative ? -maxab : maxab;
      return;

      carry_exit:
      if (pc != a) {
         do {
            *(++pc) = *(++a);
            i--;
         } while (i);
      }
      *c = anegative ? -maxab : maxab;
   }
   else {
      /* signs a and b are different...use _ntl_zsub */

      if (anegative) {
         a[0] = -sa;
         _ntl_zsub(b, a, cc);
         if (!a_alias) a[0] = sa;
      }
      else {
         b[0] = -sb;
         _ntl_zsub(a, b, cc);
         if (!b_alias) b[0] = sb;
      }
   }
}

void
_ntl_zsub(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *cc)
{
   long sa;
   long sb;
   long anegative;
   long a_alias, b_alias;
   _ntl_verylong c;

   if (!b) {
      if (a)
         _ntl_zcopy(a, cc);
      else
         _ntl_zzero(cc);
      return;
   }

   if (!a) {
      _ntl_zcopy(b, cc);
      _ntl_znegate(cc);
      return;
   }

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   if ((anegative = ((sa = a[0]) < 0)) == ((sb = b[0]) < 0)) {
      /* signs agree */

      long i, carry, *pc;

      if (anegative) {
         sa = -sa;
         sb = -sb;
      }

      carry = sa - sb;
      if (!carry) {
         long *aa = a + sa;
         long *bb = b + sa;
         
         i = sa;
         while (i && !(carry = (*aa - *bb))) {
            aa--; bb--; i--;
         }
      }

      if (!carry) {
         _ntl_zzero(cc);
         return;
      }

      if (carry < 0) {
         { long t = sa; sa = sb; sb = t; }
         { long t = a_alias; a_alias = b_alias; b_alias = t; }
         { long *t = a; a = b; b = t; }
         anegative = !anegative;
      }

      if (MustAlloc(c, sa)) {
         _ntl_zsetlength(&c, sa);
         /* must have !a_alias */
         if (b_alias) b = c;
         *cc = c;
      }

      i = sb;
      carry = 0;
      pc = c;

      do {
#if (!NTL_ARITH_RIGHT_SHIFT || defined(NTL_CLEAN_INT))
         long t = (*(++a)) - (*(++b)) - carry;
         carry = (t < 0);
#else
         long t = (*(++a)) - (*(++b)) + carry;
         carry = t >> NTL_NBITS;
#endif
         *(++pc) = t & NTL_RADIXM;
         i--;
      } while (i);

      i = sa-sb;
      while (carry) {
         long t = (*(++a)) - 1;
#if (!NTL_ARITH_RIGHT_SHIFT || defined(NTL_CLEAN_INT))
         carry = (t < 0);
#else
         carry = t >> NTL_NBITS;
#endif
         *(++pc) = t & NTL_RADIXM;
         i--;
      }

      if (i) {
         if (pc != a) {
            do {
               *(++pc) = *(++a);
               i--;
            } while (i);
         }
      }
      else {
         while (sa > 1 && *pc == 0) { sa--; pc--; } 
      }

      if (anegative) sa = -sa;
      *c = sa;
   }
   else {
      /* signs of a and b are different...use _ntl_zadd */

      if (anegative) {
         a[0] = -sa;
         _ntl_zadd(a, b, cc);
         if (!a_alias) a[0] = sa;
         c = *cc;
         c[0] = -c[0];
      }
      else {
         b[0] = -sb;
         _ntl_zadd(a, b, cc);
         if (!b_alias) b[0] = sb;
      }
   }
}


void
_ntl_zsmul(_ntl_verylong a, long d, _ntl_verylong *bb)
{
   long sa;
   long anegative, bnegative;
   _ntl_verylong b;
   long a_alias;


   if (d == 2) {
      _ntl_z2mul(a, bb);
      return;
   }


   if ((d >= NTL_RADIX) || (d <= -NTL_RADIX)) {
      static _ntl_verylong x = 0;
      _ntl_zintoz(d,&x);
      _ntl_zmul(a, x, bb);
      return;
   }

   if (!a || (a[0] == 1 && a[1] == 0)) {
      _ntl_zzero(bb);
      return;
   }

   if (!d) {
      _ntl_zzero(bb);
      return;
   }

   /* both inputs non-zero */

   anegative = 0;
   bnegative = 0;

   if ((sa = a[0]) < 0) {
      anegative = 1;
      a[0] = sa = (-sa);
      if (d < 0)
         d = (-d);
      else
         bnegative = 1;
   }
   else if (bnegative = (d < 0))
      d = (-d);

   b = *bb;
   a_alias = (a == b);

   if (MustAlloc(b, sa + 1)) {
      _ntl_zsetlength(&b, sa + 1);
      if (a_alias) a = b;
      *bb = b;
   }

   zsxmul(d, b+1, a);

   sa++;
   while ((sa > 1) && (!(b[sa])))
      sa--;
   b[0] = sa;

   if (bnegative)
      b[0] = (-b[0]);

   if (anegative && !a_alias)
      a[0] = -a[0];
}

void _ntl_zsubpos(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *cc)
{
   long sa;
   long sb;

   long *c, *pc;
   long i, carry;

   long b_alias;

   if (!b) {
      if (a)
         _ntl_zcopy(a, cc);
      else
         _ntl_zzero(cc);
      return;
   }

   if (!a) {
      _ntl_zzero(cc);
      return;
   }

   sa = a[0];
   sb = b[0];

   c = *cc;
   b_alias = (b == c);

   if (MustAlloc(c, sa)) {
      _ntl_zsetlength(&c, sa);
      if (b_alias) b = c;
      *cc = c;
   }

   i = sb;
   carry = 0;
   pc = c;

   while (i) {
#if (!NTL_ARITH_RIGHT_SHIFT || defined(NTL_CLEAN_INT))
      long t = (*(++a)) - (*(++b)) - carry;
      carry = (t < 0);
#else
      long t = (*(++a)) - (*(++b)) + carry;
      carry = t >> NTL_NBITS;
#endif
      *(++pc) = t & NTL_RADIXM;
      i--;
   }

   i = sa-sb;
   while (carry) {
      long t = (*(++a)) - 1;
#if (!NTL_ARITH_RIGHT_SHIFT || defined(NTL_CLEAN_INT))
      carry = (t < 0);
#else
      carry = t >> NTL_NBITS;
#endif
      *(++pc) = t & NTL_RADIXM;
      i--;
   }

   if (i) {
      if (pc != a) {
         do {
            *(++pc) = *(++a);
            i--;
         } while (i);
      }
   }
   else {
      while (sa > 1 && *pc == 0) { sa--; pc--; } 
   }

   *c = sa;
}





static long *kmem = 0;     /* globals for Karatsuba */
static long max_kmem = 0;


/* These cross-over points were estimated using
   a Sparc-10, a Sparc-20, and a Pentium-90. */

#if (defined(NTL_SINGLE_MUL))
#define KARX (18)
#else
#define KARX (16) 
#endif

/* Auxilliary routines for Karatsuba */


static
void kar_fold(long *T, long *b, long hsa)
{
   long sb, *p2, *p3, i, carry;

   sb = *b;
   p2 = b + hsa;
   p3 = T;
   carry = 0;

   for (i = sb-hsa; i>0; i--) {
      long t = (*(++b)) + (*(++p2)) + carry;
      carry = t >> NTL_NBITS;
      *(++p3) = t & NTL_RADIXM;
   }

   for (i = (hsa << 1) - sb; i>0; i--) {
      long t = (*(++b)) + carry;
      carry = t >> NTL_NBITS;
      *(++p3) = t & NTL_RADIXM;
   }

   if (carry) {
      *(++p3) = carry;
      *T = hsa + 1;
   }
   else
      *T = hsa;
}

static
void kar_sub(long *T, long *c)
{
   long i, carry;

   i = *c;
   carry = 0;

   while (i>0) {
#if (!NTL_ARITH_RIGHT_SHIFT || defined(NTL_CLEAN_INT))
      long t = (*(++T)) - (*(++c)) - carry;
      carry = (t < 0);
#else
      long t = (*(++T)) - (*(++c)) + carry;
      carry = t >> NTL_NBITS;
#endif
      *T = t & NTL_RADIXM;
      i--;
   }

   while (carry) {
      long t = (*(++T)) - 1;
#if (!NTL_ARITH_RIGHT_SHIFT || defined(NTL_CLEAN_INT))
      carry = (t < 0);
#else
      carry = t >> NTL_NBITS;
#endif
      *T = t & NTL_RADIXM;
   }
}


static
void kar_add(long *c, long *T, long hsa)
{
   long i, carry;

   c += hsa;
   i = *T;
   while (T[i] == 0 && i > 0) i--;
   carry = 0;

   while (i>0) {
      long t = (*(++c)) + (*(++T)) + carry;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
      i--;
   }

   while (carry) {
      long t = (*(++c)) + 1;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
   }
}

static
void kar_fix(long *c, long *T, long hsa)
{
   long i, carry, s;

   s = *T;

   i = hsa;
   while (i>0) {
      *(++c) = *(++T);
      i--;
   }


   i = s - hsa;
   carry = 0;

   while (i > 0) {
      long t = (*(++c)) + (*(++T)) + carry;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
      i--;
   }

   while (carry) {
      long t = (*(++c)) + 1;
      carry = t >> NTL_NBITS;
      *c = t & NTL_RADIXM;
   }
}
      
   

static
void kar_mul(long *c, long *a, long *b, long *stk)
{
   long sa, sb, sc;

   if (*a < *b) { long *t = a; a = b; b = t; }

   sa = *a;
   sb = *b;
   sc = sa + sb;

   if (sb < KARX) {
      /* classic algorithm */

      long *pc, i, *pb;

      pc = c;
      for (i = sc; i; i--) {
         pc++;
         *pc = 0;
      }
   
      pc = c;
      pb = b;
      for (i = sb; i; i--) {
         pb++;
         pc++;
         zaddmul(*pb, pc, a);
      }
   }
   else {
      long hsa = (sa + 1) >> 1;

      if (hsa < sb) {
         /* normal case */

         long *T1, *T2, *T3;

         /* allocate space */

         T1 = stk;  stk += hsa + 2;
         T2 = stk;  stk += hsa + 2;
         T3 = stk;  stk += (hsa << 1) + 3;

         if (stk-kmem > max_kmem) zhalt("internal error: kmem overflow");

         /* compute T1 = a_lo + a_hi */

         kar_fold(T1, a, hsa);

         /* compute T2 = b_lo + b_hi */

         kar_fold(T2, b, hsa);
         
         /* recursively compute T3 = T1 * T2 */

         kar_mul(T3, T1, T2, stk);

         /* recursively compute a_hi * b_hi into high part of c */
         /* and subtract from T3 */

         {
            long olda, oldb;

            olda = a[hsa];  a[hsa] = sa-hsa;
            oldb = b[hsa];  b[hsa] = sb-hsa;

            kar_mul(c + (hsa << 1), a + hsa, b + hsa, stk);
            kar_sub(T3, c + (hsa << 1));

            a[hsa] = olda;
            b[hsa] = oldb;
         }

         /* recursively compute a_lo*b_lo into low part of c */
         /* and subtract from T3 */

         *a = hsa;
         *b = hsa;

         kar_mul(c, a, b, stk);
         kar_sub(T3, c);

         *a = sa;
         *b = sb;

         /* finally, add T3 * NTL_RADIX^{hsa} to c */

         kar_add(c, T3, hsa);
      }
      else {
         /* degenerate case */

         long *T;
         
         T = stk;  stk += sb + hsa + 1;

         if (stk-kmem > max_kmem) zhalt("internal error: kmem overflow");

         /* recursively compute b*a_hi into high part of c */
         {
            long olda;

            olda = a[hsa];  a[hsa] = sa-hsa;
            kar_mul(c + hsa, a + hsa, b, stk);
            a[hsa] = olda;
         }

         /* recursively compute b*a_lo into T */

         *a = hsa;
         kar_mul(T, a, b, stk);
         *a = sa;

         /* fix-up result */

         kar_fix(c, T, hsa);
      }
   }

   /* normalize result */

   while (c[sc] == 0 && sc > 1) sc--;
   *c = sc;
}



#if (defined(NTL_SINGLE_MUL))
#define KARSX (36)
#else
#define KARSX (32) 
#endif

static
void kar_sq(long *c, long *a, long *stk)
{
   long sa, sc;

   sa = *a;
   sc = sa << 1;

   if (sa < KARSX) {
      /* classic algorithm */

      long carry, i, j, *pc;

      pc = c;
      for (i = sc; i; i--) {
         pc++;
         *pc = 0;
      }

      carry = 0;
      i = 0;
      for (j = 1; j <= sa; j++) {
         unsigned long uncar;
         long t;

         i += 2;
         uncar = ((unsigned long) carry) + (((unsigned long) c[i-1]) << 1);
         t = uncar & NTL_RADIXM;
         zaddmulpsq(t, a[j], carry);
         c[i-1] = t;
         zaddmulsq(sa-j, c+i, a+j);
         uncar =  (uncar >> NTL_NBITS) + (((unsigned long) c[i]) << 1);
         uncar += ((unsigned long) carry);
         carry = uncar >> NTL_NBITS;
         c[i] = uncar & NTL_RADIXM;
      }
   }
   else {
      long hsa = (sa + 1) >> 1;
      long *T1, *T2, olda;

      T1 = stk;  stk += hsa + 2;
      T2 = stk;  stk += (hsa << 1) + 3;

      if (stk-kmem > max_kmem) zhalt("internal error: kmem overflow");

      kar_fold(T1, a, hsa);
      kar_sq(T2, T1, stk);

      olda = a[hsa];  a[hsa] = sa - hsa;
      kar_sq(c + (hsa << 1), a + hsa, stk);
      kar_sub(T2, c + (hsa << 1));
      a[hsa] = olda;

      *a = hsa;
      kar_sq(c, a, stk);
      kar_sub(T2, c);
      *a = sa;

      kar_add(c, T2, hsa);
   }

   while (c[sc] == 0 && sc > 1) sc--;
   *c = sc;
}


void _ntl_zmul(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *cc)
{
   static _ntl_verylong mem = 0;
   _ntl_verylong c = *cc;

   if (!a || (a[0] == 1 && a[1] == 0) || !b || (b[0] == 1 && b[1] == 0)) {
      _ntl_zzero(cc);
      return;
   }

   if (a == b) {
      if (a == c) {
         _ntl_zcopy(a, &mem);
         a = mem;
      }

      _ntl_zsq(a, cc);
   }
   else {
      long aneg, bneg, sa, sb, sc;

      if (a == c) {
         _ntl_zcopy(a, &mem);
         a = mem;
      }
      else if (b == c) {
         _ntl_zcopy(b, &mem);
         b = mem;
      }

      sa = *a;
      if (sa < 0) {
         *a = sa = -sa;
         aneg = 1;
      }
      else
         aneg = 0;

      sb = *b;
      if (*b < 0) {
         *b = sb =  -sb;
         bneg = 1;
      }
      else
         bneg = 0;

      sc = sa + sb;
      if (MustAlloc(c, sc)) {
         _ntl_zsetlength(&c, sc);
         *cc = c;
      }

      /* we optimize for *very* small numbers,
       * avoiding all function calls and loops */

      if (sa <= 3 && sb <= 3) {
         long carry, d;

         switch (sa) {
         case 1:
            switch (sb) {
            case 1:
               carry = 0;
               zxmulp(c[1], a[1], b[1], carry);
               c[2] = carry;
               break;
            case 2:
               carry = 0;
               d = a[1];
               zxmulp(c[1], b[1], d, carry);
               zxmulp(c[2], b[2], d, carry);
               c[3] = carry;
               break;
            case 3:
               carry = 0;
               d = a[1];
               zxmulp(c[1], b[1], d, carry);
               zxmulp(c[2], b[2], d, carry);
               zxmulp(c[3], b[3], d, carry);
               c[4] = carry;
               break;
            }
            break;

         case 2:
            switch (sb) {
            case 1:
               carry = 0;
               d = b[1];
               zxmulp(c[1], a[1], d, carry);
               zxmulp(c[2], a[2], d, carry);
               c[3] = carry;
               break;
            case 2:
               carry = 0;
               d = b[1];
               zxmulp(c[1], a[1], d, carry);
               zxmulp(c[2], a[2], d, carry);
               c[3] = carry;
               carry = 0;
               d = b[2];
               zaddmulp(c[2], a[1], d, carry);
               zaddmulp(c[3], a[2], d, carry);
               c[4] = carry;
               break;
             case 3:
               carry = 0;
               d = a[1];
               zxmulp(c[1], b[1], d, carry);
               zxmulp(c[2], b[2], d, carry);
               zxmulp(c[3], b[3], d, carry);
               c[4] = carry;
               carry = 0;
               d = a[2];
               zaddmulp(c[2], b[1], d, carry);
               zaddmulp(c[3], b[2], d, carry);
               zaddmulp(c[4], b[3], d, carry);
               c[5] = carry;
               break;
            }
            break;

         case 3:
            switch (sb) {
            case 1:
               carry = 0;
               d = b[1];
               zxmulp(c[1], a[1], d, carry);
               zxmulp(c[2], a[2], d, carry);
               zxmulp(c[3], a[3], d, carry);
               c[4] = carry;
               break;
            case 2:
               carry = 0;
               d = b[1];
               zxmulp(c[1], a[1], d, carry);
               zxmulp(c[2], a[2], d, carry);
               zxmulp(c[3], a[3], d, carry);
               c[4] = carry;
               carry = 0;
               d = b[2];
               zaddmulp(c[2], a[1], d, carry);
               zaddmulp(c[3], a[2], d, carry);
               zaddmulp(c[4], a[3], d, carry);
               c[5] = carry;
               break;
            case 3:
               carry = 0;
               d = b[1];
               zxmulp(c[1], a[1], d, carry);
               zxmulp(c[2], a[2], d, carry);
               zxmulp(c[3], a[3], d, carry);
               c[4] = carry;
               carry = 0;
               d = b[2];
               zaddmulp(c[2], a[1], d, carry);
               zaddmulp(c[3], a[2], d, carry);
               zaddmulp(c[4], a[3], d, carry);
               c[5] = carry;
               carry = 0;
               d = b[3];
               zaddmulp(c[3], a[1], d, carry);
               zaddmulp(c[4], a[2], d, carry);
               zaddmulp(c[5], a[3], d, carry);
               c[6] = carry;
               break;
            }
         }

         if (c[sc] == 0) sc--;
         if (aneg != bneg) sc = -sc;
         *c = sc;
         if (aneg) *a = -sa;
         if (bneg) *b = -sb;

         return;
         
      }


#ifdef NTL_GMP_HACK

      if (_ntl_gmp_hack && sa >= NTL_GMP_MUL_CROSS && sb >= NTL_GMP_MUL_CROSS) {
         long gsa = L_TO_G(sa); 
         long gsb = L_TO_G(sb);

         mp_limb_t *a_p, *b_p, *c_p, *end_p;

         NTL_GSPACE((gsa + gsb) << 1);
         a_p = gspace_data;
         b_p = a_p + gsa;
         c_p = b_p + gsb;
         end_p = c_p + (gsa + gsb);

         lip_to_gmp(a+1, a_p, sa);
         lip_to_gmp(b+1, b_p, sb); 

         if (!a_p[gsa-1]) gsa--;
         if (!b_p[gsb-1]) gsb--;
         end_p[-1] = 0;
         end_p[-2] = 0;

         if (gsa >= gsb)
            mpn_mul(c_p, a_p, gsa, b_p, gsb);
         else
            mpn_mul(c_p, b_p, gsb, a_p, gsa);

         gmp_to_lip(c+1, c_p, sc);

         if (!c[sc]) sc--; 
         if (aneg != bneg) sc = -sc;
         c[0] = sc;
         if (aneg) *a = - *a;
         if (bneg) *b = - *b;
         return;
      }



#endif

      if (*a < KARX || *b < KARX) {
         /* classic algorithm */

         long i, *pc;

         pc = c;
         for (i = sc; i; i--) {
            pc++;
            *pc = 0;
         }
   
         pc = c;
   
         if (*a >= *b) {
            long *pb = b;
            for (i = *pb; i; i--) {
               pb++;
               pc++;
               zaddmul(*pb, pc, a);
            }
         }
         else {
            long *pa = a;
            for (i = *pa; i; i--) {
               pa++; 
               pc++;
               zaddmul(*pa, pc, b);
            }
         }

         while (c[sc] == 0 && sc > 1) sc--;
         if (aneg != bneg) sc = -sc;
         c[0] = sc;
         if (aneg) *a = - *a;
         if (bneg) *b = - *b;
      }
      else {
         /* karatsuba */

         long n, hn, sp;

         if (*a < *b)
            n = *b;
         else
            n = *a;

         sp = 0;
         do {
            hn = (n + 1) >> 1;
            sp += (hn << 2) + 7;
            n = hn+1;
         } while (n >= KARX);

         if (sp > max_kmem) {
            if (max_kmem == 0) 
               kmem = (long *) NTL_MALLOC(sp, sizeof(long), 0);
            else
               kmem = (long *) NTL_REALLOC(kmem, sp, sizeof(long), 0);

            max_kmem = sp;
            if (!kmem) zhalt("out of memory in karatsuba");
         }

         kar_mul(c, a, b, kmem);
         if (aneg != bneg) *c = - *c;
         if (aneg) *a = - *a;
         if (bneg) *b = - *b;
      }
   }
}


void _ntl_zsq(_ntl_verylong a, _ntl_verylong *cc)
{
   static _ntl_verylong mem = 0;
   _ntl_verylong c = *cc;
   long sa, aneg, sc;

   if (!a || (a[0] == 1 && a[1] == 0)) {
      _ntl_zzero(cc);
      return;
   }

   if (a == c) {
      _ntl_zcopy(a, &mem);
      a = mem;
   }

   sa = *a;

   if (*a < 0) {
      *a = sa = -sa;
      aneg = 1;
   }
   else
      aneg = 0;

   sc = (sa) << 1;
   if (MustAlloc(c, sc)) {
      _ntl_zsetlength(&c, sc);
      *cc = c;
   }

   if (sa <= 3) {
      long carry, d;

      switch (sa) {
      case 1:
         carry = 0;
         zxmulp(c[1], a[1], a[1], carry);
         c[2] = carry;
         break;

      case 2:
         carry = 0;
         d = a[1];
         zxmulp(c[1], a[1], d, carry);
         zxmulp(c[2], a[2], d, carry);
         c[3] = carry;
         carry = 0;
         d = a[2];
         zaddmulp(c[2], a[1], d, carry);
         zaddmulp(c[3], a[2], d, carry);
         c[4] = carry;
         break;

      case 3:
         carry = 0;
         d = a[1];
         zxmulp(c[1], a[1], d, carry);
         zxmulp(c[2], a[2], d, carry);
         zxmulp(c[3], a[3], d, carry);
         c[4] = carry;
         carry = 0;
         d = a[2];
         zaddmulp(c[2], a[1], d, carry);
         zaddmulp(c[3], a[2], d, carry);
         zaddmulp(c[4], a[3], d, carry);
         c[5] = carry;
         carry = 0;
         d = a[3];
         zaddmulp(c[3], a[1], d, carry);
         zaddmulp(c[4], a[2], d, carry);
         zaddmulp(c[5], a[3], d, carry);
         c[6] = carry;
         break;
      }

      if (c[sc] == 0) sc--;
      *c = sc;
      if (aneg) *a = -sa;

      return;
   }

#ifdef NTL_GMP_HACK

      if (_ntl_gmp_hack && sa >= NTL_GMP_SQR_CROSS) {
         long gsa = L_TO_G(sa); 

         mp_limb_t *a_p, *c_p, *end_p;;

         NTL_GSPACE(gsa + gsa + gsa);
         a_p = gspace_data;
         c_p = a_p + gsa;
         end_p = c_p + (gsa + gsa);

         lip_to_gmp(a+1, a_p, sa);

         if (!a_p[gsa-1]) gsa--;
         end_p[-1] = 0;
         end_p[-2] = 0;

         mpn_mul(c_p, a_p, gsa, a_p, gsa);

         gmp_to_lip(c+1, c_p, sc);

         if (!c[sc]) sc--; 
         c[0] = sc;
         if (aneg) *a = - *a;
         return;
      }


#endif

   if (sa < KARSX) {
      /* classic algorithm */

      long carry, i, j, *pc;

      pc = c;
      for (i = sc; i; i--) {
         pc++;
         *pc = 0;
      }

      carry = 0;
      i = 0;
      for (j = 1; j <= sa; j++) {
         unsigned long uncar;
         long t;

         i += 2;
         uncar = ((unsigned long) carry) + (((unsigned long) c[i-1]) << 1);
         t = uncar & NTL_RADIXM;
         zaddmulpsq(t, a[j], carry);
         c[i-1] = t;
         zaddmulsq(sa-j, c+i, a+j);
         uncar =  (uncar >> NTL_NBITS) + (((unsigned long) c[i]) << 1);
         uncar += ((unsigned long) carry);
         carry = uncar >> NTL_NBITS;
         c[i] = uncar & NTL_RADIXM;
      }


      while (c[sc] == 0 && sc > 1) sc--;
      c[0] = sc;
      if (aneg) *a = - *a;
   }
   else {
      /* karatsuba */

      long n, hn, sp;

      n = *a;

      sp = 0;
      do {
         hn = (n + 1) >> 1;
         sp += hn + hn + hn + 5;
         n = hn+1;
      } while (n >= KARSX);

      if (sp > max_kmem) {
         if (max_kmem == 0) 
            kmem = (long *) NTL_MALLOC(sp, sizeof(long), 0);
         else
            kmem = (long *) NTL_REALLOC(kmem, sp, sizeof(long), 0);

         max_kmem = sp;
         if (!kmem) zhalt("out of memory in karatsuba");
      }

      kar_sq(c, a, kmem);
      if (aneg) *a = - *a;
   }
}




long _ntl_zsdiv(_ntl_verylong a, long d, _ntl_verylong *bb)
{
   long sa;
   _ntl_verylong b = *bb;

   if (!d) {
      zhalt("division by zero in _ntl_zsdiv");
   }

   if (!a) {
      _ntl_zzero(bb);
      return (0);
   }


   if (d == 2) {
      long is_odd = a[1] & 1;
      long fix = (a[0] < 0) & is_odd;
      _ntl_zrshift(a, 1, bb);
      if (fix) _ntl_zsadd(*bb, -1, bb);
      return is_odd;
   } 


   if ((sa = a[0]) < 0)
      sa = (-sa);

   /* if b aliases a, then b won't move */
   _ntl_zsetlength(&b, sa);
   *bb = b;

   if ((d >= NTL_RADIX) || (d <= -NTL_RADIX)) {
      static _ntl_verylong zd = 0;
      static _ntl_verylong zb = 0;

      _ntl_zintoz(d, &zb);
      _ntl_zdiv(a, zb, &b, &zd);
      *bb = b;
      return (_ntl_ztoint(zd));
   }
   else {
      long den = d;
      double deninv;
      long carry = 0;
      long i;
      long flag = (*a < 0 ? 2 : 0) | (den < 0 ? 1 : 0);

      if (den < 0)
         den = -den;
      deninv = 1.0 / ((double) den);

      if (a[sa] < den && sa > 1)
         carry = a[sa--];

      for (i = sa; i; i--) {
         zdiv21(carry, a[i], den, deninv, b[i]);
      }

      while ((sa > 1) && (!(b[sa])))
         sa--;
      b[0] = sa;

      if (flag) {
         if (flag <= 2) {
            if (!carry)
               _ntl_znegate(&b);
            else {
               _ntl_zsadd(b, 1, &b);
               b[0] = -b[0];
               if (flag == 1)
                  carry = carry - den;
               else
                  carry = den - carry;
               *bb = b;
            }
         }
         else
            carry = -carry;
      }

      return (carry);
   }
}

long _ntl_zsmod(_ntl_verylong a, long d)
{
   long sa;

   if (!a) {
      return (0);
   }

   if (d == 2) return (a[1] & 1);

   if (!d) {
      zhalt("division by zero in _ntl_zsdiv");
   }

   if ((sa = a[0]) < 0)
      sa = (-sa);

   if ((d >= NTL_RADIX) || (d <= -NTL_RADIX)) {
      static _ntl_verylong zd = 0;
      static _ntl_verylong zb = 0;

      _ntl_zintoz(d, &zb);
      _ntl_zmod(a, zb, &zd);
      return (_ntl_ztoint(zd));
   }
   else {
      long den = d;
      double deninv;
      long carry = 0;
      long i;
      long flag = (*a < 0 ? 2 : 0) | (den < 0 ? 1 : 0);

      if (den < 0)
         den = -den;
      deninv = 1.0 / ((double) den);

      if (a[sa] < den && sa > 1)
         carry = a[sa--];

      for (i = sa; i; i--) {
         zrem21(carry, a[i], den, deninv);
      }

      if (flag) {
         if (flag <= 2) {
            if (carry) {
               if (flag == 1)
                  carry = carry - den;
               else
                  carry = den - carry;
            }
         }
         else
            carry = -carry;
      }

      return (carry);
   }
}

void _ntl_zmultirem(_ntl_verylong a, long n, long* dd, long *rr)
{
   long j;
   long sa;

   if (!a || (a[0] == 1 && a[1] == 0)) {
      for (j = 0; j < n; j++) rr[j] = 0;
      return;
   }

   sa = a[0];

   for (j = 0; j < n; j++) {
      long den = dd[j];
      double deninv;
      long carry = 0;
      long i;
      long lsa = sa;

      deninv = 1.0 / ((double) den);

      if (a[lsa] < den && lsa > 1)
         carry = a[lsa--];

      for (i = lsa; i; i--) {
         zrem21(carry, a[i], den, deninv);
      }

      rr[j] = carry;
   }
}




#if (defined(NTL_SINGLE_MUL))

#define REDUCE_CROSS 8

void _ntl_zmultirem2(_ntl_verylong a, long n, long* dd, double **ttbl, long *rr)

{
   long d; 
   double *tbl;
   long ac0, ac1, ac2;
   long *ap; 
   double *tp;
   long sa, i, j, k;
   unsigned long pc0, pc1, lo_wd, hi_wd;
   long carry;
   double xx, yy, dinv;

   if (!a ||  a[0] < REDUCE_CROSS || a[0] >= NTL_RADIX) {
      _ntl_zmultirem(a, n, dd, rr);
      return;
   }

   sa = a[0];

   for (i = 0; i < n; i++) {
      d = dd[i];
      tbl = ttbl[i];

      ac0 = a[1];
      ac1 = 0;
      ac2 = 0;

      ap = &a[2];
      tp = &tbl[1];

      k = sa-1;
      while (k) {
         j = k;
         if (j > 64) j = 64;
         k -= j;

         pc0 = 0;
         pc1 = 0;

         xx =  ((double) (*ap++))* (*tp++) + DENORMALIZE;
         for (; j > 1; j--) {
            yy =  ((double) (*ap++))*(*tp++) 
                  + DENORMALIZE;
            NTL_FetchHiLo(hi_wd, lo_wd, xx);
            pc0 +=  lo_wd & 0x3FFFFFF;
            pc1 += ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
            xx = yy;
         }
         NTL_FetchHiLo(hi_wd, lo_wd, xx);
         pc0 +=  lo_wd & 0x3FFFFFF;
         pc1 += ((hi_wd<<6)|(lo_wd>>26)) & 0x3FFFFFF;
         pc1 += (pc0 >> 26);

         ac0 += (pc0 & 0x3FFFFFF);
         if (ac0 > NTL_RADIX) {
            ac0 -= NTL_RADIX;
            carry = 1;
         }
         else
            carry = 0;

         ac1 += carry + (pc1 & 0x3FFFFFF);
         if (ac1 > NTL_RADIX) {
            ac1 -= NTL_RADIX;
            carry = 1;
         }
         else
            carry = 0;

         ac2 += carry + (pc1 >> 26);
      }

      carry = 0;
      dinv = ((double) 1)/((double) d);
      if (ac2 >= d) {
         zrem21(carry, ac2, d, dinv);
      }
      else
         carry = ac2;

      zrem21(carry, ac1, d, dinv);
      zrem21(carry, ac0, d, dinv);

      rr[i] = carry;
   }
}


#elif (defined(NTL_TBL_REM))

#if (defined(NTL_LONG_LONG))

/* This version uses the double-word long type directly.
 * It's a little faster that the other one.
 * It accumlates 8 double-word products before stepping
 * a higher-level accumulator.
 */



void _ntl_zmultirem3(_ntl_verylong a, long n, long* dd, long **ttbl, long *rr)
{
   long sa, i, j, d, *tbl, ac0, ac1, ac2, *ap, *tp, k, carry;
   double dinv;
   NTL_LL_TYPE acc;

   if (!a || a[0] < 8 || a[0] >= NTL_RADIX) {
      _ntl_zmultirem(a, n, dd, rr);
      return;
   }

   sa = a[0];
 
   for (i = 0; i < n; i++) {
      d = dd[i];
      tbl = ttbl[i];
      acc = a[1];
      ac2 = 0;
      ap = &a[2];
      tp = &tbl[1];

      k = sa - 7;

      for (j = 0; j < k; j += 7) {
         acc += ((NTL_LL_TYPE) ap[j+0]) * ((NTL_LL_TYPE) tp[j+0]);
         acc += ((NTL_LL_TYPE) ap[j+1]) * ((NTL_LL_TYPE) tp[j+1]);
         acc += ((NTL_LL_TYPE) ap[j+2]) * ((NTL_LL_TYPE) tp[j+2]);
         acc += ((NTL_LL_TYPE) ap[j+3]) * ((NTL_LL_TYPE) tp[j+3]);
         acc += ((NTL_LL_TYPE) ap[j+4]) * ((NTL_LL_TYPE) tp[j+4]);
         acc += ((NTL_LL_TYPE) ap[j+5]) * ((NTL_LL_TYPE) tp[j+5]);
         acc += ((NTL_LL_TYPE) ap[j+6]) * ((NTL_LL_TYPE) tp[j+6]);
         ac2 += (long) (acc >> (2*NTL_NBITS));
         acc &= (((NTL_LL_TYPE) 1) << (2*NTL_NBITS)) - ((NTL_LL_TYPE) 1);
      }

      k = sa - 1;

      for (; j < k; j++)
         acc += ((NTL_LL_TYPE) ap[j+0]) * ((NTL_LL_TYPE) tp[j+0]);

      ac2 += (long) (acc >> (2*NTL_NBITS));
      acc &= (((NTL_LL_TYPE) 1) << (2*NTL_NBITS)) - ((NTL_LL_TYPE) 1);

      ac0 = (long) (acc & ( (((NTL_LL_TYPE) 1) << (NTL_NBITS)) - ((NTL_LL_TYPE) 1) ));
      ac1 = (long) (acc >> NTL_NBITS);
      

      carry = 0;
      dinv = ((double) 1)/((double) d);
      if (ac2 >= d) {
         zrem21(carry, ac2, d, dinv);
      }
      else
         carry = ac2;

      zrem21(carry, ac1, d, dinv);
      zrem21(carry, ac0, d, dinv);

      rr[i] = carry;
   }
}

#else

void _ntl_zmultirem3(_ntl_verylong a, long n, long* dd, long **ttbl, long *rr)
{
   long sa, i, d, *tbl, ac0, ac1, ac2, *ap, *tp, k, t, carry;
   double dinv;

   if (!a || a[0] < 8 || a[0] >= NTL_RADIX) {
      _ntl_zmultirem(a, n, dd, rr);
      return;
   }

   sa = a[0];
 
   for (i = 0; i < n; i++) {
      d = dd[i];
      tbl = ttbl[i];
      ac0 = a[1];
      ac1 = 0;
      ac2 = 0;
      ap = &a[2];
      tp = &tbl[1];
      k = sa-1;
      
      while (k) {
         zxmulp(t, *ap, *tp, ac0);
         ac1 += ac0;
         ac2 += ac1 >> NTL_NBITS;
         ac1 &= NTL_RADIXM;
         ac0 = t;
         k--;
         ap++;
         tp++;
      }

      carry = 0;
      dinv = ((double) 1)/((double) d);
      if (ac2 >= d) {
         zrem21(carry, ac2, d, dinv);
      }
      else
         carry = ac2;

      zrem21(carry, ac1, d, dinv);
      zrem21(carry, ac0, d, dinv);

      rr[i] = carry;
   }
}




#endif



#endif




long _ntl_zsfastrem(_ntl_verylong a, long d)
/* assumes a >= 0, and 0 < d < NTL_RADIX */
/* computes a % d */

{
   long sa;

   if (!a || (a[0] == 1 && a[1] == 0)) {
      return 0;
   }

   sa = a[0];

   {
      long den = d;
      double deninv = ((double)1)/((double)den);
      long carry = 0;
      long i;
      long lsa = sa;

      if (a[lsa] < den && lsa > 1)
         carry = a[lsa--];
      for (i = lsa; i; i--) {
         zrem21(carry, a[i], den, deninv);
      }
      return carry;
   }
}



void _ntl_zdiv(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *qq, _ntl_verylong *rr)
{
   long sa, sb, sq, i;
   long sign;
   long q1;
   long *rp;
   double btopinv, aux;
   static _ntl_verylong q=0, r=0;

   if (!b || (((sb=b[0]) == 1) && (!b[1]))) {
      zhalt("division by zero in _ntl_zdiv");
   }

   if (!a || (((sa=a[0]) == 1) && (!a[1]))) {
      _ntl_zzero(qq);
      if (rr) _ntl_zzero(rr);
      return;
   }


   if (sb == 1) {
      long t1 = _ntl_zsdiv(a, b[1], qq);
      if (rr) _ntl_zintoz(t1, rr);
      return;
   }

   if (sb == -1) {
      long t1 = _ntl_zsdiv(a, -b[1], qq);
      if (rr) _ntl_zintoz(t1, rr);
      return;
   }

   sign = 0;
   if (sa < 0) {
      a[0] = sa = -sa;
      sign = 2;
   }

   if (sb < 0) {
      b[0] = sb = -sb;
      sign |= 1;
   }


   sq = sa-sb+1;

   if (sq <= 0) {
      _ntl_zcopy(a, &r);
      _ntl_zzero(&q);
      goto done;
   }

#ifdef NTL_GMP_HACK
   if (_ntl_gmp_hack && sb >= NTL_GMP_DIV_CROSS && sq >= NTL_GMP_DIV_CROSS) {
      GT_INIT

      lip_to_mpz(a, gt_1);
      lip_to_mpz(b, gt_2);
      mpz_fdiv_qr(gt_3, gt_4, gt_1, gt_2);
      mpz_to_lip(&q, gt_3);
      mpz_to_lip(&r, gt_4);

      goto done;
   }
#endif

   _ntl_zsetlength(&q, sq);
   _ntl_zsetlength(&r, sa+1);
 
   _ntl_zcopy(a, &r);
   rp = &r[sa+1];
   *rp = 0;

   r[0] = 0; /* this streamlines the last evaluation of aux */

   btopinv = b[sb]*NTL_FRADIX + b[sb-1];
   if (sb > 2) 
      btopinv = NTL_FRADIX / (btopinv*NTL_FRADIX + b[sb-2]);
   else
      btopinv = 1.0 / btopinv;

   
   aux = btopinv*(rp[-1]*NTL_FRADIX + rp[-2]);
   if (aux >= NTL_FRADIX)
      aux = NTL_FRADIX-1;

   for (i = sq; i >= 1; i--, rp--) {
      q1 = (long) aux;
      if (q1) {
         zsubmul(q1, &r[i], b);
      }

      while (rp[0] < 0) {
         zaddmulone(&r[i], b);
         q1--;
      }

      while (rp[0] > 0) {
         zsubmulone(&r[i], b);
         q1++;
      }

      aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
      while (aux > NTL_FRADIX - 16) {
         /* q1 might be too small */
         if (aux >= NTL_FRADIX)
            aux = NTL_FRADIX-1;


         zsubmulone(&r[i], b);
         if (rp[0] < 0) {
            /* oops...false alarm! */
            zaddmulone(&r[i], b);
            break;
         }
         else {
            q1++;
            aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
         }
      }

      q[i] = q1;
   }

   while (sq > 1 && q[sq] == 0) sq--;
   q[0] = sq;

   i = sb;
   while (i > 1 && r[i] == 0) i--;
   r[0] = i;

done:
   if (sign)
   {
      if (sign <= 2)
      {
         if (!(r[1]) && (r[0] == 1))
            _ntl_znegate(&q);
         else
         {
            _ntl_zsadd(q, 1, &q);
            _ntl_znegate(&q);
            if (sign == 1)
               _ntl_zsub(r, b, &r);
            else
               _ntl_zsub(b, r, &r);
         }
      }
      else
         _ntl_znegate(&r);

      if (sign & 2)
         a[0] = -sa;

      if (sign & 1)
         b[0] = -sb;
   }

   _ntl_zcopy(q, qq);
   if (rr) _ntl_zcopy(r, rr);
}

void
_ntl_zmod(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *rr)
{
   long sa, sb, sq, i;
   long sign;
   long q1;
   long *rp;
   double btopinv, aux;
   static _ntl_verylong r=0;

   if (!b || (((sb=b[0]) == 1) && (!b[1]))) {
      zhalt("division by zero in _ntl_zdiv");
   }

   if (!a || (((sa=a[0]) == 1) && (!a[1]))) {
      _ntl_zzero(rr);
      return;
   }


   if (sb == 1) {
      _ntl_zintoz(_ntl_zsmod(a, b[1]), rr);
      return;
   }

   if (sb == -1) {
      _ntl_zintoz(_ntl_zsmod(a, -b[1]), rr);
      return;
   }

   sign = 0;
   if (sa < 0) {
      a[0] = sa = -sa;
      sign = 2;
   }

   if (sb < 0) {
      b[0] = sb = -sb;
      sign |= 1;
   }


   sq = sa-sb+1;

   if (sq <= 0) {
      _ntl_zcopy(a, &r);
      goto done;
   }

#ifdef NTL_GMP_HACK
   if (_ntl_gmp_hack && sb >= NTL_GMP_REM_CROSS && sq >= NTL_GMP_REM_CROSS) {
      GT_INIT

      lip_to_mpz(a, gt_1);
      lip_to_mpz(b, gt_2);
      mpz_fdiv_r(gt_4, gt_1, gt_2);
      mpz_to_lip(&r, gt_4);

      goto done;
   }
#endif

   _ntl_zsetlength(&r, sa+1);
 
   _ntl_zcopy(a, &r);
   rp = &r[sa+1];
   *rp = 0;

   r[0] = 0; /* this streamlines the last evaluation of aux */

   btopinv = b[sb]*NTL_FRADIX + b[sb-1];
   if (sb > 2) 
      btopinv = NTL_FRADIX / (btopinv*NTL_FRADIX + b[sb-2]);
   else
      btopinv = 1.0 / btopinv;

   
   aux = btopinv*(rp[-1]*NTL_FRADIX + rp[-2]);
   if (aux >= NTL_FRADIX)
      aux = NTL_FRADIX-1;

   for (i = sq; i >= 1; i--, rp--) {
      q1 = (long) aux;
      if (q1) {
         zsubmul(q1, &r[i], b);
      }

      while (rp[0] < 0) {
         zaddmulone(&r[i], b);
      }

      while (rp[0] > 0) {
         zsubmulone(&r[i], b);
      }

      aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
      while (aux > NTL_FRADIX - 16) {
         /* q1 might be too small */
         if (aux >= NTL_FRADIX)
            aux = NTL_FRADIX-1;


         zsubmulone(&r[i], b);
         if (rp[0] < 0) {
            /* oops...false alarm! */
            zaddmulone(&r[i], b);
            break;
         }
         else {
            aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
         }
      }
   }

   i = sb;
   while (i > 1 && r[i] == 0) i--;
   r[0] = i;

done:
   if (sign)
   {
      if (sign <= 2)
      {
         if (!(r[1]) && (r[0] == 1))
            /* no op */;
         else
         {
            if (sign == 1)
               _ntl_zsub(r, b, &r);
            else
               _ntl_zsub(b, r, &r);
         }
      }
      else
         _ntl_znegate(&r);

      if (sign & 2)
         a[0] = -sa;

      if (sign & 1)
         b[0] = -sb;
   }

   _ntl_zcopy(r, rr);
}

void
_ntl_zquickmod(_ntl_verylong *rr, _ntl_verylong b)
{
   long sa, sb, sq, i;
   long q1;
   long *rp;
   double btopinv, aux;
   _ntl_verylong r;

   sb = b[0];

   r = *rr;

   if (!r || (((sa=r[0]) == 1) && (!r[1]))) {
      _ntl_zzero(rr);
      return;
   }


   if (sb == 1) {
      _ntl_zintoz(_ntl_zsmod(r, b[1]), rr);
      return;
   }

   sq = sa-sb+1;

   if (sq <= 0) {
      return;
   } 

#ifdef NTL_GMP_HACK
   if (_ntl_gmp_hack && sb >= NTL_GMP_QREM_CROSS && sq >= NTL_GMP_QREM_CROSS) {
      GT_INIT

      lip_to_mpz(r, gt_1);
      lip_to_mpz(b, gt_2);
      mpz_fdiv_r(gt_4, gt_1, gt_2);
      mpz_to_lip(rr, gt_4);
      return;
   }
#endif

   _ntl_zsetlength(rr, sa+1);
   r = *rr;
 
   rp = &r[sa+1];
   *rp = 0;

   r[0] = 0; /* this streamlines the last evaluation of aux */

   btopinv = b[sb]*NTL_FRADIX + b[sb-1];
   if (sb > 2) 
      btopinv = NTL_FRADIX / (btopinv*NTL_FRADIX + b[sb-2]);
   else
      btopinv = 1.0 / btopinv;

   
   aux = btopinv*(rp[-1]*NTL_FRADIX + rp[-2]);
   if (aux >= NTL_FRADIX)
      aux = NTL_FRADIX-1;

   for (i = sq; i >= 1; i--, rp--) {
      q1 = (long) aux;
      if (q1) {
         zsubmul(q1, &r[i], b);
      }

      while (rp[0] < 0) {
         zaddmulone(&r[i], b);
      }

      while (rp[0] > 0) {
         zsubmulone(&r[i], b);
      }

      aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
      while (aux > NTL_FRADIX - 16) {
         /* q1 might be too small */
         if (aux >= NTL_FRADIX)
            aux = NTL_FRADIX-1;


         zsubmulone(&r[i], b);
         if (rp[0] < 0) {
            /* oops...false alarm! */
            zaddmulone(&r[i], b);
            break;
         }
         else {
            aux = btopinv*((rp[-1]*NTL_FRADIX + rp[-2])*NTL_FRADIX + rp[-3]);
         }
      }
   }

   i = sb;
   while (i > 1 && r[i] == 0) i--;
   r[0] = i;
}


void
_ntl_zaddmod(
	_ntl_verylong a,
	_ntl_verylong b,
	_ntl_verylong n,
	_ntl_verylong *c
	)
{
	if (*c != n) {
		_ntl_zadd(a, b, c);
		if (_ntl_zcompare(*c, n) >= 0)
			_ntl_zsubpos(*c, n, c);
	}
	else {
		static _ntl_verylong mem = 0;

		_ntl_zadd(a, b, &mem);
		if (_ntl_zcompare(mem, n) >= 0)
			_ntl_zsubpos(mem, n, c);
		else
			_ntl_zcopy(mem, c);
	}
}

void
_ntl_zsubmod(
	_ntl_verylong a,
	_ntl_verylong b,
	_ntl_verylong n,
	_ntl_verylong *c
	)
{
	static _ntl_verylong mem = 0;
	long cmp;

	if ((cmp=_ntl_zcompare(a, b)) < 0) {
		_ntl_zadd(n, a, &mem);
		_ntl_zsubpos(mem, b, c);
	} else if (!cmp) 
		_ntl_zzero(c);
	else 
		_ntl_zsubpos(a, b, c);
}

void
_ntl_zsmulmod(
	_ntl_verylong a,
	long d,
	_ntl_verylong n,
	_ntl_verylong *c
	)
{
	static _ntl_verylong mem = 0;

	_ntl_zsmul(a, d, &mem);
	_ntl_zquickmod(&mem, n);
	_ntl_zcopy(mem, c);
}




void
_ntl_zmulmod(
	_ntl_verylong a,
	_ntl_verylong b,
	_ntl_verylong n,
	_ntl_verylong *c
	)
{
	static _ntl_verylong mem = 0;

	_ntl_zmul(a, b, &mem);
	_ntl_zquickmod(&mem, n);
	_ntl_zcopy(mem, c);
}

void
_ntl_zsqmod(
	_ntl_verylong a,
	_ntl_verylong n,
	_ntl_verylong *c
	)
{
	static _ntl_verylong mem = 0;

	_ntl_zsq(a, &mem);
	_ntl_zquickmod(&mem, n);
	_ntl_zcopy(mem, c);
}


void
_ntl_zinvmod(
	_ntl_verylong a,
	_ntl_verylong n,
	_ntl_verylong *c
	)
{
	if (_ntl_zinv(a, n, c))
		zhalt("undefined inverse in _ntl_zinvmod");
}


static long 
zxxeucl(
   _ntl_verylong ain,
   _ntl_verylong nin,
   _ntl_verylong *invv,
   _ntl_verylong *uu
   )
{
   static _ntl_verylong a = 0;
   static _ntl_verylong n = 0;
   static _ntl_verylong q = 0;
   static _ntl_verylong w = 0;
   static _ntl_verylong x = 0;
   static _ntl_verylong y = 0;
   static _ntl_verylong z = 0;
   _ntl_verylong inv = *invv;
   _ntl_verylong u = *uu;
   long diff;
   long ilo;
   long sa;
   long sn;
   long temp;
   long e;
   long fast;
   long parity;
   long gotthem;
   _ntl_verylong pin;
   _ntl_verylong p;
   long i;
   long try11;
   long try12;
   long try21;
   long try22;
   long got11;
   long got12;
   long got21;
   long got22;
   double hi;
   double lo;
   double dt;
   double fhi, fhi1;
   double flo, flo1;
   double num;
   double den;
   double dirt;


   _ntl_zsetlength(&a, (e = (ain[0] > nin[0] ? ain[0] : nin[0])));
   _ntl_zsetlength(&n, e);
   _ntl_zsetlength(&q, e);
   _ntl_zsetlength(&w, e);
   _ntl_zsetlength(&x, e);
   _ntl_zsetlength(&y, e);
   _ntl_zsetlength(&z, e);
   _ntl_zsetlength(&inv, e);
   *invv = inv;
   _ntl_zsetlength(&u, e);
   *uu = u;

   fhi1 = 1.0 + ((double) 32.0)/NTL_FDOUBLE_PRECISION;
   flo1 = 1.0 - ((double) 32.0)/NTL_FDOUBLE_PRECISION;

   fhi = 1.0 + ((double) 8.0)/NTL_FDOUBLE_PRECISION;
   flo = 1.0 - ((double) 8.0)/NTL_FDOUBLE_PRECISION;

   pin = &ain[0];
   p = &a[0];
   for (i = (*pin); i >= 0; i--)
      *p++ = *pin++;
   pin = &nin[0];
   p = &n[0];
   for (i = (*pin); i >= 0; i--)
      *p++ = *pin++;
   inv[0] = 1;
   inv[1] = 1;
   w[0] = 1;
   w[1] = 0;
   while (n[0] > 1 || n[1] > 0)
   {
      gotthem = 0;
      sa = a[0];
      sn = n[0];
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = a[0];
         p = &a[sa];
         num = ((double) (*p)) * NTL_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_FRADIX;
         if (sa > 2)
            num += (*(p - 1));
         sn = n[0];
         p = &n[sn];
         den = (double) (*p) * NTL_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_FRADIX;
         if (sn > 2)
            den += (*(p - 1));
         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_FRADIX;
            lo *= NTL_FRADIX;
         }
         try11 = 1;
         try12 = 0;
         try21 = 0;
         try22 = 1;
         parity = 1;
         fast = 1; 
         while (fast > 0)
         {
            parity = 1 - parity;
            if (hi >= NTL_FRADIX)
               fast = 0;
            else
            {
               ilo = (long)lo;
               dirt = hi - ilo;
               if (dirt < 1.0/NTL_FDOUBLE_PRECISION || !ilo || ilo < (long)hi)
                  fast = 0;
               else
               {
                  dt = lo-ilo;
                  lo = flo / dirt;
                  if (dt > 1.0/NTL_FDOUBLE_PRECISION)
                     hi = fhi / dt;
                  else
                     hi = NTL_FRADIX;
                  temp = try11;
                  try11 = try21;
                  if ((NTL_RADIX - temp) / ilo < try21)
                     fast = 0;
                  else
                     try21 = temp + ilo * try21;
                  temp = try12;
                  try12 = try22;
                  if ((NTL_RADIX - temp) / ilo < try22)
                     fast = 0;
                  else
                     try22 = temp + ilo * try22;
                  if ((fast > 0) && (parity > 0))
                  {
                     gotthem = 1;
                     got11 = try11;
                     got12 = try12;
                     got21 = try21;
                     got22 = try22;
                  }
               }
            }
         }
      }
      if (gotthem)
      {
         _ntl_zsmul(inv, got11, &x);
         _ntl_zsmul(w, got12, &y);
         _ntl_zsmul(inv, got21, &z);
         _ntl_zsmul(w, got22, &w);
         _ntl_zadd(x, y, &inv);
         _ntl_zadd(z, w, &w);
         _ntl_zsmul(a, got11, &x);
         _ntl_zsmul(n, got12, &y);
         _ntl_zsmul(a, got21, &z);
         _ntl_zsmul(n, got22, &n);
         _ntl_zsub(x, y, &a);
         _ntl_zsub(n, z, &n);
      }
      else
      {
         _ntl_zdiv(a, n, &q, &a);
         _ntl_zmul(q, w, &x);
         _ntl_zadd(inv, x, &inv);
         if (a[0] > 1 || a[1] > 0)
         {
            _ntl_zdiv(n, a, &q, &n);
            _ntl_zmul(q, inv, &x);
            _ntl_zadd(w, x, &w);
         }
         else
         {
            p = &a[0];
            pin = &n[0];
            for (i = (*pin); i >= 0; i--)
               *p++ = *pin++;
            n[0] = 1;
            n[1] = 0;
            _ntl_zcopy(w, &inv);
            _ntl_znegate(&inv);
         }
      }
   }

   if ((a[0] == 1) && (a[1] == 1))
      e = 0;
   else 
      e = 1;

   p = &u[0];
   pin = &a[0];
   for (i = (*pin); i >= 0; i--)
      *p++ = *pin++;
   *invv = inv;
   *uu = u;

   return (e);
}

long 
_ntl_zinv(
	_ntl_verylong ain,
	_ntl_verylong nin,
	_ntl_verylong *invv
	)
{
	static _ntl_verylong u = 0;
	static _ntl_verylong v = 0;
	long sgn;


	if (_ntl_zscompare(nin, 1) <= 0) {
		zhalt("InvMod: second input <= 1");
	}

	sgn = _ntl_zsign(ain);
	if (sgn < 0) {
		zhalt("InvMod: first input negative");
	}

	if (_ntl_zcompare(ain, nin) >= 0) {
		zhalt("InvMod: first input too big");
	}


	if (sgn == 0) {
		_ntl_zcopy(nin, invv);
		return 1;
	}

#ifdef NTL_GMP_HACK
#if (__GNU_MP_VERSION >= 3)
/* versions of gmp prior to version 3 have really bad xgcd */
	if (_ntl_gmp_hack && nin[0] >= NTL_GMP_INVMOD_CROSS) {
		GT_INIT

		lip_to_mpz(ain, gt_1);
		lip_to_mpz(nin, gt_2);
		mpz_gcdext(gt_3, gt_4, 0, gt_1, gt_2);
		mpz_to_lip(&u, gt_3);
		mpz_to_lip(&v, gt_4);

		if (u && u[0] == 1 && u[1] == 1) {

			/* 
			 * We make sure v is in range 0..n-1,
			 * just in case GMP is sloppy.
		         */

                        _ntl_zmod(v, nin, &v);
			_ntl_zcopy(v, invv);
			return 0;
		}
		else {
			_ntl_zcopy(u, invv);
			return 1;
		}
	}
#endif
#endif

	if (!(zxxeucl(ain, nin, &v, &u))) {
		if (_ntl_zsign(v) < 0) _ntl_zadd(v, nin, &v);
		_ntl_zcopy(v, invv);
		return 0;
	}

	_ntl_zcopy(u, invv);
	return 1;
}

void
_ntl_zexteucl(
	_ntl_verylong aa,
	_ntl_verylong *xa,
	_ntl_verylong bb,
	_ntl_verylong *xb,
	_ntl_verylong *d
	)
{
	static _ntl_verylong modcon = 0;
	static _ntl_verylong a=0, b=0;
	long anegative = 0;
	long bnegative = 0;

	_ntl_zcopy(aa, &a);
	_ntl_zcopy(bb, &b);

	if (anegative = (a[0] < 0))
		a[0] = -a[0];
	if (bnegative = (b[0] < 0))
		b[0] = -b[0];

	if (!b[1] && (b[0] == 1))
	{
		_ntl_zone(xa);
		_ntl_zzero(xb);
		_ntl_zcopy(a, d);
		goto done;
	}

	if (!a[1] && (a[0] == 1))
	{
		_ntl_zzero(xa);
		_ntl_zone(xb);
		_ntl_zcopy(b, d);
		goto done;
	}

#ifdef NTL_GMP_HACK
#if (__GNU_MP_VERSION >= 3)
/* versions of gmp prior to version 3 have really bad xgcd */
   if (_ntl_gmp_hack && 
       a[0] >= NTL_GMP_XGCD_CROSS && b[0] >= NTL_GMP_XGCD_CROSS) {

      static _ntl_verylong xa1 = 0, xb1 = 0, d1 = 0, b_red = 0;
      static _ntl_verylong tmp1 = 0, tmp2 = 0, tmp3 = 0;

      GT_INIT

      

      lip_to_mpz(aa, gt_1);
      lip_to_mpz(bb, gt_2);
      mpz_gcdext(gt_3, gt_4, gt_5, gt_1, gt_2);
      mpz_to_lip(&d1, gt_3);
      mpz_to_lip(&xa1, gt_4);
      mpz_to_lip(&xb1, gt_5);

      /* normalize...this ensures results agree with
         classical Euclid...not very efficient... */

      if (_ntl_zcompare(a, b) >= 0) {
         _ntl_zdiv(a, d1, &a, 0);
         _ntl_zdiv(b, d1, &b, 0);
         _ntl_zdiv(xa1, b, &tmp1, &xa1);
         _ntl_zlshift(xa1, 1, &tmp2);
         if (anegative) _ntl_zsadd(tmp2, 1, &tmp2);
         if (_ntl_zcompare(tmp2, b) > 0) {
            _ntl_zsub(xa1, b, &xa1);
            _ntl_zsadd(tmp1, 1, &tmp1);
         }
         _ntl_zmul(tmp1, a, &tmp2);
         if (anegative != bnegative) 
            _ntl_zsub(xb1, tmp2, &xb1);
         else
            _ntl_zadd(xb1, tmp2, &xb1);
      }
      else {
         _ntl_zdiv(a, d1, &a, 0);
         _ntl_zdiv(b, d1, &b, 0);
         _ntl_zdiv(xb1, a, &tmp1, &xb1);
         _ntl_zlshift(xb1, 1, &tmp2);
         if (bnegative) _ntl_zsadd(tmp2, 1, &tmp2);
         if (_ntl_zcompare(tmp2, a) > 0) {
            _ntl_zsub(xb1, a, &xb1);
            _ntl_zsadd(tmp1, 1, &tmp1);
         }
         _ntl_zmul(tmp1, b, &tmp2);
         if (anegative != bnegative) 
            _ntl_zsub(xa1, tmp2, &xa1);
         else
            _ntl_zadd(xa1, tmp2, &xa1);
      }


      /* end normalize */

      _ntl_zcopy(d1, d);
      _ntl_zcopy(xa1, xa);
      _ntl_zcopy(xb1, xb);
      return;
   }
#endif
#endif


	zxxeucl(a, b, xa, d);
	_ntl_zmul(a, *xa, xb);
	_ntl_zsub(*d, *xb, xb);
	_ntl_zdiv(*xb, b, xb, &modcon);

	if ((modcon[1]) || (modcon[0] != 1))
	{
		zhalt("non-zero remainder in _ntl_zexteucl   BUG");
	}
done:
	if (anegative)
	{
		_ntl_znegate(xa);
	}
	if (bnegative)
	{
		_ntl_znegate(xb);
	}
}



/* I've adapted LIP's extended euclidean algorithm to
 * do rational reconstruction.  -- VJS.
 */


long 
_ntl_zxxratrecon(
   _ntl_verylong ain,
   _ntl_verylong nin,
   _ntl_verylong num_bound,
   _ntl_verylong den_bound,
   _ntl_verylong *num_out,
   _ntl_verylong *den_out
   )
{
   static _ntl_verylong a = 0;
   static _ntl_verylong n = 0;
   static _ntl_verylong q = 0;
   static _ntl_verylong w = 0;
   static _ntl_verylong x = 0;
   static _ntl_verylong y = 0;
   static _ntl_verylong z = 0;
   static _ntl_verylong inv = 0;
   static _ntl_verylong u = 0;
   static _ntl_verylong a_bak = 0;
   static _ntl_verylong n_bak = 0;
   static _ntl_verylong inv_bak = 0;
   static _ntl_verylong w_bak = 0;

   _ntl_verylong p;

   long diff;
   long ilo;
   long sa;
   long sn;
   long snum;
   long sden;
   long e;
   long fast;
   long temp;
   long parity;
   long gotthem;
   long try11;
   long try12;
   long try21;
   long try22;
   long got11;
   long got12;
   long got21;
   long got22;

   double hi;
   double lo;
   double dt;
   double fhi, fhi1;
   double flo, flo1;
   double num;
   double den;
   double dirt;

   if (_ntl_zsign(num_bound) < 0)
      zhalt("rational reconstruction: bad numerator bound");

   if (!num_bound)
      snum = 1;
   else
      snum = num_bound[0];

   if (_ntl_zsign(den_bound) <= 0)
      zhalt("rational reconstruction: bad denominator bound");

   sden = den_bound[0];

   if (_ntl_zsign(nin) <= 0)
      zhalt("rational reconstruction: bad modulus");

   if (_ntl_zsign(ain) < 0 || _ntl_zcompare(ain, nin) >= 0)
      zhalt("rational reconstruction: bad residue");

      
   e = nin[0];
   _ntl_zsetlength(&a, e);
   _ntl_zsetlength(&n, e);
   _ntl_zsetlength(&q, e);
   _ntl_zsetlength(&w, e);
   _ntl_zsetlength(&x, e);
   _ntl_zsetlength(&y, e);
   _ntl_zsetlength(&z, e);
   _ntl_zsetlength(&inv, e);
   _ntl_zsetlength(&u, e);
   _ntl_zsetlength(&a_bak, e);
   _ntl_zsetlength(&n_bak, e);
   _ntl_zsetlength(&inv_bak, e);
   _ntl_zsetlength(&w_bak, e);

   fhi1 = 1.0 + ((double) 32.0)/NTL_FDOUBLE_PRECISION;
   flo1 = 1.0 - ((double) 32.0)/NTL_FDOUBLE_PRECISION;

   fhi = 1.0 + ((double) 8.0)/NTL_FDOUBLE_PRECISION;
   flo = 1.0 - ((double) 8.0)/NTL_FDOUBLE_PRECISION;

   _ntl_zcopy(ain, &a);
   _ntl_zcopy(nin, &n);

   _ntl_zone(&inv);
   _ntl_zzero(&w);

   while (1)
   {
      if (w[0] >= sden && _ntl_zcompare(w, den_bound) > 0) break;
      if (n[0] <= snum && _ntl_zcompare(n, num_bound) <= 0) break;

      _ntl_zcopy(a, &a_bak);
      _ntl_zcopy(n, &n_bak);
      _ntl_zcopy(w, &w_bak);
      _ntl_zcopy(inv, &inv_bak);

      gotthem = 0;
      sa = a[0];
      sn = n[0];
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = a[0];
         p = &a[sa];
         num = (double) (*p) * NTL_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_FRADIX;
         if (sa > 2)
            num += (*(p - 1));
         sn = n[0];
         p = &n[sn];
         den = (double) (*p) * NTL_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_FRADIX;
         if (sn > 2)
            den += (*(p - 1));
         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_FRADIX;
            lo *= NTL_FRADIX;
         }
         try11 = 1;
         try12 = 0;
         try21 = 0;
         try22 = 1;
         parity = 1;
         fast = 1; 
         while (fast > 0)
         {
            parity = 1 - parity;
            if (hi >= NTL_FRADIX)
               fast = 0;
            else
            {
               ilo = (long)lo;
               dirt = hi - ilo;
               if (dirt < 1.0/NTL_FDOUBLE_PRECISION || !ilo || ilo < (long)hi)
                  fast = 0;
               else
               {
                  dt = lo-ilo;
                  lo = flo / dirt;
                  if (dt > 1.0/NTL_FDOUBLE_PRECISION)
                     hi = fhi / dt;
                  else
                     hi = NTL_FRADIX;
                  temp = try11;
                  try11 = try21;
                  if ((NTL_RADIX - temp) / ilo < try21)
                     fast = 0;
                  else
                     try21 = temp + ilo * try21;
                  temp = try12;
                  try12 = try22;
                  if ((NTL_RADIX - temp) / ilo < try22)
                     fast = 0;
                  else
                     try22 = temp + ilo * try22;
                  if ((fast > 0) && (parity > 0))
                  {
                     gotthem = 1;
                     got11 = try11;
                     got12 = try12;
                     got21 = try21;
                     got22 = try22;
                  }
               }
            }
         }
      }
      if (gotthem)
      {
         _ntl_zsmul(inv, got11, &x);
         _ntl_zsmul(w, got12, &y);
         _ntl_zsmul(inv, got21, &z);
         _ntl_zsmul(w, got22, &w);
         _ntl_zadd(x, y, &inv);
         _ntl_zadd(z, w, &w);
         _ntl_zsmul(a, got11, &x);
         _ntl_zsmul(n, got12, &y);
         _ntl_zsmul(a, got21, &z);
         _ntl_zsmul(n, got22, &n);
         _ntl_zsub(x, y, &a);
         _ntl_zsub(n, z, &n);
      }
      else
      {
         _ntl_zdiv(a, n, &q, &a);
         _ntl_zmul(q, w, &x);
         _ntl_zadd(inv, x, &inv);
         if (a[0] > 1 || a[1] > 0)
         {
            _ntl_zdiv(n, a, &q, &n);
            _ntl_zmul(q, inv, &x);
            _ntl_zadd(w, x, &w);
         }
         else
         {
            break;
         }
      }
   }

   _ntl_zcopy(a_bak, &a);
   _ntl_zcopy(n_bak, &n);
   _ntl_zcopy(w_bak, &w);
   _ntl_zcopy(inv_bak, &inv);

   _ntl_znegate(&w);

   while (1)
   {
      sa = w[0];
      if (sa < 0) w[0] = -sa;
      if (w[0] >= sden && _ntl_zcompare(w, den_bound) > 0) return 0;
      w[0] = sa;

      if (n[0] <= snum && _ntl_zcompare(n, num_bound) <= 0) break;
      
      fast = 0;
      sa = a[0];
      sn = n[0];
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = a[0];
         p = &a[sa];
         num = (double) (*p) * NTL_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_FRADIX;
         if (sa > 2)
            num += (*(p - 1));
         sn = n[0];
         p = &n[sn];
         den = (double) (*p) * NTL_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_FRADIX;
         if (sn > 2)
            den += (*(p - 1));
         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_FRADIX;
            lo *= NTL_FRADIX;
         }
         if (hi < NTL_FRADIX)
         {
            ilo = (long)lo;
            if (ilo == (long)hi)
               fast = 1;
         }
      }

      if (fast) 
      {
         if (ilo != 0) {
            if (ilo == 1) {
               _ntl_zsub(inv, w, &inv);
               _ntl_zsubpos(a, n, &a);
            }
            else if (ilo == 2) {
               _ntl_z2mul(w, &x);
               _ntl_zsub(inv, x, &inv);
               _ntl_z2mul(n, &x);
               _ntl_zsubpos(a, x, &a);
            }
            else if (ilo ==3) {
               _ntl_z2mul(w, &x);
               _ntl_zadd(w, x, &x);
               _ntl_zsub(inv, x, &inv);
               _ntl_z2mul(n, &x);
               _ntl_zadd(n, x, &x);
               _ntl_zsubpos(a, x, &a);
            }
            else if (ilo == 4) {
               _ntl_zlshift(w, 2, &x);
               _ntl_zsub(inv, x, &inv);
               _ntl_zlshift(n, 2, &x);
               _ntl_zsubpos(a, x, &a);
            }
            else {
               _ntl_zsmul(w, ilo, &x);
               _ntl_zsub(inv, x, &inv);
               _ntl_zsmul(n, ilo, &x);
               _ntl_zsubpos(a, x, &a);
            }
         }
      }
      else {
         _ntl_zdiv(a, n, &q, &a);
         _ntl_zmul(q, w, &x);
         _ntl_zsub(inv, x, &inv);
      }

      _ntl_zswap(&a, &n);
      _ntl_zswap(&inv, &w);
   }

   if (_ntl_zsign(w) < 0) {
      _ntl_znegate(&w);
      _ntl_znegate(&n);
   }

   _ntl_zcopy(n, num_out);
   _ntl_zcopy(w, den_out);

   return 1;
}



static
long OptWinSize(long n)
/* finds k that minimizes n/(k+1) + 2^{k-1} */

{
   long k;
   double v, v_new;


   v = n/2.0 + 1.0;
   k = 1;

   for (;;) {
      v_new = n/((double)(k+2)) + ((double)(1L << k));
      if (v_new >= v) break;
      v = v_new;
      k++;
   }

   return k;
}

      

static
void _ntl_zsppowermod(long a, _ntl_verylong e, _ntl_verylong n, 
                      _ntl_verylong *x)
{
   _ntl_verylong res;
   long i, k;

   if (_ntl_ziszero(e)) {
      _ntl_zone(x);
      return;
   }

   res = 0;
   _ntl_zsetlength(&res, n[0]);
   _ntl_zone(&res);

   k = _ntl_z2log(e);

   for (i = k - 1; i >= 0; i--) {
      _ntl_zsqmod(res, n, &res);
      if (_ntl_zbit(e, i))
         _ntl_zsmulmod(res, a, n, &res);
   }

   if (_ntl_zsign(e) < 0) _ntl_zinvmod(res, n, &res);

   _ntl_zcopy(res, x);
   _ntl_zfree(&res);
}



static
void _ntl_ztwopowermod( _ntl_verylong e, _ntl_verylong n, 
                      _ntl_verylong *x)
{
   _ntl_verylong res;
   long i, k;

   if (_ntl_ziszero(e)) {
      _ntl_zone(x);
      return;
   }

   res = 0;
   _ntl_zsetlength(&res, n[0]);
   _ntl_zone(&res);

   k = _ntl_z2log(e);

   for (i = k - 1; i >= 0; i--) {
      _ntl_zsqmod(res, n, &res);
      if (_ntl_zbit(e, i))
         _ntl_zaddmod(res, res, n, &res);
   }

   if (_ntl_zsign(e) < 0) _ntl_zinvmod(res, n, &res);

   _ntl_zcopy(res, x);
   _ntl_zfree(&res);
}


void _ntl_zpowermod(_ntl_verylong g, _ntl_verylong e, _ntl_verylong F,
                    _ntl_verylong *h)

/* h = g^e mod f using "sliding window" algorithm

   remark: the notation (h, g, e, F) is strange, because I
   copied the code from BB.c.
*/

{
   _ntl_verylong res, *v, t;
   long n, i, k, val, cnt, m;

   if (_ntl_zsign(g) < 0 || _ntl_zcompare(g, F) >= 0 || 
       _ntl_zscompare(F, 1) <= 0) 
      zhalt("PowerMod: bad args");

#ifdef NTL_GMP_HACK

   if (_ntl_gmp_hack) {
      res = 0;
      _ntl_gmp_powermod(&res, g, e, F);

      /* careful...don't compute directly into h, as it could lead
       * to an allocation error in some cases.
       */

      _ntl_zcopy(res, h);
      _ntl_zfree(&res);
      return;
   }

#endif

   if (!g || g[0] == 1 || g[0] == -1) {
      long gg = _ntl_ztoint(g);
      if (gg == 2)
         _ntl_ztwopowermod(e, F, h);
      else
         _ntl_zsppowermod(gg, e, F, h);
      return;
   }

   if (_ntl_zscompare(e, 0) == 0) {
      _ntl_zone(h);
      return;
   }

   if (_ntl_zscompare(e, 1) == 0) {
      _ntl_zcopy(g, h);
      return;
   }

   if (_ntl_zscompare(e, -1) == 0) {
      _ntl_zinvmod(g, F, h);
      return;
   }

   if (_ntl_zscompare(e, 2) == 0) {
      _ntl_zsqmod(g, F, h);
      return;
   }

   if (_ntl_zscompare(e, -2) == 0) {
      res = 0;
      _ntl_zsqmod(g, F, &res);
      _ntl_zinvmod(res, F, h);
      _ntl_zfree(&res);
      return;
   }

   n = _ntl_z2log(e);

   res = 0;
   _ntl_zone(&res);

   if (n < 16) {
      /* plain square-and-multiply algorithm */

      for (i = n - 1; i >= 0; i--) {
         _ntl_zsqmod(res, F, &res);
         if (_ntl_zbit(e, i))
            _ntl_zmulmod(res, g, F, &res);
      }

      if (_ntl_zsign(e) < 0) _ntl_zinvmod(res, F, &res);

      _ntl_zcopy(res, h);
      _ntl_zfree(&res);
      return;
   }

   k = OptWinSize(n);

   if (k > 5) k = 5;

   v = (_ntl_verylong *) NTL_MALLOC(1L << (k-1), sizeof(_ntl_verylong), 0);
   if (!v) zhalt("out of memory");
   for (i = 0; i < (1L << (k-1)); i++)
      v[i] = 0; 

   _ntl_zcopy(g, &v[0]);
 
   if (k > 1) {
      t = 0;
      _ntl_zsqmod(g, F, &t);

      for (i = 1; i < (1L << (k-1)); i++)
         _ntl_zmulmod(v[i-1], t, F, &v[i]);

      _ntl_zfree(&t);
   }


   val = 0;
   for (i = n-1; i >= 0; i--) {
      val = (val << 1) | _ntl_zbit(e, i); 
      if (val == 0)
         _ntl_zsqmod(res, F, &res);
      else if (val >= (1L << (k-1)) || i == 0) {
         cnt = 0;
         while ((val & 1) == 0) {
            val = val >> 1;
            cnt++;
         }

         m = val;
         while (m > 0) {
            _ntl_zsqmod(res, F, &res);
            m = m >> 1;
         }

         _ntl_zmulmod(res, v[val >> 1], F, &res);

         while (cnt > 0) {
            _ntl_zsqmod(res, F, &res);
            cnt--;
         }

         val = 0;
      }
   }

   if (_ntl_zsign(e) < 0) _ntl_zinvmod(res, F, &res);

   _ntl_zcopy(res, h);

   _ntl_zfree(&res);
   for (i = 0; i < (1L << (k-1)); i++)
      _ntl_zfree(&v[i]);
   free(v);
}






void
_ntl_zexp(
	_ntl_verylong a,
	long e,
	_ntl_verylong *bb
	)
{
	long k;
	long len_a;
	long sa;
	static _ntl_verylong res = 0;

	if (!a)
		sa = 0;
	else {
		sa = a[0];
		if (sa < 0) sa = -sa;
	}

	if (sa <= 1) {
		_ntl_zexps(_ntl_ztoint(a), e, bb);
		return;
	}


	if (!e)
	{
		_ntl_zone(bb);
		return;
	}

	if (e < 0)
		zhalt("negative exponent in _ntl_zexp");

	if (_ntl_ziszero(a))
	{
		_ntl_zzero(bb);
		return;
	}

	len_a = _ntl_z2log(a);
	if (len_a > (NTL_MAX_LONG-(NTL_NBITS-1))/e)
		zhalt("overflow in _ntl_zexp");

	_ntl_zsetlength(&res, (len_a*e+NTL_NBITS-1)/NTL_NBITS);

	_ntl_zcopy(a, &res);
	k = 1;
	while ((k << 1) <= e)
		k <<= 1;
	while (k >>= 1) {
		_ntl_zsq(res, &res);
		if (e & k)
			_ntl_zmul(a, res, &res);
	}

	_ntl_zcopy(res, bb);
}

void
_ntl_zexps(
	long a,
	long e,
	_ntl_verylong *bb
	)
{
	long k;
	long len_a;
	static _ntl_verylong res = 0;

	if (!e)
	{
		_ntl_zone(bb);
		return;
	}

	if (e < 0)
		zhalt("negative exponent in _ntl_zexps");

	if (!a)
	{
		_ntl_zzero(bb);
		return;
	}

	if (a >= NTL_RADIX || a <= -NTL_RADIX) {
		_ntl_zintoz(a, &res);
		_ntl_zexp(res, e, &res);
		return;
	}

	len_a = _ntl_z2logs(a);
	if (len_a > (NTL_MAX_LONG-(NTL_NBITS-1))/e)
		zhalt("overflow in _ntl_zexps");

	_ntl_zsetlength(&res, (len_a*e+NTL_NBITS-1)/NTL_NBITS);

	_ntl_zintoz(a, &res);
	k = 1;
	while ((k << 1) <= e)
		k <<= 1;
	while (k >>= 1) {
		_ntl_zsq(res, &res);
		if (e & k)
			_ntl_zsmul(res, a, &res);
	}

	_ntl_zcopy(res, bb);
}


void
_ntl_z2mul(
	_ntl_verylong n,
	_ntl_verylong *rres
	)
{
	long sn;
	long i;
	long n_alias;
	long carry;
	_ntl_verylong res;

	if (!n)
	{
		_ntl_zzero(rres);
		return;
	}


	if ((!n[1]) && (n[0] == 1))
	{
		_ntl_zzero(rres);
		return;
	}

	if ((sn = n[0]) < 0)
		sn = -sn;

	res = *rres;
	n_alias = (n == res);

	_ntl_zsetlength(&res, sn + 1);
	if (n_alias) n = res;
	*rres = res;

	carry = 0;

	for (i = 1; i <= sn; i++)
	{
		if ((res[i] = (n[i] << 1) + carry) >= NTL_RADIX)
		{
			res[i] -= NTL_RADIX;
			carry = 1;
		}
		else
			carry = 0;
	}

	if (carry)
		res[++sn] = 1;

	if (n[0] < 0)
		res[0] = -sn;
	else
		res[0] = sn;
}


long 
_ntl_z2div(
	_ntl_verylong n,
	_ntl_verylong *rres
	)
{
	long sn;
	long i;
	long result;
	_ntl_verylong res = *rres;

	if ((!n) || ((!n[1]) && (n[0] == 1)))
	{
		_ntl_zzero(rres);
		return (0);
	}

	if ((sn = n[0]) < 0)
		sn = -sn;

	/* n won't move if res aliases n */
	_ntl_zsetlength(&res, sn);
	*rres = res;

	result = n[1] & 1;
	for (i = 1; i < sn; i++)
	{
		res[i] = (n[i] >> 1);
		if (n[i + 1] & 1)
			res[i] += (NTL_RADIX >> 1);
	}

	if (res[sn] = (n[sn] >> 1))
		res[0] = n[0];
	else if (sn == 1)
	{
		res[0] = 1;
	}
	else if (n[0] < 0)
		res[0] = -sn + 1;
	else
		res[0] = sn - 1;

	return (result);
}


void
_ntl_zlshift(
	_ntl_verylong n,
	long k,
	_ntl_verylong *rres
	)
{
	long big;
	long small;
	long sn;
	long i;
	long cosmall;
	long n_alias;
	_ntl_verylong res;


	if (!n)
	{
		_ntl_zzero(rres);
		return;
	}

	if ((!n[1]) && (n[0] == 1))
	{
		_ntl_zzero(rres);
		return;
	}

	res = *rres;
	n_alias = (n == res);
	

	if (!k)
	{
		if (!n_alias)
			_ntl_zcopy(n, rres);
		return;
	}

	if (k < 0)
	{
		if (k < -NTL_MAX_LONG) 
			_ntl_zzero(rres); 
		else
			_ntl_zrshift(n, -k, rres);
		return;
	}
	if (k == 1)
	{
		_ntl_z2mul(n, rres);
		return;
	}

	if ((sn = n[0]) < 0)
		sn = -sn;

	i = sn + (big = k / NTL_NBITS);
	if (small = k - big * NTL_NBITS)
	{
		_ntl_zsetlength(&res, i + 1);
		if (n_alias) n = res;
		*rres = res;

		res[i + 1] = n[sn] >> (cosmall = NTL_NBITS - small);
		for (i = sn; i > 1; i--)
			res[i + big] = ((((unsigned long) n[i]) << small) & NTL_RADIXM) + (n[i - 1] >> cosmall);
		res[big + 1] = (((unsigned long) n[1]) << small) & NTL_RADIXM;
		for (i = big; i; i--)
			res[i] = 0;
		if (res[sn + big + 1])
			big++;
	}
	else
	{
		_ntl_zsetlength(&res, i);
		if (n_alias) n = res;
		*rres = res;

		for (i = sn; i; i--)
			res[i + big] = n[i];
		for (i = big; i; i--)
			res[i] = 0;
	}
	if (n[0] > 0)
		res[0] = n[0] + big;
	else
		res[0] = n[0] - big;
}


void
_ntl_zrshift(
	_ntl_verylong n,
	long k,
	_ntl_verylong *rres
	)
{
	long big;
	long small;
	long sn;
	long i;
	long cosmall;
	_ntl_verylong res;

	if (!n)
	{
		_ntl_zzero(rres);
		return;
	}

	if ((!n[1]) && (n[0] == 1))
	{
		_ntl_zzero(rres);
		return;
	}

	res = *rres;

	if (!k)
	{
		if (n != res)
			_ntl_zcopy(n, rres);
		return;
	}

	if (k < 0)
	{
		if (k < -NTL_MAX_LONG) zhalt("overflow in _ntl_zrshift");
		_ntl_zlshift(n, -k, rres);
		return;
	}

	if (k == 1)
	{
		_ntl_z2div(n, rres);
		return;
	}

	big = k / NTL_NBITS;
	small = k - big * NTL_NBITS;

	if ((sn = n[0]) < 0)
		sn = -sn;

	if ((big >= sn) || 
            ((big == sn - 1) && small && (!(n[sn] >> small))))
        /* The microsoft optimizer generates bad code without
           the above test for small != 0 */
	{
		_ntl_zzero(rres);
		return;
	}

	sn -= big;

	/* n won't move if res aliases n */
	_ntl_zsetlength(&res, sn);
	*rres = res;

	if (small)
	{
		cosmall = NTL_NBITS - small;
		for (i = 1; i < sn; i++)
			res[i] = (n[i + big] >> small) +
				((((unsigned long) n[i + big + 1]) << cosmall) & NTL_RADIXM);
		if (!(res[sn] = (n[sn + big] >> small)))
			sn--;
	}
	else
		for (i = 1; i <= sn; i++)
			res[i] = n[i + big];
	if (n[0] > 0)
		res[0] = sn;
	else
		res[0] = -sn;
}


long 
_ntl_zmakeodd(
	_ntl_verylong *nn
	)
{
	_ntl_verylong n = *nn;
	long i;
	long shift = 1;

	if (!n || (!n[1] && (n[0] == 1)))
		return (0);
	while (!(n[shift]))
		shift++;
	i = n[shift];
	shift = NTL_NBITS * (shift - 1);
	while (!(i & 1))
	{
		shift++;
		i >>= 1;
	}
	_ntl_zrshift(n, shift, &n);
	return (shift);
}


long 
_ntl_znumtwos(
	_ntl_verylong n
	)
{
	long i;
	long shift = 1;

	if (!n || (!n[1] && (n[0] == 1)))
		return (0);
	while (!(n[shift]))
		shift++;
	i = n[shift];
	shift = NTL_NBITS * (shift - 1);
	while (!(i & 1))
	{
		shift++;
		i >>= 1;
	}
	return (shift);
}


long 
_ntl_zsqrts(
	long n
	)
{
	long a;
	long ndiva;
	long newa;
	static _ntl_verylong ln=0;
	static _ntl_verylong rr=0;

	if (n < 0) 
		zhalt("_ntl_zsqrts: negative argument");

	if (n <= 0)
		return (0);
	if (n <= 3)
		return (1);
	if (n <= 8)
		return (2);
	if (n >= NTL_RADIX)
	{
		_ntl_zintoz(n,&ln);
		_ntl_zsqrt(ln,&rr);
		return(_ntl_ztoint(rr));
	}
	newa = 3L << (2 * (NTL_NBITSH - 1));
	a = 1L << NTL_NBITSH;
	while (!(n & newa))
	{
		newa >>= 2;
		a >>= 1;
	}
	while (1)
	{
		newa = ((ndiva = n / a) + a) / 2;
		if (newa - ndiva <= 1)
		{
			if (newa * newa <= n)
				return (newa);
			else
				return (ndiva);
		}
		a = newa;
	}
}


void _ntl_zsqrt(_ntl_verylong n, _ntl_verylong *rr)
{
	static _ntl_verylong a = 0;
	static _ntl_verylong ndiva = 0;
	static _ntl_verylong diff = 0;
	static _ntl_verylong r = 0;
	long i;

	if (!n) {
		_ntl_zzero(rr);
		return;
	}

	if ((i = n[0]) < 0)
		zhalt("negative argument in _ntl_zsqrt");

	if (i == 1) {
		_ntl_zintoz(_ntl_zsqrts(n[1]), rr);
		return;
	}

#ifdef NTL_GMP_HACK
   if (_ntl_gmp_hack && i >= NTL_GMP_SQRT_CROSS) {
      GT_INIT

      lip_to_mpz(n, gt_1);
      mpz_sqrt(gt_2, gt_1);
      mpz_to_lip(&r, gt_2);
      _ntl_zcopy(r, rr);
      return;
   }
#endif

	_ntl_zsetlength(&a, i);
	_ntl_zsetlength(&ndiva, i);
	_ntl_zsetlength(&diff, i);

	a[(a[0] = (i + 1) / 2)] = _ntl_zsqrts(n[i]) + 1;
	if (!(i & 1))
		a[a[0]] <<= NTL_NBITSH;

	if (a[a[0]] & NTL_RADIX) {
		a[a[0]] = 0;
		a[0]++;
		a[a[0]] = 1;
	}

	for (i = a[0] - 1; i; i--)
		a[i] = 0;

	while (1) {
		_ntl_zdiv(n, a, &ndiva, &r);
		_ntl_zadd(a, ndiva, &r);
		_ntl_zrshift(r, 1, &r);
		if (_ntl_zcompare(r, ndiva) <= 0) 
			goto done;

		_ntl_zsubpos(r, ndiva, &diff);
		if ((diff[0] == 1) && (diff[1] <= 1)) {
			_ntl_zsq(r, &diff);
			if (_ntl_zcompare(diff, n) > 0)
				_ntl_zcopy(ndiva, &r);

			goto done;
		}
		_ntl_zcopy(r, &a);
	}
done:
	_ntl_zcopy(r, rr);
}



void
_ntl_zgcd(
	_ntl_verylong mm1,
	_ntl_verylong mm2,
	_ntl_verylong *rres
	)
{
	long agrb;
	long shibl;
	static _ntl_verylong aa = 0;
	static _ntl_verylong bb = 0;
	static _ntl_verylong cc = 0;
	_ntl_verylong a;
	_ntl_verylong b;
	_ntl_verylong c;
	_ntl_verylong d;
	long m1negative;
	long m2negative;

	/* _ntl_ziszero is necessary here and below to fix an
	   an aliasing bug in LIP */

	if (_ntl_ziszero(mm1))
	{
		if (mm2 != *rres)
			_ntl_zcopy(mm2,rres);
		_ntl_zabs(rres);
		return;
	}

	if (_ntl_ziszero(mm2))
	{
		if (mm1 != *rres)
			_ntl_zcopy(mm1,rres);
		_ntl_zabs(rres);
		return;
	}

	if (mm1 == mm2)
	{
		if (mm1 != *rres)
			_ntl_zcopy(mm1, rres);
		_ntl_zabs(rres);
		return;
	}

#ifdef NTL_GMP_HACK
   {
      long s1 = mm1[0];
      long s2 = mm2[0];

      if (s1 < 0) s1 = -s1;
      if (s2 < 0) s2 = -s2;

      if (_ntl_gmp_hack && s1 >= NTL_GMP_GCD_CROSS && s2 >= NTL_GMP_GCD_CROSS) {
         GT_INIT 

         lip_to_mpz(mm1, gt_1);
         lip_to_mpz(mm2, gt_2);
         mpz_gcd(gt_3, gt_1, gt_2);
         mpz_to_lip(&aa, gt_3);
         _ntl_zcopy(aa, rres);
         return;
      }
   }
#endif

	if (m1negative = (mm1[0] < 0))
		mm1[0] = -mm1[0];
	if (m2negative = (mm2[0] < 0))
		mm2[0] = -mm2[0];

	if ((agrb = mm1[0]) < mm2[0])
		agrb = mm2[0];

	_ntl_zsetlength(&aa, agrb+1); 
	_ntl_zsetlength(&bb, agrb+1);
	_ntl_zsetlength(&cc, agrb+1);

	if (mm1[0] != mm2[0])
	{
		if (mm1[0] > mm2[0])
		{
			_ntl_zcopy(mm2, &aa);
			_ntl_zmod(mm1, aa, &bb);
		}
		else
		{
			_ntl_zcopy(mm1, &aa);
			_ntl_zmod(mm2, aa, &bb);
		}
		if (!(bb[1]) && (bb[0] == 1))
		{
			a = aa;
			goto done;
		}
	}
	else
	{
		_ntl_zcopy(mm1, &aa);
		_ntl_zcopy(mm2, &bb);
	}
	if ((agrb = _ntl_zmakeodd(&aa)) < (shibl = _ntl_zmakeodd(&bb)))
		shibl = agrb;
	if (!(agrb = _ntl_zcompare(aa, bb)))
	{
		a = aa;
		goto endshift;
	}
	else if (agrb < 0)
	{
		a = bb;
		b = aa;
	}
	else
	{
		a = aa;
		b = bb;
	}
	c = cc;
	_ntl_zsubpos(a, b, &c);
	do
	{
		_ntl_zmakeodd(&c);
		if (!(agrb = _ntl_zcompare(b, c)))
		{
			a = b;
			goto endshift;
		}
		else if (agrb > 0)
		{
			a = b;
			b = c;
			c = a;
		}
		else
		{
			d = a;
			a = c;
			c = d;
		}
		_ntl_zsubpos(a, b, &c);
	} while (c[1] || c[0] != 1);
endshift:
	_ntl_zlshift(a, shibl, &a);
done:
	if (m1negative)
		mm1[0] = -mm1[0];
	if (m2negative)
		mm2[0] = -mm2[0];
	_ntl_zcopy(a, rres);
}


long _ntl_zsign(_ntl_verylong a)
{
	if (!a)
	{
		return (0);
	}
	if (a[0] < 0)
		return (-1);
	if (a[0] > 1)
		return (1);
	if (a[1])
		return (1);
	return (0);
}

void _ntl_zabs(_ntl_verylong *pa)
{
	_ntl_verylong a = *pa;

	if (!a)
		return;
	if (a[0] < 0)
		a[0] = (-a[0]);
}

long 
_ntl_z2logs(
	long aa
	)
{
	long i = 0;
	unsigned long a;

	if (aa < 0)
		a = - ((unsigned long) aa);
	else
		a = aa;

	while (a>=256)
		i += 8, a >>= 8;
	if (a >=16)
		i += 4, a >>= 4;
	if (a >= 4)
		i += 2, a >>= 2;
	if (a >= 2)
		i += 2;
	else if (a >= 1)
		i++;
	return (i);
}

long 
_ntl_z2log(
	_ntl_verylong a
	)
{
	long la;

	if (!a)
		return (0);
	la = (a[0] > 0 ? a[0] : -a[0]);
	return ( NTL_NBITS * (la - 1) + _ntl_z2logs(a[la]) );
}



long
_ntl_zscompare(
	_ntl_verylong a,
	long b
	)
{
	if (!b) 
		return _ntl_zsign(a);
	else {
		static _ntl_verylong c = 0;
		_ntl_zintoz(b, &c);
		return (_ntl_zcompare(a, c));
	}
}

void
_ntl_zswap(
	_ntl_verylong *a,
	_ntl_verylong *b
	)
{
	_ntl_verylong c;

	if ((*a && ((*a)[-1] & 1)) || (*b && ((*b)[-1] & 1))) {
		static _ntl_verylong t = 0; 
		_ntl_zcopy(*a, &t);
		_ntl_zcopy(*b, a);
		_ntl_zcopy(t, b);
		return;
	}

	c = *a;
	*a = *b;
	*b = c;
}

long
_ntl_ziszero(
	_ntl_verylong a
	)
{
	if (!a) return (1);
	if (a[1]) return (0);
	if (a[0]==1) return (1);
	return (0);
}

long
_ntl_zodd(
	_ntl_verylong a
	)
{
	if (!a) return (0);
	return (a[1]&1);
}


long
_ntl_zbit(
	_ntl_verylong a,
	long p
	)
{
        long bl;
        long wh;
        long sa;

	if (p < 0 || !a) return 0;
	bl = (p/NTL_NBITS);
        wh = 1L << (p - NTL_NBITS*bl);
	bl ++;
        sa = a[0];
        if (sa < 0) sa = -sa;
        if (sa < bl) return (0);
	if (a[bl] & wh) return (1);
	return (0);
}


void
_ntl_zlowbits(
	_ntl_verylong a,
	long b,
        _ntl_verylong *cc
        )
{
        _ntl_verylong c;
	long bl;
	long wh;
	long sa;

	if (_ntl_ziszero(a) || (b<=0)) {
		_ntl_zzero(cc);
		return;
	}


        bl = b/NTL_NBITS;
        wh = b - NTL_NBITS*bl;
	if (wh != 0) 
		bl++;
	else
		wh = NTL_NBITS;

	sa = a[0];
	if (sa < 0) sa = -sa;

	if (sa < bl) {
		_ntl_zcopy(a,cc);
		_ntl_zabs(cc);
		return;
	}

	c = *cc;

	/* a won't move if c aliases a */
        _ntl_zsetlength(&c, bl);
        *cc = c;

	for (sa=1; sa<bl; sa++)
		c[sa] = a[sa];
	c[bl] = a[bl]&((1L<<wh)-1);
	while ((bl>1) && (!c[bl]))
		bl --;
	c[0] = bl;
}


long _ntl_zslowbits(_ntl_verylong a, long p)
{
   static _ntl_verylong x = 0;

   if (p > NTL_BITS_PER_LONG)
      p = NTL_BITS_PER_LONG;

   _ntl_zlowbits(a, p, &x);

   return _ntl_ztoint(x);
}




long
_ntl_zweights(
	long aa
	)
{
	unsigned long a;
	long res = 0;

	if (aa < 0) 
		a = - ((unsigned long) aa);
	else
		a = aa;
   
	while (a) {
		if (a & 1) res ++;
		a >>= 1;
	}
	return (res);
}

long
_ntl_zweight(
        _ntl_verylong a
        )
{
	long i;
	long res = 0;

	if (!a) return (0);
	i = a[0];
	if (i<0) i = -i;
	for (;i;i--)
		res += _ntl_zweights(a[i]);
	return (res);
}



void
_ntl_zand(
	_ntl_verylong a,
	_ntl_verylong b,
	_ntl_verylong *cc
	)

{
	_ntl_verylong c;
	long sa;
	long sb;
	long sm;
	long a_alias;
	long b_alias;

	if (_ntl_ziszero(a) || _ntl_ziszero(b)) {
		_ntl_zzero(cc);
		return;
	}

	c = *cc;
	a_alias = (a == c);
	b_alias = (b == c);

	sa = a[0];
	if (sa < 0) sa = -sa;

	sb = b[0];
	if (sb < 0) sb = -sb;

	sm = (sa > sb ? sb : sa );

	_ntl_zsetlength(&c, sm);
	if (a_alias) a = c;
	if (b_alias) b = c;
	*cc = c;

	for (sa = 1; sa <= sm; sa ++) 
		c[sa] = a[sa] & b[sa];

	while ((sm > 1) && (!(c[sm])))
		sm --;
	c[0] = sm;
}

void
_ntl_zxor(
	_ntl_verylong a,
	_ntl_verylong b,
	_ntl_verylong *cc
	)
{
        _ntl_verylong c;
        long sa;
        long sb;
        long sm;
        long la;
        long i;
	long a_alias;
	long b_alias;

	if (_ntl_ziszero(a)) {
		_ntl_zcopy(b,cc);
		_ntl_zabs(cc);
		return;
	}

	if (_ntl_ziszero(b)) {
		_ntl_zcopy(a,cc);
		_ntl_zabs(cc);
		return;
	}

	c = *cc;
	a_alias = (a == c);
	b_alias = (b == c);

        sa = a[0];
        if (sa < 0) sa = -sa;

        sb = b[0];
        if (sb < 0) sb = -sb;

	if (sa > sb) {
		la = sa;
		sm = sb;
	} else {
		la = sb;
		sm = sa;
	}

        _ntl_zsetlength(&c, la);
        if (a_alias) a = c;
        if (b_alias) b = c;
        *cc = c;

        for (i = 1; i <= sm; i ++)
                c[i] = a[i] ^ b[i];

	if (sa > sb)
		for (;i <= la; i++) c[i] = a[i];
	else
		for (;i <= la; i++) c[i] = b[i];

        while ((la > 1) && (!(c[la])))
                la --;
        c[0] = la;
}

void
_ntl_zor(
	_ntl_verylong a,
	_ntl_verylong b,
        _ntl_verylong *cc
        )
{
        _ntl_verylong c;
        long sa;
        long sb;
        long sm;
        long la;
        long i;
	long a_alias;
	long b_alias;

	if (_ntl_ziszero(a)) {
		_ntl_zcopy(b,cc);
		_ntl_zabs(cc);
		return;
	}

	if (_ntl_ziszero(b)) {
		_ntl_zcopy(a,cc);
		_ntl_zabs(cc);
		return;
	}

	c = *cc;
	a_alias = (a == c);
	b_alias = (b == c);

        sa = a[0];
        if (sa < 0) sa = -sa;

        sb = b[0];
        if (sb < 0) sb = -sb;

        if (sa > sb) {
                la = sa;
                sm = sb;
        } else {
                la = sb;
                sm = sa;
        }

        _ntl_zsetlength(&c, la);
        if (a_alias) a = c;
        if (b_alias) b = c;
        *cc = c;

        for (i = 1; i <= sm; i ++)
                c[i] = a[i] | b[i];

        if (sa > sb)
                for (;i <= la; i++) c[i] = a[i];
        else
                for (;i <= la; i++) c[i] = b[i];

        c[0] = la;
}

long
_ntl_zsetbit(
        _ntl_verylong *a,
        long b
        )
{
	long bl;
	long wh;
	long sa;

	if (b<0) zhalt("_ntl_zsetbit: negative index");

	if (_ntl_ziszero(*a)) {
		_ntl_zintoz(1,a);
		_ntl_zlshift(*a,b,a);
		return (0);
	}

	bl = (b/NTL_NBITS);
	wh = 1L << (b - NTL_NBITS*bl);
	bl ++;
	sa = (*a)[0];
	if (sa<0) sa = -sa;
	if (sa >= bl) {
		sa = (*a)[bl] & wh;
		(*a)[bl] |= wh;
		if (sa) return (1);
		return (0);
	} else {
		_ntl_zsetlength(a,bl);
		sa ++;
		for (;sa<=bl;sa++) (*a)[sa]=0;
		if ((*a)[0] < 0)
			(*a)[0] = -bl;
		else (*a)[0] = bl;
		(*a)[bl] |= wh;
		return (0);
	}
}

long
_ntl_zswitchbit(
        _ntl_verylong *a,
        long p
        )
{
        long bl;
        long wh;
        long sa;

        if (p < 0) zhalt("_ntl_zswitchbit: negative index");

        if (_ntl_ziszero(*a)) {
		_ntl_zintoz(1,a);
		_ntl_zlshift(*a,p,a);
		return (0);
	}

        bl = (p/NTL_NBITS);
        wh = 1L << (p - NTL_NBITS*bl);
        bl ++;
        sa = (*a)[0];
        if (sa < 0) sa = -sa;
        if ((sa < bl) || (!((*a)[bl] & wh))) {
		_ntl_zsetbit(a,p);
		return (0);
	}
	(*a)[bl] ^= wh;
        while ((sa>1) && (!(*a)[sa]))
		sa --;
	if ((*a)[0] > 0) (*a)[0] = sa;
	else (*a)[0] = -sa;
        return (1);
}


long _ntl_zsize(_ntl_verylong rep)
{
   if (!rep || (rep[0] == 1 && rep[1] == 0))
      return 0;
   else if (rep[0] < 0)
      return -rep[0];
   else
      return rep[0];
}


long _ntl_zdigit(_ntl_verylong rep, long i) 
{
   long sa;
   if (i < 0 || !rep) return 0;

   sa = rep[0];
   if (sa < 0) sa = -sa;
   if (i >= sa) return 0;
   return rep[i+1];
}

long _ntl_zisone(_ntl_verylong rep)
{
   return rep != 0 && rep[0] == 1 && rep[1] == 1;
}

long _ntl_zsptest(_ntl_verylong rep)
{
   return !rep || rep[0] == 1 || rep[0] == -1;
}

long _ntl_zwsptest(_ntl_verylong rep)
{
   return !rep || rep[0] == 1 || rep[0] == -1;
}

long _ntl_zcrtinrange(_ntl_verylong g, _ntl_verylong a)
{
   long sa, sg, carry, i, diff;

   if (!a || a[0] < 0 || (a[0] == 1 && a[1] == 0)) return 0;

   sa = a[0];

   if (!g) return 1;

   sg = g[0];

   if (sg == 1 && g[1] == 0) return 1;

   if (sg < 0) sg = -sg;

   if (sa-sg > 1) return 1;

   if (sa-sg < 0) return 0;

   carry=0;

   if (sa-sg == 1) {
      if (a[sa] > 1) return 1;
      carry = 1;
   }

   i = sg;
   diff = 0;
   while (i > 0 && diff == 0) {
      diff = (carry << (NTL_NBITS-1)) + (a[i] >> 1) - g[i];
      carry = (a[i] & 1);
      i--;
   }

   if (diff == 0) {
      if (carry) return 1;
      return (g[0] > 0);
   }
   else
      return (diff > 0);
}



void _ntl_zfrombytes(_ntl_verylong *x, const unsigned char *p, long n)
{
   long sz;
   long i;
   _ntl_verylong a;
   long bitpos, wordpos, bitoffset, diff;

   if (n <= 0) {
      _ntl_zzero(x);
      return;
   }

   if (n > (NTL_MAX_LONG-(NTL_NBITS-1))/8)
      zhalt("ZZFromBytes: excessive length");

   sz = (n*8 + NTL_NBITS-1)/NTL_NBITS;

   _ntl_zsetlength(x, sz);

   a = *x;

   for (i = 1; i <= sz; i++)
      a[i] = 0;

   for (i = 0; i < n; i++) {
      bitpos = i*8;
      wordpos = bitpos/NTL_NBITS;
      bitoffset = bitpos - wordpos*NTL_NBITS;
      diff = NTL_NBITS-bitoffset;

      if (diff < 8) {
         a[wordpos+1] |= 
            ((( ((unsigned long)(p[i])) & 255UL ) << bitoffset) & NTL_RADIXM);
         a[wordpos+2] = ( ((long)(p[i])) & 255 ) >> diff;
      }
      else {
         a[wordpos+1] |= (( ((long)(p[i])) & 255 ) << bitoffset);
      }
   }

   while (sz > 1 && a[sz] == 0) sz--;
   a[0] = sz;
}


void _ntl_zbytesfromz(unsigned char *p, _ntl_verylong a, long nn)
{
   long k = _ntl_z2log(a);
   long n = (k+7)/8;
   long sz = _ntl_zsize(a);
   long min_n =  ((n < nn) ? n : nn);

   long i;

   for (i = 0; i < min_n; i++) {
      long bitpos = i*8;
      long wordpos = bitpos/NTL_NBITS;
      long bitoffset = bitpos - wordpos*NTL_NBITS;
      long diff;

      p[i] = (a[wordpos+1] >> bitoffset) & 255;

      diff = NTL_NBITS - bitoffset;

      if (diff < 8 && wordpos < sz-1) {
         long msk = (1L << (8-diff))-1;
         p[i] |= ((a[wordpos+2] & msk) << diff);
      }
   }

   for (i = min_n; i < nn; i++)
      p[i] = 0;
}



long _ntl_zblock_construct_alloc(_ntl_verylong *x, long d, long n)
{
   long nwords, nbytes, AllocAmt, m, *p, *q, j;


   /* check n value */

   if (n <= 0)
      zhalt("block construct: n must be positive");

   /* check d value */

   if (d <= 0) 
      zhalt("block construct: d must be positive");

   if (NTL_OVERFLOW(d, NTL_NBITS, NTL_NBITS) || 
       NTL_OVERFLOW(d, sizeof(long), 3*sizeof(long)))
      zhalt("block construct: d too large");

#ifdef L_TO_G_CHECK_LEN
   /* this makes sure that numbers don't get too big for GMP */
   if (d >= (1L << (NTL_BITS_PER_INT-4)))
      zhalt("size too big for GMP");
#endif

   nwords = d + 3;
   nbytes = nwords*sizeof(long);
   
   AllocAmt = (NTL_MAX_ALLOC_BLOCK - sizeof(long)) / nbytes;
   if (AllocAmt == 0) AllocAmt = 1;

   if (AllocAmt < n)
      m = AllocAmt;
   else
      m = n;

   p = (long *) NTL_MALLOC(m, nbytes, sizeof(long));
   if (!p) zhalt("out of memory in block construct");

   *p = m;

   q = p+2;
   *x = q;
   
   for (j = 0; j < m; j++) {
      q[-1] = ((d+1) << 1) | 1;
      q[0] = 1;
      q[1] = 0;
      q += nwords;
   }

   return m;
}

void _ntl_zblock_construct_set(_ntl_verylong x, _ntl_verylong *y, long i)
{
   long d, size;

   d = (x[-1] >> 1) - 1;
   size = d + 3;

   *y = x + i*size;
}


long _ntl_zblock_destroy(_ntl_verylong x)
{
   long m, *p;

   p = x - 2;
   m = *p;
   free(p);
   return m;
}


long _ntl_zblock_storage(long d)
{
   long size = d+3;
   return size * (sizeof (long)) + sizeof(_ntl_verylong);
}




/* The following routines provide special support for ZZ_pX
 * arithmetic. */



/* this is a generic single-precision mul mod that will work
 * on any platform */


#define SP_MUL_MOD(r, a, b, n)  \
{  \
   long l__a = (a);  \
   long l__b = (b);  \
   long l__n = (n);  \
   long l__q;  \
   unsigned long l__res;  \
  \
   l__q  = (long) ((((double) l__a) * ((double) l__b)) / ((double) l__n));  \
   l__res = ((unsigned long) l__a)*((unsigned long) l__b) - \
            ((unsigned long) l__q)*((unsigned long) l__n);  \
   if (l__res >> (NTL_BITS_PER_LONG-1))  \
      l__res += l__n;  \
   else if (((long) l__res) >= l__n)  \
      l__res -= l__n;  \
  \
   r = (long) l__res;  \
}




static
void sp_ext_eucl(long *dd, long *ss, long *tt, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) zhalt("integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) zhalt("integer overflow");
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

   *dd = u;
   *ss = u1;
   *tt = v1;
}


static
long sp_inv_mod(long a, long n)
{
   long d, s, t;

   sp_ext_eucl(&d, &s, &t, a, n);
   if (d != 1) zhalt("inverse undefined");
   if (s < 0)
      return s + n;
   else
      return s;
}


/* Data structures and algorithms for fast Chinese Remaindering */


struct crt_body_lip {
   _ntl_verylong *v;
   long sbuf;
   long n;
};

#ifdef NTL_GMP_HACK

struct crt_body_gmp {
   MP_INT *v;
   long sbuf;
   long n;
   mpz_t buf;
};

struct crt_body_gmp1 {
   long n;
   long levels;
   long *primes;
   long *inv_vec;
   long *val_vec;
   long *index_vec;
   MP_INT *prod_vec;
   MP_INT *rem_vec;
   MP_INT *coeff_vec;
   MP_INT temps[2];
   mpz_t modulus;
};

#endif


struct crt_body {
   long strategy;

   union {
      struct crt_body_lip L;

#ifdef NTL_GMP_HACK
      struct crt_body_gmp G;
      struct crt_body_gmp1 G1;
#endif
   } U;
};


long _ntl_crt_struct_special(void *crt_struct)
{
   struct crt_body *c = (struct crt_body *) crt_struct;
   return (c->strategy == 2);
}


void _ntl_crt_struct_init(void **crt_struct, long n, _ntl_verylong p,
                          const long *primes)
{
   struct crt_body *c;

   c = (struct crt_body *) NTL_MALLOC(1, sizeof(struct crt_body), 0);
   if (!c) zhalt("out of memory");


#ifdef NTL_GMP_HACK
   if (_ntl_gmp_hack && n >= 600) {
      struct crt_body_gmp1 *C = &c->U.G1;
      long *q;
      long i, j;
      long levels, vec_len;
      long *val_vec, *inv_vec;
      long *index_vec;
      MP_INT *prod_vec, *rem_vec, *coeff_vec;
      MP_INT *temps;

      mpz_init(C->modulus);
      lip_to_mpz(p, C->modulus);

      temps = &C->temps[0];

      mpz_init(&temps[0]);
      mpz_init(&temps[1]);
   
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) zhalt("out of memory");

      val_vec = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!val_vec) zhalt("out of memory");

      inv_vec = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!inv_vec) zhalt("out of memory");

      for (i = 0; i < n; i++)
         q[i] = primes[i];

      levels = 0;
      while ((n >> levels) >= 16) levels++;

      vec_len = (1L << levels) - 1;

      index_vec = (long *) NTL_MALLOC((vec_len+1), sizeof(long), 0);
      if (!index_vec) zhalt("out of memory");

      prod_vec = (MP_INT *) NTL_MALLOC(vec_len, sizeof(MP_INT), 0);
      if (!prod_vec) zhalt("out of memory");

      rem_vec = (MP_INT *) NTL_MALLOC(vec_len, sizeof(MP_INT), 0);
      if (!rem_vec) zhalt("out of memory");

      coeff_vec = (MP_INT *) NTL_MALLOC(n, sizeof(MP_INT), 0);
      if (!coeff_vec) zhalt("out of memory");

      for (i = 0; i < vec_len; i++)
         mpz_init(&prod_vec[i]);

      for (i = 0; i < vec_len; i++)
         mpz_init(&rem_vec[i]);

      for (i = 0; i < n; i++)
         mpz_init(&coeff_vec[i]);

      index_vec[0] = 0;
      index_vec[1] = n;

      for (i = 0; i <= levels-2; i++) {
         long start = (1L << i) - 1;
         long finish = (1L << (i+1)) - 2;
         for (j = finish; j >= start; j--) {
            index_vec[2*j+2] = index_vec[j] + (index_vec[j+1] - index_vec[j])/2;
            index_vec[2*j+1] = index_vec[j];
         }
         index_vec[2*finish+3] = n;
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         /* multiply primes index_vec[i]..index_vec[i+1]-1 into 
          * prod_vec[i]
          */

         mpz_set_ui(&prod_vec[i], 1);
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            mpz_mul_ui(&prod_vec[i], &prod_vec[i], q[j]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            mpz_fdiv_q_ui(&coeff_vec[j], &prod_vec[i], q[j]);
      }

      for (i = (1L << (levels-1)) - 2; i >= 0; i--)
         mpz_mul(&prod_vec[i], &prod_vec[2*i+1], &prod_vec[2*i+2]);


      /* the following is asymptotically the bottleneck...but it
       * it probably doesn't matter. */

      for (i = 0; i < n; i++) {
         long tt;
         mpz_fdiv_q_ui(&temps[0], &prod_vec[0], q[i]);
         tt = mpn_mod_1(temps[0]._mp_d, temps[0]._mp_size, q[i]);
         inv_vec[i] = sp_inv_mod(tt, q[i]);
      }

      c->strategy = 2;
      C->n = n;
      C->primes = q;
      C->val_vec = val_vec;
      C->inv_vec = inv_vec;
      C->levels = levels;
      C->index_vec = index_vec;
      C->prod_vec = prod_vec;
      C->rem_vec = rem_vec;
      C->coeff_vec = coeff_vec;

      *crt_struct = (void *) c;
      return;
   }

   if (_ntl_gmp_hack) {
      struct crt_body_gmp *C = &c->U.G;
      long i;
      c->strategy = 1;

      C->n = n;
      C->v = (MP_INT *) NTL_MALLOC(n, sizeof(MP_INT), 0);
      if (!C->v) zhalt("out of memory");

      for (i = 0; i < n; i++)
         mpz_init(&C->v[i]);

      C->sbuf = L_TO_G(p[0])+2;

      mpz_init(C->buf);
      _mpz_realloc(C->buf, C->sbuf);

      *crt_struct = (void *) c;
      return;
   }
   
#endif

   {
      struct crt_body_lip *C = &c->U.L;
      long i;
      c->strategy = 0;
   
      C->n = n;
   
      C->v = (_ntl_verylong *) NTL_MALLOC(n, sizeof(_ntl_verylong), 0);
      if (!C->v) zhalt("out of memory");
   
      for (i = 0; i < n; i++) 
         C->v[i] = 0;
   
      C->sbuf = p[0]+3;
   
      *crt_struct = (void *) c;
   }
}

void _ntl_crt_struct_insert(void *crt_struct, long i, _ntl_verylong m)
{
   struct crt_body *c = (struct crt_body *) crt_struct;

   switch (c->strategy) {

   case 0: {
      _ntl_zcopy(m, &c->U.L.v[i]);
      break;
   }

#ifdef NTL_GMP_HACK
   case 1: {
      lip_to_mpz(m, &c->U.G.v[i]);
      break;
   }
#endif

   default:
      zhalt("_ntl_crt_struct_insert: inconsistent strategy");
      
   } /* end switch */
}

void _ntl_crt_struct_free(void *crt_struct)
{
   struct crt_body *c = (struct crt_body *) crt_struct;

   switch (c->strategy) {

   case 0: {
      struct crt_body_lip *C = &c->U.L;
      long i, n;

      n = C->n;
   
      for (i = 0; i < n; i++) 
         _ntl_zfree(&C->v[i]);

      free(C->v);
   
      free(c);

      break;
   }

#ifdef NTL_GMP_HACK

   case 1: {
      struct crt_body_gmp *C = &c->U.G;
      long i, n;

      n = C->n;
   
      for (i = 0; i < n; i++) 
         mpz_clear(&C->v[i]);

      mpz_clear(C->buf);

      free(C->v);
   
      free(c);
      break;
   }

   case 2: { 
      struct crt_body_gmp1 *C = &c->U.G1;
      long n = C->n;
      long levels = C->levels;
      long *primes = C->primes;
      long *inv_vec = C->inv_vec;
      long *val_vec = C->val_vec;
      long *index_vec = C->index_vec;
      MP_INT *prod_vec = C->prod_vec;
      MP_INT *rem_vec = C->rem_vec;
      MP_INT *coeff_vec = C->coeff_vec;
      MP_INT *temps = C->temps;
      MP_INT *modulus = C->modulus;
      long vec_len = (1L << levels) - 1;

      long i;

      for (i = 0; i < vec_len; i++)
         mpz_clear(&prod_vec[i]);

      for (i = 0; i < vec_len; i++)
         mpz_clear(&rem_vec[i]);

      for (i = 0; i < n; i++)
         mpz_clear(&coeff_vec[i]);

      mpz_clear(&temps[0]);
      mpz_clear(&temps[1]);

      mpz_clear(modulus);

      free(primes);
      free(inv_vec);
      free(val_vec);
      free(index_vec);
      free(prod_vec);
      free(rem_vec);
      free(coeff_vec);

      free(c);
      break;
   }

#endif

   default:

      zhalt("_ntl_crt_struct_free: inconsistent strategy");

   } /* end case */
}



#ifdef NTL_GMP_HACK


static
void mpz_add_mul_many(MP_INT *res, MP_INT *a, long *b, 
                      long n, long sz)

{
   mp_limb_t *xx, *yy; 
   long i, sx;
   long sy;
   mp_limb_t carry;

   sx = sz + 2;
   if (sx > res->_mp_alloc)
      _mpz_realloc(res, sx);

   xx = res->_mp_d;

   for (i = 0; i < sx; i++)
      xx[i] = 0;

   for (i = 0; i < n; i++) {
      yy = a[i]._mp_d;
      sy = a[i]._mp_size; 

      if (!sy || !b[i]) continue;

      carry = mpn_addmul_1(xx, yy, sy, b[i]);
      yy = xx + sy;
      *yy += carry;

      if (*yy < carry) { /* unsigned comparison! */
         do {
            yy++;
            *yy += 1;
         } while (*yy == 0);
      }
   }

   while (sx > 0 && xx[sx-1] == 0) sx--;
   res->_mp_size = sx;
}

#endif

void _ntl_crt_struct_eval(void *crt_struct, _ntl_verylong *x, const long *b)
{
   struct crt_body *c = (struct crt_body *) crt_struct;

   switch (c->strategy) {

   case 0: {
      struct crt_body_lip *C = &c->U.L;

      _ntl_verylong xx, yy, *a;
      long i, sx, n;
   
      n = C->n;
      sx = C->sbuf;
   
      _ntl_zsetlength(x, sx);
      xx = *x;
   
      a = C->v;
   
   
      for (i = 1; i <= sx; i++)
         xx[i] = 0;
   
      xx++;
   
      for (i = 0; i < n; i++) {
         yy = a[i];
   
         if (!yy || !b[i]) continue;
   
         zaddmul(b[i], xx, yy);
         yy = xx + yy[0];
     
         if ((*yy) >= NTL_RADIX) {
            (*yy) -= NTL_RADIX;
            yy++;
            while ((*yy) == NTL_RADIX-1) {
               *yy = 0;
               yy++;
            }
            (*yy)++;
         }
      }
   
      xx--;
      while (sx > 1 && xx[sx] == 0) sx--;
      xx[0] = sx;

      break;
   }

#ifdef NTL_GMP_HACK

   case 1: {
      struct crt_body_gmp *C = &c->U.G;

      mp_limb_t *xx, *yy; 
      MP_INT *a;
      long i, sx, n;
      long sy;
      mp_limb_t carry;
   
      n = C->n;
      sx = C->sbuf;
   
      xx = C->buf->_mp_d;

      for (i = 0; i < sx; i++)
         xx[i] = 0;
   
      a = C->v;
   
      for (i = 0; i < n; i++) {
         yy = a[i]._mp_d;
         sy = a[i]._mp_size; 
   
         if (!sy || !b[i]) continue;
   
         carry = mpn_addmul_1(xx, yy, sy, b[i]);
         yy = xx + sy;
         *yy += carry;

         if (*yy < carry) { /* unsigned comparison! */
            do {
               yy++;
               *yy += 1;
            } while (*yy == 0);
         }
      }
   
      while (sx > 0 && xx[sx-1] == 0) sx--;
      C->buf->_mp_size = sx;
      mpz_to_lip(x, C->buf);
      break;
   }

   case 2: {
      struct crt_body_gmp1 *C = &c->U.G1;

      long n = C->n;
      long levels = C->levels;
      long *primes = C->primes;
      long *inv_vec = C->inv_vec;
      long *val_vec = C->val_vec;
      long *index_vec = C->index_vec;
      MP_INT *prod_vec = C->prod_vec;
      MP_INT *rem_vec = C->rem_vec;
      MP_INT *coeff_vec = C->coeff_vec;
      MP_INT *temps = C->temps;
      long vec_len = (1L << levels) - 1;

      long i, j;

      for (i = 0; i < n; i++) {
         SP_MUL_MOD(val_vec[i], b[i], inv_vec[i], primes[i]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         long j1 = index_vec[i];
         long j2 = index_vec[i+1];
         mpz_add_mul_many(&rem_vec[i], &coeff_vec[j1], &val_vec[j1], j2-j1, 
                          prod_vec[i]._mp_size);
      }

      for (i = (1L << (levels-1)) - 2; i >= 0; i--) {
         mpz_mul(&temps[0], &prod_vec[2*i+1], &rem_vec[2*i+2]);
         mpz_mul(&temps[1], &rem_vec[2*i+1], &prod_vec[2*i+2]);
         mpz_add(&rem_vec[i], &temps[0], &temps[1]);
      }

      /* temps[0] = rem_vec[0] mod prod_vec[0] (least absolute residue) */
      mpz_fdiv_r(&temps[0], &rem_vec[0], &prod_vec[0]);
      mpz_sub(&temps[1], &temps[0], &prod_vec[0]);
      if (mpz_cmpabs(&temps[0], &temps[1]) > 0)
         mpz_set(&temps[0], &temps[1]);

      mpz_fdiv_r(&temps[1], &temps[0], C->modulus);
      mpz_to_lip(x, &temps[1]);

      break;
   }
#endif

   default:

      zhalt("_crt_struct_eval: inconsistent strategy");

   } /* end case */
}






/* Data structures and algorithms for multi-modulus remaindering */


struct rem_body_lip {
   long n;
   long *primes;
};

#ifdef NTL_GMP_HACK

struct rem_body_gmp {
   long n;
   long levels;
   long *primes;
   long *index_vec;
   MP_INT *prod_vec;
   MP_INT *rem_vec;
};

#endif



#ifdef NTL_SINGLE_MUL 

struct rem_body_sm {
   long n;
   long *primes;
   double **tbl;
};

#endif


#if (defined(NTL_TBL_REM) && !defined(NTL_SINGLE_MUL))

struct rem_body_tbl {
   long n;
   long *primes;
   long **tbl;
};

#endif

struct rem_body {
   long strategy;

   union {
      struct rem_body_lip L;
#ifdef NTL_GMP_HACK
      struct rem_body_gmp G;
#endif
#ifdef NTL_SINGLE_MUL 
      struct rem_body_sm S;
#endif
#if (defined(NTL_TBL_REM) && !defined(NTL_SINGLE_MUL))
      struct rem_body_tbl T;
#endif
   } U;
};



void _ntl_rem_struct_init(void **rem_struct, long n, _ntl_verylong modulus, 
                          const long *p)
{
   struct rem_body *r;

   r = (struct rem_body *) NTL_MALLOC(1, sizeof(struct rem_body), 0);
   if (!r) zhalt("out of memory");
   
#ifdef NTL_GMP_HACK
   if (_ntl_gmp_hack && n >= 32) {
      struct rem_body_gmp *R = &r->U.G;

      long *q;
      long i, j;
      long levels, vec_len;
      long *index_vec;
      MP_INT *prod_vec, *rem_vec;
   
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) zhalt("out of memory");
   
      for (i = 0; i < n; i++)
         q[i] = p[i];

      levels = 0;
      while ((n >> levels) >= 4) levels++;

      vec_len = (1L << levels) - 1;

      index_vec = (long *) NTL_MALLOC((vec_len+1), sizeof(long), 0);
      if (!index_vec) zhalt("out of memory");

      prod_vec = (MP_INT *) NTL_MALLOC(vec_len, sizeof(MP_INT), 0);
      if (!prod_vec) zhalt("out of memory");

      rem_vec = (MP_INT *) NTL_MALLOC(vec_len, sizeof(MP_INT), 0);
      if (!rem_vec) zhalt("out of memory");

      for (i = 0; i < vec_len; i++)
         mpz_init(&prod_vec[i]);

      for (i = 0; i < vec_len; i++)
         mpz_init(&rem_vec[i]);

      index_vec[0] = 0;
      index_vec[1] = n;

      for (i = 0; i <= levels-2; i++) {
         long start = (1L << i) - 1;
         long finish = (1L << (i+1)) - 2;
         for (j = finish; j >= start; j--) {
            index_vec[2*j+2] = index_vec[j] + (index_vec[j+1] - index_vec[j])/2;
            index_vec[2*j+1] = index_vec[j];
         }
         index_vec[2*finish+3] = n;
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         /* multiply primes index_vec[i]..index_vec[i+1]-1 into 
          * prod_vec[i]
          */

         mpz_set_ui(&prod_vec[i], 1);
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            mpz_mul_ui(&prod_vec[i], &prod_vec[i], q[j]);
      }

      for (i = (1L << (levels-1)) - 2; i >= 3; i--)
         mpz_mul(&prod_vec[i], &prod_vec[2*i+1], &prod_vec[2*i+2]);

      r->strategy = 1;
      R->n = n;
      R->primes = q;
      R->levels = levels;
      R->index_vec = index_vec;
      R->prod_vec = prod_vec;
      R->rem_vec = rem_vec;

      *rem_struct = (void *) r;

      return;
   }
#endif

#ifdef NTL_SINGLE_MUL
   {
      struct rem_body_sm *R = &r->U.S;

      long *qq;
      long i;
      double **tbl;
      long t, t1, j, q;
      long sz = modulus[0];
   
      r->strategy = 2;
      R->n = n;
      qq = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!qq) zhalt("out of memory");
      R->primes = qq;
   
      for (i = 0; i < n; i++)
         qq[i] = p[i];

      tbl = (double **) NTL_MALLOC(n, sizeof(double *), 0);
      if (!tbl) zhalt("out of space");

      for (i = 0; i < n; i++) {
         tbl[i] = (double *) NTL_MALLOC(sz, sizeof(double), 0);
         if (!tbl[i]) zhalt("out of space");
      }
      R->tbl = tbl;
   
   
      for (i = 0; i < n; i++) {
         q = qq[i];
         t = (((long)1) << NTL_NBITS) % q;
         t1 = 1;
         tbl[i][0] = (double) 1;
         for (j = 1; j < sz; j++) {
            SP_MUL_MOD(t1, t1, t, q);
            tbl[i][j] = (double) t1;
         }
      }
      
      *rem_struct = (void *) r;
      return;
   }

#endif

#if (defined(NTL_TBL_REM) && !defined(NTL_SINGLE_MUL))

   {
      struct rem_body_tbl *R = &r->U.T;

      long *qq;
      long i;
      long **tbl;
      long t, t1, j, q;
      long sz = modulus[0];
   
      r->strategy = 3;
      R->n = n;
      qq = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!qq) zhalt("out of memory");
      R->primes = qq;
   
      for (i = 0; i < n; i++)
         qq[i] = p[i];

      tbl = (long **) NTL_MALLOC(n, sizeof(long *), 0);
      if (!tbl) zhalt("out of space");

      for (i = 0; i < n; i++) {
         tbl[i] = (long *) NTL_MALLOC(sz, sizeof(long), 0);
         if (!tbl[i]) zhalt("out of space");
      }
      R->tbl = tbl;
   
   
      for (i = 0; i < n; i++) {
         q = qq[i];
         t = (((long)1) << NTL_NBITS) % q;
         t1 = 1;
         tbl[i][0] = 1;
         for (j = 1; j < sz; j++) {
            SP_MUL_MOD(t1, t1, t, q);
            tbl[i][j] = t1;
         }
      }
      
      *rem_struct = (void *) r;
      return;
   }

#endif


   {
      struct rem_body_lip *R = &r->U.L;

      long *q;
      long i;
   
      r->strategy = 0;
      R->n = n;
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) zhalt("out of memory");
      R->primes = q;
   
      for (i = 0; i < n; i++)
         q[i] = p[i];
   
      *rem_struct = (void *) r;
   }
}

void _ntl_rem_struct_free(void *rem_struct)
{
   struct rem_body *r = (struct rem_body *) rem_struct;

   switch (r->strategy) {

   case 0: {
      free(r->U.L.primes);
      free(r);
      break;
   }

#ifdef NTL_GMP_HACK

   case 1: {
      struct rem_body_gmp *R = &r->U.G;

      long levels = R->levels;
      long vec_len = (1L << levels) - 1;
      long i;

      for (i = 0; i < vec_len; i++) 
         mpz_clear(&R->prod_vec[i]);

      for (i = 0; i < vec_len; i++) 
         mpz_clear(&R->rem_vec[i]);

      free(R->primes);
      free(R->index_vec);
      free(R->prod_vec);
      free(R->rem_vec);
      free(r);
      break;
   }

#endif


#ifdef NTL_SINGLE_MUL

   case 2: {
      struct rem_body_sm *R = &r->U.S;

      long n = R->n;
      double **tbl = R->tbl;
      long i;

      for (i = 0; i < n; i++)
         free(tbl[i]);

      free(tbl);

      free(R->primes);
      free(r);
      break;
   }

#endif


#if (defined(NTL_TBL_REM) && !defined(NTL_SINGLE_MUL))
   
   case 3: {
      struct rem_body_tbl *R = &r->U.T;

      long n = R->n;
      long **tbl = R->tbl;
      long i;

      for (i = 0; i < n; i++)
         free(tbl[i]);

      free(tbl);

      free(R->primes);
      free(r);
      break;
   }


#endif

   default:
      zhalt("_ntl_rem_struct_free: inconsistent strategy");

   } /* end switch */

}


void _ntl_rem_struct_eval(void *rem_struct, long *x, _ntl_verylong a)
{
   struct rem_body *r = (struct rem_body *) rem_struct;

   switch (r->strategy) {

   case 0: {
      struct rem_body_lip *R = &r->U.L;

      _ntl_zmultirem(a, R->n, R->primes, x);
      return;
   }

#ifdef NTL_GMP_HACK

   case 1: {
      struct rem_body_gmp *R = &r->U.G;

      long n = R->n;
      long levels = R->levels;
      long *q = R->primes;
      long *index_vec = R->index_vec;
      MP_INT *prod_vec = R->prod_vec;
      MP_INT *rem_vec = R->rem_vec;
      long vec_len = (1L << levels) - 1;

      long i, j;

      lip_to_mpz(a, &rem_vec[1]);
      mpz_set(&rem_vec[2], &rem_vec[1]);

      for (i = 1; i < (1L << (levels-1)) - 1; i++) {
         mpz_mod(&rem_vec[2*i+1], &rem_vec[i], &prod_vec[2*i+1]);
         mpz_mod(&rem_vec[2*i+2], &rem_vec[i], &prod_vec[2*i+2]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         long lo = index_vec[i];
         long hi = index_vec[i+1];
         mp_limb_t *s1p = rem_vec[i]._mp_d;
         long s1size = rem_vec[i]._mp_size;
         if (s1size == 0) {
            for (j = lo; j <hi; j++)
               x[j] = 0;
         }
         else {
            for (j = lo; j < hi; j++)
               x[j] = mpn_mod_1(s1p, s1size, q[j]);
         }
      }

      break;
   }

#endif


#ifdef NTL_SINGLE_MUL

   case 2: {
      struct rem_body_sm *R = &r->U.S;

      _ntl_zmultirem2(a, R->n, R->primes, R->tbl, x);
      break;
   }

#endif

#if (defined(NTL_TBL_REM) && !defined(NTL_SINGLE_MUL))

   case 3: {
      struct rem_body_tbl *R = &r->U.T;

      _ntl_zmultirem3(a, R->n, R->primes, R->tbl, x);
      break;
   }

#endif

   default:
      zhalt("_ntl_rem_struct_eval: inconsistent strategy");


   } /* end switch */
}



void
_ntl_zaorsmul_1(_ntl_verylong x, long y, long sub,  _ntl_verylong *ww)
{
   static _ntl_verylong tmp = 0;

   if (y == 0) return;

   if (y == 1) {
      if (sub)
         _ntl_zsub(*ww, x, ww);
      else
         _ntl_zadd(*ww, x, ww);
      return;
   }

   if (y == -1) {
      if (!sub)
         _ntl_zsub(*ww, x, ww);
      else
         _ntl_zadd(*ww, x, ww);
      return;
   }

   _ntl_zsmul(x, y, &tmp);
   if (sub)
      _ntl_zsub(*ww, tmp, ww);
   else
      _ntl_zadd(*ww, tmp, ww);
}


void
_ntl_zsaddmul(_ntl_verylong x, long y,  _ntl_verylong *ww)
{
  _ntl_zaorsmul_1(x, y, 0, ww);
}

void
_ntl_zssubmul(_ntl_verylong x, long y,  _ntl_verylong *ww)
{
  _ntl_zaorsmul_1(x, y, 1, ww);
}






void
_ntl_zaorsmul(_ntl_verylong x, _ntl_verylong y, long sub,  _ntl_verylong *ww)
{
   static _ntl_verylong tmp = 0;

   _ntl_zmul(x, y, &tmp);
   if (sub)
      _ntl_zsub(*ww, tmp, ww);
   else
      _ntl_zadd(*ww, tmp, ww);
}


void
_ntl_zaddmul(_ntl_verylong x, _ntl_verylong y,  _ntl_verylong *ww)
{
  _ntl_zaorsmul(x, y, 0, ww);
}

void
_ntl_zsubmul(_ntl_verylong x, _ntl_verylong y,  _ntl_verylong *ww)
{
  _ntl_zaorsmul(x, y, 1, ww);
}


