
/*
 * This is a "wrapper" layer that builds on top of the "mpn" layer of gmp.
 * This layer provides much of the same functionality of the "mpz"
 * layer of gmp, but the interface it provides is much more like
 * the interface provided by lip.
 *
 * This layer was written under the following assumptions about gmp:
 *  1) mp_limb_t is an unsigned integral type
 *  2) sizeof(mp_limb_t) == sizeof(long) or sizeof(mp_limb_t) == 2*sizeof(long)
 *  3) the number of bits of an mp_limb_t is equal to that of a long,
 *     or twice that of a long
 *  4) the number of bits of a gmp radix is equal to the number of bits
 *     of an mp_limb_t
 *
 * Except for assumption (1), these assumptions are verified in the
 * installation script, and they should be universally satisfied in practice,
 * except when gmp is built using the proposed, new "nail" fetaure
 * (in which some bits of an mp_limb_t are unused).
 * The code here will not work properly with the "nail" feature;
 * however, I have (attempted to) identify all such problem spots,
 * and any other places where assumptions (2-4) are made,
 * with a comment labeled "DIRT".
 */



#include <NTL/lip.h>

#include <NTL/ctools.h>


#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include <gmp.h>

typedef mp_limb_t *_ntl_limb_t_ptr;

int _ntl_gmp_hack = 0;

#if (__GNU_MP_VERSION < 3)

#error "You have to use GNP version >= 3.1"

#endif

#if ((__GNU_MP_VERSION == 3) && (__GNU_MP_VERSION_MINOR < 1))

#error "You have to use GNP version >= 3.1"

#endif



/* v 3.1 is supposed  mpn_tdiv_qr defined, but it doesn't.
   Here's a workaround */

#if ((__GNU_MP_VERSION == 3) && (__GNU_MP_VERSION_MINOR == 1) && (__GNU_MP_VERSION_PATCHLEVEL == 0))

#define mpn_tdiv_qr __MPN(tdiv_qr)


#ifdef __cplusplus
extern "C" 
#endif
void mpn_tdiv_qr(mp_ptr, mp_ptr, mp_size_t, mp_srcptr, mp_size_t, 
                 mp_srcptr, mp_size_t);

#endif


#if (defined(NTL_CXX_ONLY) && !defined(__cplusplus))
#error "CXX_ONLY flag set...must use C++ compiler"
#endif


union gbigint_header {

   long info[2];
   mp_limb_t alignment;

};

/* A bigint is represented as two long's, ALLOC and SIZE, followed by a 
 * vector DATA of mp_limb_t's.  
 * 
 * ALLOC is of the form
 *    (alloc << 2) | continue_flag | frozen_flag
 * where 
 *    - alloc is the number of allocated mp_limb_t's,
 *    - continue flag is either 2 or 0,
 *    - frozen_flag is either 1 or 0.
 * If frozen_flag is set, then the space for this bigint is *not*
 * managed by the _ntl_gsetlength and _ntl_gfree routines,
 * but are instead managed by the vec_ZZ_p and ZZVec routines.
 * The continue_flag is only set when the frozen_flag is set.
 * 
 * SIZE is the number of mp_limb_t's actually
 * used by the bigint, with the sign of SIZE having
 * the sign of the bigint.
 * Note that the zero bigint is represented as SIZE=0.
 * 
 * Bigint's are accessed through a handle, which is pointer to void.
 * A null handle logically represents the bigint zero.
 * This is done so that the interface presented to higher level
 * routines is essentially the same as that of NTL's traditional
 * long integer package.
 * 
 * The components ALLOC, SIZE, and DATA are all accessed through
 * macros using pointer casts.  While all of may seem a bit dirty, 
 * it should be quite portable: objects are never referenced
 * through pointers of different types, and no alignmement
 * problems should arise.
 * 
 * Actually, mp_limb_t is usually the type unsigned long.
 * However, on some 64-bit platforms, the type long is only 32 bits,
 * and gmp makes mp_limb_t unsigned long long in this case.
 * This is fairly rare, as the industry standard for Unix is to
 * have 64-bit longs on 64-bit machines.
 */ 

#if 1

#define ALLOC(p) (((long *) (p))[0])
#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

#define STORAGE(len) ((long)(2*sizeof(long) + (len)*sizeof(mp_limb_t)))

/* DIRT: STORAGE computes the number of bytes to allocate for a bigint
 * of maximal SIZE len.  This should be computed so that one
 * can store several such bigints in a contiguous array
 * of memory without breaking any alignment requirements.
 * Currently, it is assumed (and explicitly checked in the NTL installation
 * script) that sizeof(mp_limb_t) is either sizeof(long) or
 * 2*sizeof(long), and therfore, nothing special needs to
 * be done to enfoce alignment requirements.  If this assumption
 * should change, then the storage layout for bigints must be
 * re-designed.   
 */

#define MustAlloc(c, len)  (!(c) || (ALLOC(c) >> 2) < (len))



#define GET_SIZE_NEG(sz, neg, p)  \
do  \
{   \
   long _s;   \
   _s = SIZE(p);   \
   if (_s < 0) {  \
      sz = -_s;  \
      neg = 1;  \
   }  \
   else {  \
      sz = _s;  \
      neg = 0;  \
   }  \
}  \
while (0)

#define STRIP(sz, p)  \
do  \
{  \
   long _i;  \
   _i = sz - 1;  \
   while (_i >= 0 && p[_i] == 0) _i--;  \
   sz = _i + 1;  \
}  \
while (0) 

#define ZEROP(p) (!p || !SIZE(p))

#define ONEP(p) (p && SIZE(p) == 1 && DATA(p)[0] == 1)

#define SWAP_BIGINT(a, b)  \
do  \
{  \
   _ntl_gbigint _t;  \
   _t = a;  \
   a = b;  \
   b = _t;  \
}  \
while (0) 

#define SWAP_LONG(a, b)  \
do  \
{  \
   long _t;  \
   _t = a;  \
   a = b;  \
   b = _t;  \
}  \
while (0) 

#define SWAP_LIMB_PTR(a, b)  \
do  \
{  \
   _ntl_limb_t_ptr _t;  \
   _t = a;  \
   a = b;  \
   b = _t;  \
}  \
while (0) 

#define COUNT_BITS(cnt, a)  \
do  \
{  \
   long _i = 0;  \
   mp_limb_t _a = (a);  \
  \
   while (_a>=256)  \
      _i += 8, _a >>= 8;  \
   if (_a >=16)  \
      _i += 4, _a >>= 4;  \
   if (_a >= 4)  \
      _i += 2, _a >>= 2;  \
   if (_a >= 2)  \
      _i += 2;  \
   else if (_a >= 1)  \
      _i++;  \
  \
   cnt = _i;  \
}  \
while (0) 

#else

/* These are C++ inline functions that are equivalent to the above
 * macros.  They are mainly intended as a debugging aid.
 */


inline long& ALLOC(_ntl_gbigint p) 
   { return (((long *) p)[0]); }

inline long& SIZE(_ntl_gbigint p) 
   { return (((long *) p)[1]); }

inline mp_limb_t * DATA(_ntl_gbigint p) 
   { return ((mp_limb_t *) (((long *) (p)) + 2)); }

inline long STORAGE(long len)
   { return ((long)(2*sizeof(long) + (len)*sizeof(mp_limb_t))); }

inline long MustAlloc(_ntl_gbigint c, long len)  
   { return (!(c) || (ALLOC(c) >> 2) < (len)); }


inline void GET_SIZE_NEG(long& sz, long& neg, _ntl_gbigint p)
{ 
   long s; 
   s = SIZE(p); 
   if (s < 0) {
      sz = -s;
      neg = 1;
   }
   else {
      sz = s;
      neg = 0;
   }
}

inline void STRIP(long& sz, mp_limb_t *p)
{
   long i;
   i = sz - 1;
   while (i >= 0 && p[i] == 0) i--;
   sz = i + 1;
}

inline long ZEROP(_ntl_gbigint p)
{
   return !p || !SIZE(p);
}

inline long ONEP(_ntl_gbigint p)
{
   return p && SIZE(p) == 1 && DATA(p)[0] == 1;
}

inline void SWAP_BIGINT(_ntl_gbigint& a, _ntl_gbigint& b)
{
   _ntl_gbigint t;
   t = a;
   a = b;
   b = t;
}

inline void SWAP_LONG(long& a, long& b)
{
   long t;
   t = a;
   a = b;
   b = t;
}

inline void SWAP_LIMB_PTR(_ntl_limb_t_ptr& a, _ntl_limb_t_ptr& b)
{
   _ntl_limb_t_ptr t;
   t = a;
   a = b;
   b = t;
}


inline void COUNT_BITS(long& cnt, mp_limb_t a)
{
   long i = 0;

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

   cnt = i;
}

#endif

#define STORAGE_OVF(len) NTL_OVERFLOW(len, sizeof(mp_limb_t), 2*sizeof(long))


/* ForceNormal ensures a normalized bigint */

static 
void ForceNormal(_ntl_gbigint x)
{
   long sx, xneg;
   mp_limb_t *xdata;

   if (!x) return;
   GET_SIZE_NEG(sx, xneg, x);
   xdata = DATA(x);
   STRIP(sx, xdata);
   if (xneg) sx = -sx;
   SIZE(x) = sx;
}


static 
void ghalt(char *c)
{
   fprintf(stderr,"fatal error:\n   %s\nexit...\n",c);
   fflush(stderr);
   _ntl_abort();
}


long _ntl_gmaxalloc(_ntl_gbigint x)
{
   if (!x)
      return 0;
   else
      return (ALLOC(x) >> 2) - 1;
}

#define MIN_SETL	(4)
   /* _ntl_gsetlength allocates a multiple of MIN_SETL digits */



void _ntl_gsetlength(_ntl_gbigint *v, long len)
{
   _ntl_gbigint x = *v;

   if (len < 0)
      ghalt("negative size allocation in _ntl_zgetlength");

   if (NTL_OVERFLOW(len, NTL_ZZ_NBITS, 0))
      ghalt("size too big in _ntl_gsetlength");

#ifdef NTL_SMALL_MP_SIZE_T
   /* this makes sure that numbers don't get too big for GMP */
   if (len >= (1L << (NTL_BITS_PER_INT-4)))
      ghalt("size too big for GMP");
#endif


   if (x) {
      long oldlen = ALLOC(x);
      long fixed = oldlen & 1;
      oldlen = oldlen >> 2;

      if (fixed) {
         if (len > oldlen) 
            ghalt("internal error: can't grow this _ntl_gbigint");
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
      if (NTL_OVERFLOW(len, NTL_ZZ_NBITS, 0))
         ghalt("size too big in _ntl_gsetlength");

      if (STORAGE_OVF(len))
         ghalt("reallocation failed in _ntl_gsetlength");

      ALLOC(x) = len << 2;
      if (!(x = (_ntl_gbigint)NTL_REALLOC((void *) x, 1, STORAGE(len), 0))) {
         ghalt("reallocation failed in _ntl_gsetlength");
      }
   }
   else {
      len++;  /* as above, always allocate one more than explicitly reqested */
      len = ((len+(MIN_SETL-1))/MIN_SETL)*MIN_SETL; 

      /* test len again */
      if (NTL_OVERFLOW(len, NTL_ZZ_NBITS, 0))
         ghalt("size too big in _ntl_gsetlength");

      if (STORAGE_OVF(len))
         ghalt("reallocation failed in _ntl_gsetlength");

      if (!(x = (_ntl_gbigint)NTL_MALLOC(1, STORAGE(len), 0))) {
         ghalt("allocation failed in _ntl_gsetlength");
      }
      ALLOC(x) = len << 2;
      SIZE(x) = 0;
   }

   *v = x;
}

void _ntl_gfree(_ntl_gbigint *xx)
{
   _ntl_gbigint x = *xx;

   if (!x)
      return;

   if (ALLOC(x) & 1)
      ghalt("Internal error: can't free this _ntl_gbigint");

   free((void*) x);
   *xx = 0;
   return;
}

void
_ntl_gswap(_ntl_gbigint *a, _ntl_gbigint *b)
{
   _ntl_gbigint c;

   if ((*a && (ALLOC(*a) & 1)) || (*b && (ALLOC(*b) & 1))) {
      static _ntl_gbigint t = 0; 
      _ntl_gcopy(*a, &t);
      _ntl_gcopy(*b, a);
      _ntl_gcopy(t, b);
      return;
   }

   c = *a;
   *a = *b;
   *b = c;
}


void _ntl_gcopy(_ntl_gbigint a, _ntl_gbigint *bb)
{
   _ntl_gbigint b;
   long sa, abs_sa, i;
   mp_limb_t *adata, *bdata;

   b = *bb;

   if (!a || (sa = SIZE(a)) == 0) {
      if (b) SIZE(b) = 0;
   }
   else {
      if (a != b) {
         if (sa >= 0)
            abs_sa = sa;
         else
            abs_sa = -sa;

         if (MustAlloc(b, abs_sa)) {
            _ntl_gsetlength(&b, abs_sa);
            *bb = b;
         }

         adata = DATA(a);
         bdata = DATA(b);

         for (i = 0; i < abs_sa; i++)
            bdata[i] = adata[i];

         SIZE(b) = sa;
      }
   }
}


void _ntl_gzero(_ntl_gbigint *aa) 
{
   _ntl_gbigint a = *aa;

   if (a) SIZE(a) = 0;
}

void _ntl_gone(_ntl_gbigint *aa)
{
   _ntl_gbigint a = *aa;
   if (!a) {
      _ntl_gsetlength(&a, 1);
      *aa = a;
   }

   SIZE(a) = 1;
   DATA(a)[0] = 1;
}

long _ntl_giszero(_ntl_gbigint a)
{
   return ZEROP(a);
}

long _ntl_godd(_ntl_gbigint a)
{
   if (ZEROP(a)) 
      return 0;
   else
      return DATA(a)[0]&1;
}

long _ntl_gbit(_ntl_gbigint a, long p)
{
   long bl;
   long sa;
   mp_limb_t wh;

   if (p < 0 || !a) return 0;

   bl = p/NTL_ZZ_NBITS;
   wh = ((mp_limb_t) 1) << (p - NTL_ZZ_NBITS*bl);

   sa = SIZE(a);
   if (sa < 0) sa = -sa;

   if (sa <= bl) return 0;
   if (DATA(a)[bl] & wh) return 1;
   return 0;
}

void _ntl_glowbits(_ntl_gbigint a, long b, _ntl_gbigint *cc)
{
   _ntl_gbigint c;

   long bl;
   long wh;
   long sa;
   long i;
   mp_limb_t *adata, *cdata;

   if (ZEROP(a) || (b<=0)) {
      _ntl_gzero(cc);
      return;
   }

   bl = b/NTL_ZZ_NBITS;
   wh = b - NTL_ZZ_NBITS*bl;
   if (wh != 0) 
      bl++;
   else
      wh = NTL_ZZ_NBITS;

   sa = SIZE(a);
   if (sa < 0) sa = -sa;

   if (sa < bl) {
      _ntl_gcopy(a,cc);
      _ntl_gabs(cc);
      return;
   }

   c = *cc;

   /* a won't move if c aliases a */
   _ntl_gsetlength(&c, bl);
   *cc = c;

   adata = DATA(a);
   cdata = DATA(c);

   for (i = 0; i < bl-1; i++)
      cdata[i] = adata[i];

   if (wh == NTL_ZZ_NBITS)
      cdata[bl-1] = adata[bl-1];
   else
      cdata[bl-1] = adata[bl-1] & ((((mp_limb_t) 1) << wh) - ((mp_limb_t) 1));

   STRIP(bl, cdata);
   SIZE(c) = bl; 
}

long _ntl_gslowbits(_ntl_gbigint a, long p)
{
   static _ntl_gbigint x = 0;

   if (p > NTL_BITS_PER_LONG)
      p = NTL_BITS_PER_LONG;

   _ntl_glowbits(a, p, &x);

   return _ntl_gtoint(x);
}

long _ntl_gsetbit(_ntl_gbigint *a, long b)
{
   long bl;
   long sa, aneg;
   long i;
   mp_limb_t wh, *adata, tmp;

   if (b<0) ghalt("_ntl_gsetbit: negative index");

   if (ZEROP(*a)) {
      _ntl_gintoz(1, a);
      _ntl_glshift(*a, b, a);
      return 0;
   }

   bl = (b/NTL_ZZ_NBITS);
   wh = ((mp_limb_t) 1) << (b - NTL_ZZ_NBITS*bl);

   GET_SIZE_NEG(sa, aneg, *a);

   if (sa > bl) {
      adata = DATA(*a);
      tmp = adata[bl] & wh;
      adata[bl] |= wh;
      if (tmp) return 1;
      return 0;
   } 
   else {
      _ntl_gsetlength(a, bl+1);
      adata = DATA(*a);
      for (i = sa; i < bl; i++)
         adata[i] = 0;
      adata[bl] = wh;

      sa = bl+1;
      if (aneg) sa = -sa;
      SIZE(*a) = sa;
      return 0;
   }
}

long _ntl_gswitchbit(_ntl_gbigint *a, long b)
{
   long bl;
   long sa, aneg;
   long i;
   mp_limb_t wh, *adata, tmp;

   if (b<0) ghalt("_ntl_gswitchbit: negative index");


   if (ZEROP(*a)) {
      _ntl_gintoz(1, a);
      _ntl_glshift(*a, b, a);
      return 0;
   }

   bl = (b/NTL_ZZ_NBITS);
   wh = ((mp_limb_t) 1) << (b - NTL_ZZ_NBITS*bl);

   GET_SIZE_NEG(sa, aneg, *a);

   if (sa > bl) {
      adata = DATA(*a);
      tmp = adata[bl] & wh;
      adata[bl] ^= wh;

      if (bl == sa-1) {
         STRIP(sa, adata);
         if (aneg) sa = -sa;
         SIZE(*a) = sa;
      }

      if (tmp) return 1;
      return 0;
   } 
   else {
      _ntl_gsetlength(a, bl+1);
      adata = DATA(*a);
      for (i = sa; i < bl; i++)
         adata[i] = 0;
      adata[bl] = wh;

      sa = bl+1;
      if (aneg) sa = -sa;
      SIZE(*a) = sa;
      return 0;
   }
}

long
_ntl_gweights(
	long aa
	)
{
	unsigned long a;
	long res = 0;
	if (aa < 0) 
		a = -((unsigned long) aa);
	else
		a = aa;
   
	while (a) {
		if (a & 1) res ++;
		a >>= 1;
	}
	return (res);
}

static long
gweights_mp_limb(
	mp_limb_t a
	)
{
	long res = 0;
   
	while (a) {
		if (a & 1) res ++;
		a >>= 1;
	}
	return (res);
}

long
_ntl_gweight(
        _ntl_gbigint a
        )
{
	long i;
	long sa;
	mp_limb_t *adata;
	long res;

	if (!a) return (0);

	sa = SIZE(a);
	if (sa < 0) sa = -sa;
	adata = DATA(a);

	res = 0;
	for (i = 0; i < sa; i++)
		res += gweights_mp_limb(adata[i]);

	return (res);
}

long 
_ntl_g2logs(
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

long _ntl_g2log(_ntl_gbigint a)
{
   long la;
   long t;

   if (!a) return 0;
   la = SIZE(a);
   if (la == 0) return 0;
   if (la < 0) la = -la;
   COUNT_BITS(t, DATA(a)[la-1]);
   return NTL_ZZ_NBITS*(la - 1) + t;
}



long _ntl_gmakeodd(_ntl_gbigint *nn)
{
   _ntl_gbigint n = *nn;
   long shift;
   mp_limb_t *ndata;
   mp_limb_t i;

   if (ZEROP(n))
      return (0);

   shift = 0;
   ndata = DATA(n);

   while (ndata[shift] == 0)
      shift++;

   i = ndata[shift];

   shift = NTL_ZZ_NBITS * shift;

   while ((i & 1) == 0) {
      shift++;
      i >>= 1;
   }
   _ntl_grshift(n, shift, &n);
   return shift;
}


long _ntl_gnumtwos(_ntl_gbigint n)
{
   long shift;
   mp_limb_t *ndata;
   mp_limb_t i;

   if (ZEROP(n))
      return (0);

   shift = 0;
   ndata = DATA(n);

   while (ndata[shift] == 0)
      shift++;

   i = ndata[shift];

   shift = NTL_ZZ_NBITS * shift;

   while ((i & 1) == 0) {
      shift++;
      i >>= 1;
   }

   return shift;
}


void _ntl_gand(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   _ntl_gbigint c;
   long sa;
   long sb;
   long sm;
   long i;
   long a_alias, b_alias;
   mp_limb_t *adata, *bdata, *cdata;

   if (ZEROP(a) || ZEROP(b)) {
      _ntl_gzero(cc);
      return;
   }

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   sa = SIZE(a);
   if (sa < 0) sa = -sa;

   sb = SIZE(b);
   if (sb < 0) sb = -sb;

   sm = (sa > sb ? sb : sa);

   _ntl_gsetlength(&c, sm);
   if (a_alias) a = c;
   if (b_alias) b = c;
   *cc = c;

   adata = DATA(a);
   bdata = DATA(b);
   cdata = DATA(c);

   for (i = 0; i < sm; i++)
      cdata[i] = adata[i] & bdata[i];

   STRIP(sm, cdata);
   SIZE(c) = sm;
}


void _ntl_gxor(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   _ntl_gbigint c;
   long sa;
   long sb;
   long sm;
   long la;
   long i;
   long a_alias, b_alias;
   mp_limb_t *adata, *bdata, *cdata;

   if (ZEROP(a)) {
      _ntl_gcopy(b,cc);
      _ntl_gabs(cc);
      return;
   }

   if (ZEROP(b)) {
      _ntl_gcopy(a,cc);
      _ntl_gabs(cc);
      return;
   }

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   sa = SIZE(a);
   if (sa < 0) sa = -sa;

   sb = SIZE(b);
   if (sb < 0) sb = -sb;

   if (sa > sb) {
      la = sa;
      sm = sb;
   } 
   else {
      la = sb;
      sm = sa;
   }

   _ntl_gsetlength(&c, la);
   if (a_alias) a = c;
   if (b_alias) b = c;
   *cc = c;

   adata = DATA(a);
   bdata = DATA(b);
   cdata = DATA(c);

   for (i = 0; i < sm; i ++)
      cdata[i] = adata[i] ^ bdata[i];

   if (sa > sb)
      for (;i < la; i++) cdata[i] = adata[i];
   else
      for (;i < la; i++) cdata[i] = bdata[i];

   STRIP(la, cdata);
   SIZE(c) = la;
}


void _ntl_gor(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   _ntl_gbigint c;
   long sa;
   long sb;
   long sm;
   long la;
   long i;
   long a_alias, b_alias;
   mp_limb_t *adata, *bdata, *cdata;

   if (ZEROP(a)) {
      _ntl_gcopy(b,cc);
      _ntl_gabs(cc);
      return;
   }

   if (ZEROP(b)) {
      _ntl_gcopy(a,cc);
      _ntl_gabs(cc);
      return;
   }

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   sa = SIZE(a);
   if (sa < 0) sa = -sa;

   sb = SIZE(b);
   if (sb < 0) sb = -sb;

   if (sa > sb) {
      la = sa;
      sm = sb;
   } 
   else {
      la = sb;
      sm = sa;
   }

   _ntl_gsetlength(&c, la);
   if (a_alias) a = c;
   if (b_alias) b = c;
   *cc = c;

   adata = DATA(a);
   bdata = DATA(b);
   cdata = DATA(c);

   for (i = 0; i < sm; i ++)
      cdata[i] = adata[i] | bdata[i];

   if (sa > sb)
      for (;i < la; i++) cdata[i] = adata[i];
   else
      for (;i < la; i++) cdata[i] = bdata[i];

   STRIP(la, cdata);
   SIZE(c) = la;
}


void _ntl_gnegate(_ntl_gbigint *aa)
{
   _ntl_gbigint a = *aa;
   if (a) SIZE(a) = -SIZE(a);
}


/*
 * DIRT: this implementation of _ntl_gintoz relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

void _ntl_gintoz(long d, _ntl_gbigint *aa)
{
   _ntl_gbigint a = *aa;

   if (d == 0) {
      if (a) SIZE(a) = 0;
   }
   else if (d > 0) {
      if (!a) {
         _ntl_gsetlength(&a, 1);
         *aa = a;
      }
   
      SIZE(a) = 1;
      DATA(a)[0] = d;
   }
   else {
      if (!a) {
         _ntl_gsetlength(&a, 1);
         *aa = a;
      }
   
      SIZE(a) = -1;
      DATA(a)[0] = -((mp_limb_t) d); /* careful! */
   }
}


/*
 * DIRT: this implementation of _ntl_guintoz relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

void _ntl_guintoz(unsigned long d, _ntl_gbigint *aa)
{
   _ntl_gbigint a = *aa;

   if (d == 0) {
      if (a) SIZE(a) = 0;
   }
   else {
      if (!a) {
         _ntl_gsetlength(&a, 1);
         *aa = a;
      }
   
      SIZE(a) = 1;
      DATA(a)[0] = d;
   }
}


long _ntl_gtoint(_ntl_gbigint a)
{
   unsigned long res = _ntl_gtouint(a);
   return NTL_ULONG_TO_LONG(res);
}

/*
 * DIRT: this implementation of _ntl_gtouint relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

unsigned long _ntl_gtouint(_ntl_gbigint a)
{
   if (ZEROP(a)) 
      return 0;

   if (SIZE(a) > 0) 
      return DATA(a)[0];

   return -DATA(a)[0];
}


long _ntl_gcompare(_ntl_gbigint a, _ntl_gbigint b)
{
   long sa, sb, cmp;
   mp_limb_t *adata, *bdata;

   if (!a) 
      sa = 0;
   else
      sa = SIZE(a);

   if (!b) 
      sb = 0;
   else
      sb = SIZE(b);

   if (sa != sb) {
      if (sa > sb)
         return 1;
      else
         return -1;
   }

   if (sa == 0)
      return 0;

   adata = DATA(a);
   bdata = DATA(b);

   if (sa > 0) {
      cmp = mpn_cmp(adata, bdata, sa);

      if (cmp > 0)
         return 1;
      else if (cmp < 0) 
         return -1;
      else
         return 0;
   }
   else {
      cmp = mpn_cmp(adata, bdata, -sa);

      if (cmp > 0)
         return -1;
      else if (cmp < 0) 
         return 1;
      else
         return 0;
   }
}


long _ntl_gsign(_ntl_gbigint a)
{
   long sa;

   if (!a) return 0;

   sa = SIZE(a);
   if (sa > 0) return 1;
   if (sa == 0) return 0;
   return -1;
}

void _ntl_gabs(_ntl_gbigint *pa)
{
   _ntl_gbigint a = *pa;

   if (!a) return;
   if (SIZE(a) < 0) SIZE(a) = -SIZE(a);
}

long _ntl_gscompare(_ntl_gbigint a, long b)
{
   if (b == 0) {
      long sa;
      if (!a) return 0;
      sa = SIZE(a);
      if (sa > 0) return 1;
      if (sa == 0) return 0;
      return -1;
   }
   else {
      static _ntl_gbigint B = 0;
      _ntl_gintoz(b, &B);
      return _ntl_gcompare(a, B);
   }
}


void _ntl_glshift(_ntl_gbigint n, long k, _ntl_gbigint *rres)
{
   _ntl_gbigint res;
   mp_limb_t *ndata, *resdata, *resdata1;
   long limb_cnt, i, sn, nneg, sres;
   long n_alias;

   if (ZEROP(n)) {
      _ntl_gzero(rres);
      return;
   }

   res = *rres;
   n_alias = (n == res);

   if (!k) {
      if (!n_alias)
         _ntl_gcopy(n, rres);
      return;
   }

   if (k < 0) {
      if (k < -NTL_MAX_LONG) 
         _ntl_gzero(rres);
      else
         _ntl_grshift(n, -k, rres);
      return;
   }

   GET_SIZE_NEG(sn, nneg, n);

   limb_cnt = k/NTL_ZZ_NBITS;
   sres = sn + limb_cnt + 1; 

   if (MustAlloc(res, sres)) {
      _ntl_gsetlength(&res, sres);
      if (n_alias) n = res;
      *rres = res;
   }

   ndata = DATA(n);
   resdata = DATA(res);
   resdata1 = resdata + limb_cnt;
   k %= NTL_ZZ_NBITS;
   sres--;

   if (k != 0) {
      mp_limb_t t = mpn_lshift(resdata1, ndata, sn, k);
      if (t != 0) {
         resdata[sres] = t;
         sres++;
      }
   }
   else {
      for (i = sn-1; i >= 0; i--)
         resdata1[i] = ndata[i];
   }

   for (i = 0; i < limb_cnt; i++)
      resdata[i] = 0;

   if (nneg) sres = -sres;
   SIZE(res) = sres;
}

void _ntl_grshift(_ntl_gbigint n, long k, _ntl_gbigint *rres)
{
   _ntl_gbigint res;
   mp_limb_t *ndata, *resdata, *ndata1;
   long limb_cnt, i, sn, nneg, sres;

   if (ZEROP(n)) {
      _ntl_gzero(rres);
      return;
   }

   if (!k) {
      if (n != *rres)
         _ntl_gcopy(n, rres);
      return;
   }

   if (k < 0) {
      if (k < -NTL_MAX_LONG) ghalt("overflow in _ntl_glshift");
      _ntl_glshift(n, -k, rres);
      return;
   }

   GET_SIZE_NEG(sn, nneg, n);

   limb_cnt = k/NTL_ZZ_NBITS;

   sres = sn - limb_cnt;

   if (sres <= 0) {
      _ntl_gzero(rres);
      return;
   }

   res = *rres;
   if (MustAlloc(res, sres)) {
      /* n won't move if res aliases n */
      _ntl_gsetlength(&res, sres);
      *rres = res;
   }

   ndata = DATA(n);
   resdata = DATA(res);
   ndata1 = ndata + limb_cnt;
   k %= NTL_ZZ_NBITS;

   if (k != 0) {
      mpn_rshift(resdata, ndata1, sres, k);
      if (resdata[sres-1] == 0)
         sres--;
   }
   else {
      for (i = 0; i < sres; i++)
         resdata[i] = ndata1[i];
   }

   if (nneg) sres = -sres;
   SIZE(res) = sres;
}
   
   

void
_ntl_gadd(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   long sa, aneg, sb, bneg, sc, cmp;
   mp_limb_t *adata, *bdata, *cdata, carry;
   _ntl_gbigint c;
   long a_alias, b_alias;

   if (ZEROP(a)) {
      _ntl_gcopy(b, cc);
      return;
   }

   if (ZEROP(b)) {
      _ntl_gcopy(a, cc);
      return;
   }

   GET_SIZE_NEG(sa, aneg, a);
   GET_SIZE_NEG(sb, bneg, b);

   if (sa < sb) {
      SWAP_BIGINT(a, b);
      SWAP_LONG(sa, sb);
      SWAP_LONG(aneg, bneg);
   }

   /* sa >= sb */

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   if (aneg == bneg) {
      /* same sign => addition */

      sc = sa + 1;
      if (MustAlloc(c, sc)) {
         _ntl_gsetlength(&c, sc);
         if (a_alias) a = c; 
         if (b_alias) b = c;
         *cc = c;
      }

      adata = DATA(a);
      bdata = DATA(b);
      cdata = DATA(c);

      carry = mpn_add(cdata, adata, sa, bdata, sb);
      if (carry) 
         cdata[sc-1] = carry;
      else
         sc--;

      if (aneg) sc = -sc;
      SIZE(c) = sc;
   }
   else {
      /* opposite sign => subtraction */

      sc = sa;
      if (MustAlloc(c, sc)) {
         _ntl_gsetlength(&c, sc);
         if (a_alias) a = c; 
         if (b_alias) b = c;
         *cc = c;
      }

      adata = DATA(a);
      bdata = DATA(b);
      cdata = DATA(c);

      if (sa > sb) 
         cmp = 1;
      else
         cmp = mpn_cmp(adata, bdata, sa);

      if (cmp == 0) {
         SIZE(c) = 0;
      }
      else {
         if (cmp < 0) cmp = 0;
         if (cmp > 0) cmp = 1;
         /* abs(a) != abs(b) && (abs(a) > abs(b) <=> cmp) */

         if (cmp)
            mpn_sub(cdata, adata, sa, bdata, sb);
         else
            mpn_sub(cdata, bdata, sb, adata, sa); /* sa == sb */

         STRIP(sc, cdata);
         if (aneg == cmp) sc = -sc;
         SIZE(c) = sc;
      }
   }
}

void
_ntl_gsadd(_ntl_gbigint a, long b, _ntl_gbigint *cc)
{
   static _ntl_gbigint B = 0;
   _ntl_gintoz(b, &B);
   _ntl_gadd(a, B, cc);
}

void
_ntl_gsub(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   long sa, aneg, sb, bneg, sc, cmp, rev;
   mp_limb_t *adata, *bdata, *cdata, carry;
   _ntl_gbigint c;
   long a_alias, b_alias;

   if (ZEROP(a)) {
      _ntl_gcopy(b, cc);
      c = *cc;
      if (c) SIZE(c) = -SIZE(c); 
      return;
   }

   if (ZEROP(b)) {
      _ntl_gcopy(a, cc);
      return;
   }

   GET_SIZE_NEG(sa, aneg, a);
   GET_SIZE_NEG(sb, bneg, b);

   if (sa < sb) {
      SWAP_BIGINT(a, b);
      SWAP_LONG(sa, sb);
      SWAP_LONG(aneg, bneg);
      rev = 1;
   }
   else 
      rev = 0;

   /* sa >= sb */

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   if (aneg != bneg) {
      /* opposite sign => addition */

      sc = sa + 1;
      if (MustAlloc(c, sc)) {
         _ntl_gsetlength(&c, sc);
         if (a_alias) a = c; 
         if (b_alias) b = c;
         *cc = c;
      }

      adata = DATA(a);
      bdata = DATA(b);
      cdata = DATA(c);

      carry = mpn_add(cdata, adata, sa, bdata, sb);
      if (carry) 
         cdata[sc-1] = carry;
      else
         sc--;

      if (aneg ^ rev) sc = -sc;
      SIZE(c) = sc;
   }
   else {
      /* same sign => subtraction */

      sc = sa;
      if (MustAlloc(c, sc)) {
         _ntl_gsetlength(&c, sc);
         if (a_alias) a = c; 
         if (b_alias) b = c;
         *cc = c;
      }

      adata = DATA(a);
      bdata = DATA(b);
      cdata = DATA(c);

      if (sa > sb) 
         cmp = 1;
      else
         cmp = mpn_cmp(adata, bdata, sa);

      if (cmp == 0) {
         SIZE(c) = 0;
      }
      else {
         if (cmp < 0) cmp = 0;
         if (cmp > 0) cmp = 1;
         /* abs(a) != abs(b) && (abs(a) > abs(b) <=> cmp) */

         if (cmp)
            mpn_sub(cdata, adata, sa, bdata, sb);
         else
            mpn_sub(cdata, bdata, sb, adata, sa); /* sa == sb */

         STRIP(sc, cdata);
         if ((aneg == cmp) ^ rev) sc = -sc;
         SIZE(c) = sc;
      }
   }
}

void
_ntl_gsubpos(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   long sa, sb, sc;
   mp_limb_t *adata, *bdata, *cdata;
   _ntl_gbigint c;
   long a_alias, b_alias;

   if (ZEROP(a)) {
      _ntl_gzero(cc);
      return;
   }

   if (ZEROP(b)) {
      _ntl_gcopy(a, cc);
      return;
   }

   sa = SIZE(a);
   sb = SIZE(b);

   c = *cc;
   a_alias = (a == c);
   b_alias = (b == c);

   sc = sa;
   if (MustAlloc(c, sc)) {
      _ntl_gsetlength(&c, sc);
      if (a_alias) a = c; 
      if (b_alias) b = c;
      *cc = c;
   }

   adata = DATA(a);
   bdata = DATA(b);
   cdata = DATA(c);

   mpn_sub(cdata, adata, sa, bdata, sb);

   STRIP(sc, cdata);
   SIZE(c) = sc;
}

void _ntl_gmul(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *cc)
{
   static _ntl_gbigint mem = 0;

   long sa, aneg, sb, bneg, alias, sc;
   mp_limb_t *adata, *bdata, *cdata, msl;
   _ntl_gbigint c;

   if (ZEROP(a) || ZEROP(b)) {
      _ntl_gzero(cc);
      return;
   }

   GET_SIZE_NEG(sa, aneg, a);
   GET_SIZE_NEG(sb, bneg, b);

   if (a == *cc || b == *cc) {
      c = mem;
      alias = 1;
   }
   else {
      c = *cc;
      alias = 0;
   }

   sc = sa + sb;
   if (MustAlloc(c, sc))
      _ntl_gsetlength(&c, sc);

   if (alias)
      mem = c;
   else
      *cc = c;

   adata = DATA(a);
   bdata = DATA(b);
   cdata = DATA(c);

   if (sa >= sb)
      msl = mpn_mul(cdata, adata, sa, bdata, sb);
   else
      msl = mpn_mul(cdata, bdata, sb, adata, sa);

   if (!msl) sc--;
   if (aneg != bneg) sc = -sc;
   SIZE(c) = sc;

   if (alias) _ntl_gcopy(mem, cc);
}

void _ntl_gsq(_ntl_gbigint a, _ntl_gbigint *cc)
{
   _ntl_gmul(a, a, cc);
   /* this is good enough...eventually, mpn_sqr_n will be called */
}


/*
 * DIRT: this implementation of _ntl_gsmul relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

void
_ntl_gsmul(_ntl_gbigint a, long d, _ntl_gbigint *bb)
{
   long sa, sb;
   long anegative, bnegative;
   _ntl_gbigint b;
   mp_limb_t *adata, *bdata;
   mp_limb_t dd, carry;
   long a_alias;

   if (ZEROP(a) || !d) {
      _ntl_gzero(bb);
      return;
   }

   GET_SIZE_NEG(sa, anegative, a);

   if (d < 0) {
      dd = - ((mp_limb_t) d); /* careful ! */
      bnegative = 1-anegative;
   }
   else {
      dd = (mp_limb_t) d;
      bnegative = anegative;
   }

   sb = sa + 1;

   b = *bb;
   a_alias = (a == b);

   if (MustAlloc(b, sb)) {
      _ntl_gsetlength(&b, sb);
      if (a_alias) a = b;
      *bb = b;
   }

   adata = DATA(a);
   bdata = DATA(b);

   if (dd == 2) 
      carry = mpn_lshift(bdata, adata, sa, 1);
   else
      carry = mpn_mul_1(bdata, adata, sa, dd);

   if (carry) 
      bdata[sa] = carry;
   else
      sb--;

   if (bnegative) sb = -sb;
   SIZE(b) = sb;
}

/*
 * DIRT: this implementation of _ntl_gsdiv relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

long _ntl_gsdiv(_ntl_gbigint a, long d, _ntl_gbigint *bb)
{
   long sa, aneg, sb, dneg;
   _ntl_gbigint b;
   mp_limb_t dd, *adata, *bdata;
   long r;

   if (!d) {
      ghalt("division by zero in _ntl_gsdiv");
   }

   if (ZEROP(a)) {
      _ntl_gzero(bb);
      return (0);
   }

   GET_SIZE_NEG(sa, aneg, a);

   if (d < 0) {
      dd = - ((mp_limb_t) d); /* careful ! */
      dneg = 1;
   }
   else {
      dd = (mp_limb_t) d;
      dneg = 0;
   }

   sb = sa;
   b = *bb;
   if (MustAlloc(b, sb)) {
      /* if b aliases a, then b won't move */
      _ntl_gsetlength(&b, sb);
      *bb = b;
   }

   adata = DATA(a);
   bdata = DATA(b);

   if (dd == 2) 
      r = mpn_rshift(bdata, adata, sa, 1) >> (NTL_ZZ_NBITS - 1);
   else
      r = mpn_divmod_1(bdata, adata, sa, dd);

   if (bdata[sb-1] == 0)
      sb--;

   SIZE(b) = sb;

   if (aneg || dneg) {
      if (aneg != dneg) {
         if (!r) {
            SIZE(b) = -SIZE(b);
         }
         else {
            _ntl_gsadd(b, 1, &b);
            SIZE(b) = -SIZE(b);
            if (dneg)
               r = r + d;
            else
               r = d - r;
            *bb = b;
         }
      }
      else
         r = -r;
   }

   return r;
}

/*
 * DIRT: this implementation of _ntl_gsmod relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

long _ntl_gsmod(_ntl_gbigint a, long d)
{
   long sa, aneg, dneg;
   mp_limb_t dd, *adata;
   long r;

   if (!d) {
      ghalt("division by zero in _ntl_gsmod");
   }

   if (ZEROP(a)) {
      return (0);
   }

   GET_SIZE_NEG(sa, aneg, a);

   if (d < 0) {
      dd = - ((mp_limb_t) d); /* careful ! */
      dneg = 1;
   }
   else {
      dd = (mp_limb_t) d;
      dneg = 0;
   }

   adata = DATA(a);

   if (dd == 2) 
      r = adata[0] & 1;
   else
      r = mpn_mod_1(adata, sa, dd);

   if (aneg || dneg) {
      if (aneg != dneg) {
         if (r) {
            if (dneg)
               r = r + d;
            else
               r = d - r;
         }
      }
      else
         r = -r;
   }

   return r;
}


void _ntl_gdiv(_ntl_gbigint a, _ntl_gbigint d, 
               _ntl_gbigint *bb, _ntl_gbigint *rr)
{
   static _ntl_gbigint b = 0, rmem = 0;

   _ntl_gbigint r;

   long sa, aneg, sb, sd, dneg, sr, in_place;
   mp_limb_t *adata, *ddata, *bdata, *rdata;

   if (ZEROP(d)) {
      ghalt("division by zero in _ntl_gdiv");
   }

   if (ZEROP(a)) {
      if (bb) _ntl_gzero(bb);
      if (rr) _ntl_gzero(rr);
      return;
   }

   GET_SIZE_NEG(sa, aneg, a);
   GET_SIZE_NEG(sd, dneg, d);

   if (!aneg && !dneg && rr && *rr != a && *rr != d) {
      in_place = 1;
      r = *rr;
   }
   else {
      in_place = 0;
      r = rmem;
   }

   if (sa < sd) {
      _ntl_gzero(&b);
      _ntl_gcopy(a, &r);
      if (aneg) SIZE(r) = -SIZE(r);
      goto done;
   }

   sb = sa-sd+1;
   if (MustAlloc(b, sb))
      _ntl_gsetlength(&b, sb);

   sr = sd;
   if (MustAlloc(r, sr))
      _ntl_gsetlength(&r, sr);

   adata = DATA(a);
   ddata = DATA(d);
   bdata = DATA(b);
   rdata = DATA(r);

   mpn_tdiv_qr(bdata, rdata, 0, adata, sa, ddata, sd);

   if (bdata[sb-1] == 0)
      sb--;
   SIZE(b) = sb;

   STRIP(sr, rdata);
   SIZE(r) = sr;

done:

   if (aneg || dneg) {
      if (aneg != dneg) {
         if (ZEROP(r)) {
            SIZE(b) = -SIZE(b);
         }
         else {
            if (bb) {
               _ntl_gsadd(b, 1, &b);
               SIZE(b) = -SIZE(b);
            }
            if (rr) {
               if (dneg)
                  _ntl_gadd(r, d, &r);
               else
                  _ntl_gsub(d, r, &r);
            }
         }
      }
      else
         SIZE(r) = -SIZE(r);
   }

   if (bb) _ntl_gcopy(b, bb);

   if (in_place)
      *rr = r;
   else {
      if (rr) _ntl_gcopy(r, rr);
      rmem = r;
   }
}


/* a simplified mod operation:  assumes a >= 0, d > 0 are non-negative,
 * that space for the result has already been allocated,
 * and that inputs do not alias output. */

static
void gmod_simple(_ntl_gbigint a, _ntl_gbigint d, _ntl_gbigint *rr)
{
   static _ntl_gbigint b = 0;

   long sa, sb, sd, sr;
   mp_limb_t *adata, *ddata, *bdata, *rdata;
   _ntl_gbigint r;

   if (ZEROP(a)) {
      _ntl_gzero(rr);
      return;
   }

   sa = SIZE(a);
   sd = SIZE(d);

   if (sa < sd) {
      _ntl_gcopy(a, rr);
      return;
   }

   sb = sa-sd+1;
   if (MustAlloc(b, sb))
      _ntl_gsetlength(&b, sb);

   sr = sd;
   r = *rr;

   adata = DATA(a);
   ddata = DATA(d);
   bdata = DATA(b);
   rdata = DATA(r);

   mpn_tdiv_qr(bdata, rdata, 0, adata, sa, ddata, sd);

   STRIP(sr, rdata);
   SIZE(r) = sr;
}


void _ntl_gmod(_ntl_gbigint a, _ntl_gbigint d, _ntl_gbigint *rr)
{
   _ntl_gdiv(a, d, 0, rr);
}

void _ntl_gquickmod(_ntl_gbigint *rr, _ntl_gbigint d)
{
   _ntl_gdiv(*rr, d, 0, rr);
}

void _ntl_gsqrt(_ntl_gbigint n, _ntl_gbigint *rr)
{
   static _ntl_gbigint r = 0;

   long sn, sr;
   mp_limb_t *ndata, *rdata;

   if (ZEROP(n)) {
      _ntl_gzero(rr);
      return;
   }

   sn = SIZE(n);
   if (sn < 0) ghalt("negative argument to _ntl_sqrt");

   sr = (sn+1)/2;
   _ntl_gsetlength(&r, sr);

   ndata = DATA(n);
   rdata = DATA(r);

   mpn_sqrtrem(rdata, 0, ndata, sn);

   STRIP(sr, rdata);
   SIZE(r) = sr;

   _ntl_gcopy(r, rr);
}

/*
 * DIRT: this implementation of _ntl_gsqrts relies crucially
 * on the assumption that the number of bits per limb_t is at least
 * equal to the number of bits per long.
 */

long _ntl_gsqrts(long n)
{
   mp_limb_t ndata, rdata;

   if (n == 0) {
      return 0;
   }

   if (n < 0) ghalt("negative argument to _ntl_sqrts");

   ndata = n;

   mpn_sqrtrem(&rdata, 0, &ndata, 1);

   return rdata;
}


void _ntl_ggcd(_ntl_gbigint m1, _ntl_gbigint m2, _ntl_gbigint *r)
{
   static _ntl_gbigint s1 = 0, s2 = 0, res = 0;

   long k1, k2, k_min, l1, l2, ss1, ss2, sres;

   _ntl_gcopy(m1, &s1);
   _ntl_gabs(&s1);

   _ntl_gcopy(m2, &s2);
   _ntl_gabs(&s2);

   if (ZEROP(s1)) {
      _ntl_gcopy(s2, r);
      return;
   }

   if (ZEROP(s2)) {
      _ntl_gcopy(s1, r);
      return;
   }

   k1 = _ntl_gmakeodd(&s1);
   k2 = _ntl_gmakeodd(&s2);

   if (k1 <= k2)
      k_min = k1;
   else
      k_min = k2;

   l1 = _ntl_g2log(s1);
   l2 = _ntl_g2log(s2);

   ss1 = SIZE(s1);
   ss2 = SIZE(s2);

   if (ss1 >= ss2)
      sres = ss1;
   else
      sres = ss2;

   /* set to max: gmp documentation is unclear on this point */

   _ntl_gsetlength(&res, sres);
   
   if (l1 >= l2)
      SIZE(res) = mpn_gcd(DATA(res), DATA(s1), ss1, DATA(s2), ss2);
   else
      SIZE(res) = mpn_gcd(DATA(res), DATA(s2), ss2, DATA(s1), ss1);

   _ntl_glshift(res, k_min, &res);

   _ntl_gcopy(res, r);
}

static long 
gxxeucl(
   _ntl_gbigint ain,
   _ntl_gbigint nin,
   _ntl_gbigint *invv,
   _ntl_gbigint *uu
   )
{
   static _ntl_gbigint a = 0;
   static _ntl_gbigint n = 0;
   static _ntl_gbigint q = 0;
   static _ntl_gbigint w = 0;
   static _ntl_gbigint x = 0;
   static _ntl_gbigint y = 0;
   static _ntl_gbigint z = 0;
   _ntl_gbigint inv = *invv;
   _ntl_gbigint u = *uu;
   long diff;
   long ilo;
   long sa;
   long sn;
   long temp;
   long e;
   long fast;
   long parity;
   long gotthem;
   mp_limb_t *p;
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

   _ntl_gsetlength(&a, (e = (SIZE(ain) > SIZE(nin) ? SIZE(ain) : SIZE(nin))));
   _ntl_gsetlength(&n, e);
   _ntl_gsetlength(&q, e);
   _ntl_gsetlength(&w, e);
   _ntl_gsetlength(&x, e);
   _ntl_gsetlength(&y, e);
   _ntl_gsetlength(&z, e);
   _ntl_gsetlength(&inv, e);
   *invv = inv;
   _ntl_gsetlength(&u, e);
   *uu = u;

   fhi1 = 1.0 + ((double) 32.0)/NTL_FDOUBLE_PRECISION;
   flo1 = 1.0 - ((double) 32.0)/NTL_FDOUBLE_PRECISION;

   fhi = 1.0 + ((double) 8.0)/NTL_FDOUBLE_PRECISION;
   flo = 1.0 - ((double) 8.0)/NTL_FDOUBLE_PRECISION;

   _ntl_gcopy(ain, &a);
   _ntl_gcopy(nin, &n);

   _ntl_gone(&inv);
   _ntl_gzero(&w);

   while (SIZE(n) > 0)
   {
      gotthem = 0;
      sa = SIZE(a);
      sn = SIZE(n);
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = SIZE(a);
         p = DATA(a) + (sa-1);
         num = (double) (*p) * NTL_ZZ_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_ZZ_FRADIX;
         if (sa > 2)
            num += (*(p - 1));

         sn = SIZE(n);
         p = DATA(n) + (sn-1);
         den = (double) (*p) * NTL_ZZ_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_ZZ_FRADIX;
         if (sn > 2)
            den += (*(p - 1));

         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_ZZ_FRADIX;
            lo *= NTL_ZZ_FRADIX;
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
            if (hi >= NTL_SP_BOUND)
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
                     hi = NTL_SP_BOUND;
                  temp = try11;
                  try11 = try21;
                  if ((NTL_WSP_BOUND - temp) / ilo < try21)
                     fast = 0;
                  else
                     try21 = temp + ilo * try21;
                  temp = try12;
                  try12 = try22;
                  if ((NTL_WSP_BOUND - temp) / ilo < try22)
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
         _ntl_gsmul(inv, got11, &x);
         _ntl_gsmul(w, got12, &y);
         _ntl_gsmul(inv, got21, &z);
         _ntl_gsmul(w, got22, &w);
         _ntl_gadd(x, y, &inv);
         _ntl_gadd(z, w, &w);
         _ntl_gsmul(a, got11, &x);
         _ntl_gsmul(n, got12, &y);
         _ntl_gsmul(a, got21, &z);
         _ntl_gsmul(n, got22, &n);
         _ntl_gsub(x, y, &a);
         _ntl_gsub(n, z, &n);
      }
      else
      {
         _ntl_gdiv(a, n, &q, &a);
         _ntl_gmul(q, w, &x);
         _ntl_gadd(inv, x, &inv);
         if (!ZEROP(a))
         {
            _ntl_gdiv(n, a, &q, &n);
            _ntl_gmul(q, inv, &x);
            _ntl_gadd(w, x, &w);
         }
         else
         {
            _ntl_gcopy(n, &a);
            _ntl_gzero(&n);
            _ntl_gcopy(w, &inv);
            _ntl_gnegate(&inv);
         }
      }
   }

   if (_ntl_gscompare(a, 1) == 0)
      e = 0;
   else 
      e = 1;

   _ntl_gcopy(a, &u);

   *invv = inv;
   *uu = u;

   return (e);
}

#if 0
void
_ntl_gexteucl(
	_ntl_gbigint aa,
	_ntl_gbigint *xa,
	_ntl_gbigint bb,
	_ntl_gbigint *xb,
	_ntl_gbigint *d
	)
{
   static _ntl_gbigint modcon = 0;
   static _ntl_gbigint a=0, b=0;
   long anegative = 0;
   long bnegative = 0;

   _ntl_gcopy(aa, &a);
   _ntl_gcopy(bb, &b);

   if (a && SIZE(a) < 0) {
      anegative = 1;
      SIZE(a) = -SIZE(a);
   }
   else
      anegative = 0;

   if (b && SIZE(b) < 0) {
      bnegative = 1;
      SIZE(b) = -SIZE(b);
   }
   else
      bnegative = 0;


   if (ZEROP(b))
   {
      _ntl_gone(xa);
      _ntl_gzero(xb);
      _ntl_gcopy(a, d);
      goto done;
   }

   if (ZEROP(a))
   {
      _ntl_gzero(xa);
      _ntl_gone(xb);
      _ntl_gcopy(b, d);
      goto done;
   }

   gxxeucl(a, b, xa, d);
   _ntl_gmul(a, *xa, xb);
   _ntl_gsub(*d, *xb, xb);
   _ntl_gdiv(*xb, b, xb, &modcon);

   if (!ZEROP(modcon))
   {
      ghalt("non-zero remainder in _ntl_gexteucl   BUG");
   }


done:
   if (anegative)
   {
      _ntl_gnegate(xa);
   }
   if (bnegative)
   {
      _ntl_gnegate(xb);
   }
}
#endif

void
_ntl_gexteucl(
	_ntl_gbigint ain,
	_ntl_gbigint *xap,
	_ntl_gbigint bin,
	_ntl_gbigint *xbp,
	_ntl_gbigint *dp
	)
{
   if (ZEROP(bin)) {
      long asign = _ntl_gsign(ain);

      _ntl_gcopy(ain, dp);
      _ntl_gabs(dp);
      _ntl_gintoz( (asign >= 0 ? 1 : -1), xap);
      _ntl_gzero(xbp);
   }
   else if (ZEROP(ain)) {
      long bsign = _ntl_gsign(bin);

      _ntl_gcopy(bin, dp);
      _ntl_gabs(dp);
      _ntl_gzero(xap);
      _ntl_gintoz(bsign, xbp); 
   }
   else {
      static _ntl_gbigint a = 0, b = 0, xa = 0, xb = 0, d = 0, tmp = 0;

      long sa, aneg, sb, bneg, rev;
      mp_limb_t *adata, *bdata, *ddata, *xadata;
      mp_size_t sxa, sd;

      GET_SIZE_NEG(sa, aneg, ain);
      GET_SIZE_NEG(sb, bneg, bin);

      _ntl_gsetlength(&a, sa+1); /* +1 because mpn_gcdext may need it */
      _ntl_gcopy(ain, &a);

      _ntl_gsetlength(&b, sb+1); /* +1 because mpn_gcdext may need it */
      _ntl_gcopy(bin, &b);


      adata = DATA(a);
      bdata = DATA(b);

      if (sa < sb || (sa == sb && mpn_cmp(adata, bdata, sa) < 0)) {
         SWAP_BIGINT(ain, bin);
         SWAP_LONG(sa, sb);
         SWAP_LONG(aneg, bneg);
         SWAP_LIMB_PTR(adata, bdata);
         rev = 1;
      }
      else 
         rev = 0;

      _ntl_gsetlength(&d, sa+1); /* +1 because mpn_gcdext may need it...
                                    documentation is unclear, but this is
                                    what is done in mpz_gcdext */
      _ntl_gsetlength(&xa, sa+1); /* ditto */

      ddata = DATA(d);
      xadata = DATA(xa);

      sd = mpn_gcdext(ddata, xadata, &sxa, adata, sa, bdata, sb);

      SIZE(d) = sd;
      SIZE(xa) = sxa;

      /* Thes two ForceNormal's are work-arounds for GMP bugs 
         in GMP 4.3.0 */
      ForceNormal(d);
      ForceNormal(xa);

      /* now we normalize xa, so that so that xa in ( -b/2d, b/2d ],
         which makes the output agree with Euclid's algorithm,
         regardless of what mpn_gcdext does */

      if (!ZEROP(xa)) {
         _ntl_gcopy(bin, &b);
         SIZE(b) = sb;
         if (!ONEP(d)) {
            _ntl_gdiv(b, d, &b, &tmp);
            if (!ZEROP(tmp)) ghalt("internal bug in _ntl_gexteucl");
         }

         if (SIZE(xa) > 0) { /* xa positive */
            if (_ntl_gcompare(xa, b) > 0) { 
               _ntl_gmod(xa, b, &xa);
            }
            _ntl_glshift(xa, 1, &tmp);
            if (_ntl_gcompare(tmp, b) > 0) {
               _ntl_gsub(xa, b, &xa);
            }
         }
         else { /* xa negative */
            SIZE(xa) = -SIZE(xa);
            if (_ntl_gcompare(xa, b) > 0) {
               SIZE(xa) = -SIZE(xa);
               _ntl_gmod(xa, b, &xa);
               _ntl_gsub(xa, b, &xa);
            }
            else {
               SIZE(xa) = -SIZE(xa);
            }
            _ntl_glshift(xa, 1, &tmp);
            SIZE(tmp) = -SIZE(tmp);
            if (_ntl_gcompare(tmp, b) >= 0) {
               _ntl_gadd(xa, b, &xa);
            }
         }
      }

      /* end normalize */
    

      if (aneg) _ntl_gnegate(&xa);

      _ntl_gmul(ain, xa, &tmp);
      _ntl_gsub(d, tmp, &tmp);
      _ntl_gdiv(tmp, bin, &xb, &tmp);

      if (!ZEROP(tmp)) ghalt("internal bug in _ntl_gexteucl");

      if (rev) SWAP_BIGINT(xa, xb);

      _ntl_gcopy(xa, xap);
      _ntl_gcopy(xb, xbp);
      _ntl_gcopy(d, dp); 
   }
}


long _ntl_ginv(_ntl_gbigint ain, _ntl_gbigint nin, _ntl_gbigint *invv)
{
   static _ntl_gbigint u = 0;
   static _ntl_gbigint d = 0;
   static _ntl_gbigint a = 0;
   static _ntl_gbigint n = 0;

   long sz; 
   long sd;
   mp_size_t su;

   if (_ntl_gscompare(nin, 1) <= 0) {
      ghalt("InvMod: second input <= 1");
   }

   if (_ntl_gsign(ain) < 0) {
      ghalt("InvMod: first input negative");
   }

   if (_ntl_gcompare(ain, nin) >= 0) {
      ghalt("InvMod: first input too big");
   }

   sz = SIZE(nin) + 2;

   if (MustAlloc(a, sz))
      _ntl_gsetlength(&a, sz);


   if (MustAlloc(n, sz))
       _ntl_gsetlength(&n, sz);


   if (MustAlloc(d, sz))
       _ntl_gsetlength(&d, sz);

   if (MustAlloc(u, sz))
       _ntl_gsetlength(&u, sz);

   _ntl_gadd(ain, nin, &a);
   _ntl_gcopy(nin, &n);

   /* We apply mpn_gcdext to (a, n) = (ain+nin, nin), because that function
    * only computes the co-factor of the larger input. This way, we avoid
    * a multiplication and a division.
    */

   sd = mpn_gcdext(DATA(d), DATA(u), &su, DATA(a), SIZE(a), DATA(n), SIZE(n));

   SIZE(d) = sd;
   SIZE(u) = su;

      /* Thes two ForceNormal's are work-arounds for GMP bugs 
         in GMP 4.3.0 */
      ForceNormal(d);
      ForceNormal(u);


   if (ONEP(d)) {

      /*
       * We make sure that u is in range 0..n-1, just in case
       * GMP is sloppy.
       */


      if (_ntl_gsign(u) < 0) {
         _ntl_gadd(u, nin, &u);
         if (_ntl_gsign(u) < 0) {
            _ntl_gmod(u, nin, &u);
         }
      }
      else if (_ntl_gcompare(u, nin) >= 0) {
         _ntl_gsub(u, nin, &u);
         if (_ntl_gcompare(u, nin) >= 0) {
             _ntl_gmod(u, nin, &u);
         }
      }

      _ntl_gcopy(u, invv);
      return 0;
   }
   else {
      _ntl_gcopy(d, invv);
      return 1;
   }
}


void
_ntl_ginvmod(
	_ntl_gbigint a,
	_ntl_gbigint n,
	_ntl_gbigint *c
	)
{
	if (_ntl_ginv(a, n, c))
		ghalt("undefined inverse in _ntl_ginvmod");
}


void
_ntl_gaddmod(
	_ntl_gbigint a,
	_ntl_gbigint b,
	_ntl_gbigint n,
	_ntl_gbigint *c
	)
{
	if (*c != n) {
		_ntl_gadd(a, b, c);
		if (_ntl_gcompare(*c, n) >= 0)
			_ntl_gsubpos(*c, n, c);
	}
	else {
		static _ntl_gbigint mem = 0;

		_ntl_gadd(a, b, &mem);
		if (_ntl_gcompare(mem, n) >= 0)
			_ntl_gsubpos(mem, n, c);
		else
			_ntl_gcopy(mem, c);
	}
}


void
_ntl_gsubmod(
	_ntl_gbigint a,
	_ntl_gbigint b,
	_ntl_gbigint n,
	_ntl_gbigint *c
	)
{
	static _ntl_gbigint mem = 0;
	long cmp;

	if ((cmp=_ntl_gcompare(a, b)) < 0) {
		_ntl_gadd(n, a, &mem);
		_ntl_gsubpos(mem, b, c);
	} else if (!cmp) 
		_ntl_gzero(c);
	else 
		_ntl_gsubpos(a, b, c);
}

void
_ntl_gsmulmod(
	_ntl_gbigint a,
	long d,
	_ntl_gbigint n,
	_ntl_gbigint *c
	)
{
	static _ntl_gbigint mem = 0;

	_ntl_gsmul(a, d, &mem);
	_ntl_gmod(mem, n, c);
}



void
_ntl_gmulmod(
	_ntl_gbigint a,
	_ntl_gbigint b,
	_ntl_gbigint n,
	_ntl_gbigint *c
	)
{
	static _ntl_gbigint mem = 0;

	_ntl_gmul(a, b, &mem);
	_ntl_gmod(mem, n, c);
}

void
_ntl_gsqmod(
	_ntl_gbigint a,
	_ntl_gbigint n,
	_ntl_gbigint *c
	)
{
	_ntl_gmulmod(a, a, n, c);
}


double _ntl_gdoub_aux(_ntl_gbigint n)
{
   double res;
   mp_limb_t *ndata;
   long i, sn, nneg;

   if (!n)
      return ((double) 0);

   GET_SIZE_NEG(sn, nneg, n);

   ndata = DATA(n);

   res = 0;
   for (i = sn-1; i >= 0; i--)
      res = res * NTL_ZZ_FRADIX + ((double) ndata[i]);

   if (nneg) res = -res;

   return res;
}

long _ntl_ground_correction(_ntl_gbigint a, long k, long residual)
{
   long direction;
   long p;
   long sgn;
   long bl;
   mp_limb_t wh;
   long i;
   mp_limb_t *adata;

   if (SIZE(a) > 0)
      sgn = 1;
   else
      sgn = -1;

   adata = DATA(a);

   p = k - 1;
   bl = (p/NTL_ZZ_NBITS);
   wh = ((mp_limb_t) 1) << (p - NTL_ZZ_NBITS*bl);

   if (adata[bl] & wh) {
      /* bit is 1...we have to see if lower bits are all 0
         in order to implement "round to even" */

      if (adata[bl] & (wh - ((mp_limb_t) 1))) 
         direction = 1;
      else {
         i = bl - 1;
         while (i >= 0 && adata[i] == 0) i--;
         if (i >= 0)
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

         /*
          * DIRT: if GMP has non-empty "nails", this won't work.
          */

         if (wh == 0) {
            wh = 1;
            bl++;
         }

         if (adata[bl] & wh)
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




double _ntl_gdoub(_ntl_gbigint n)
{
   static _ntl_gbigint tmp = 0;

   long s;
   long shamt;
   long correction;
   double x;

   s = _ntl_g2log(n);
   shamt = s - NTL_DOUBLE_PRECISION;

   if (shamt <= 0)
      return _ntl_gdoub_aux(n);

   _ntl_grshift(n, shamt, &tmp);

   correction = _ntl_ground_correction(n, shamt, 0);

   if (correction) _ntl_gsadd(tmp, correction, &tmp);

   x = _ntl_gdoub_aux(tmp);

   x = _ntl_ldexp(x, shamt);

   return x;
}


double _ntl_glog(_ntl_gbigint n)
{
   static _ntl_gbigint tmp = 0;

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

   if (_ntl_gsign(n) <= 0)
      ghalt("log argument <= 0");

   s = _ntl_g2log(n);
   shamt = s - NTL_DOUBLE_PRECISION;

   if (shamt <= 0)
      return log(_ntl_gdoub_aux(n));

   _ntl_grshift(n, shamt, &tmp);

   correction = _ntl_ground_correction(n, shamt, 0);

   if (correction) _ntl_gsadd(tmp, correction, &tmp);

   x = _ntl_gdoub_aux(tmp);

   return log(x) + shamt*log_2;
}





/* To implement _ntl_gdoubtoz, I've implemented essentially the
 * same algorithm as in LIP, processing in blocks of
 * NTL_SP_NBITS bits, rather than NTL_ZZ_NBITS.
 * This is conversion is rather delicate, and I don't want to
 * make any new assumptions about the underlying arithmetic.
 * This implementation should be quite portable. */

void _ntl_gdoubtoz(double a, _ntl_gbigint *xx)
{
   static _ntl_gbigint x = 0;

   long neg, i, t, sz;

   a = floor(a);

   if (!_ntl_IsFinite(&a))
      ghalt("_ntl_gdoubtoz: attempt to convert non-finite value");

   if (a < 0) {
      a = -a;
      neg = 1;
   }
   else
      neg = 0;

   if (a == 0) {
      _ntl_gzero(xx);
      return;
   }

   sz = 0;
   while (a >= 1) {
      a = a*(1.0/NTL_SP_FBOUND);
      sz++;
   }

   i = 0;
   _ntl_gzero(&x);

   while (a != 0) {
      i++;
      a = a*NTL_SP_FBOUND;
      t = (long) a;
      a = a - t;

      if (i == 1) {
         _ntl_gintoz(t, &x);
      }
      else {
         _ntl_glshift(x, NTL_SP_NBITS, &x);
         _ntl_gsadd(x, t, &x);
      }
   }

   if (i > sz) ghalt("bug in _ntl_gdoubtoz");

   _ntl_glshift(x, (sz-i)*NTL_SP_NBITS, xx);
   if (neg) _ntl_gnegate(xx);
}



/* I've adapted LIP's extended euclidean algorithm to
 * do rational reconstruction.  -- VJS.
 */


long 
_ntl_gxxratrecon(
   _ntl_gbigint ain,
   _ntl_gbigint nin,
   _ntl_gbigint num_bound,
   _ntl_gbigint den_bound,
   _ntl_gbigint *num_out,
   _ntl_gbigint *den_out
   )
{
   static _ntl_gbigint a = 0;
   static _ntl_gbigint n = 0;
   static _ntl_gbigint q = 0;
   static _ntl_gbigint w = 0;
   static _ntl_gbigint x = 0;
   static _ntl_gbigint y = 0;
   static _ntl_gbigint z = 0;
   static _ntl_gbigint inv = 0;
   static _ntl_gbigint u = 0;
   static _ntl_gbigint a_bak = 0;
   static _ntl_gbigint n_bak = 0;
   static _ntl_gbigint inv_bak = 0;
   static _ntl_gbigint w_bak = 0;

   mp_limb_t *p;

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

   if (_ntl_gsign(num_bound) < 0)
      ghalt("rational reconstruction: bad numerator bound");

   if (!num_bound)
      snum = 0;
   else
      snum = SIZE(num_bound);

   if (_ntl_gsign(den_bound) <= 0)
      ghalt("rational reconstruction: bad denominator bound");

   sden = SIZE(den_bound);

   if (_ntl_gsign(nin) <= 0)
      ghalt("rational reconstruction: bad modulus");

   if (_ntl_gsign(ain) < 0 || _ntl_gcompare(ain, nin) >= 0)
      ghalt("rational reconstruction: bad residue");

      
   e = SIZE(nin);

   _ntl_gsetlength(&a, e);
   _ntl_gsetlength(&n, e);
   _ntl_gsetlength(&q, e);
   _ntl_gsetlength(&w, e);
   _ntl_gsetlength(&x, e);
   _ntl_gsetlength(&y, e);
   _ntl_gsetlength(&z, e);
   _ntl_gsetlength(&inv, e);
   _ntl_gsetlength(&u, e);
   _ntl_gsetlength(&a_bak, e);
   _ntl_gsetlength(&n_bak, e);
   _ntl_gsetlength(&inv_bak, e);
   _ntl_gsetlength(&w_bak, e);

   fhi1 = 1.0 + ((double) 32.0)/NTL_FDOUBLE_PRECISION;
   flo1 = 1.0 - ((double) 32.0)/NTL_FDOUBLE_PRECISION;

   fhi = 1.0 + ((double) 8.0)/NTL_FDOUBLE_PRECISION;
   flo = 1.0 - ((double) 8.0)/NTL_FDOUBLE_PRECISION;

   _ntl_gcopy(ain, &a);
   _ntl_gcopy(nin, &n);

   _ntl_gone(&inv);
   _ntl_gzero(&w);

   while (1)
   {
      if (SIZE(w) >= sden && _ntl_gcompare(w, den_bound) > 0) break;
      if (SIZE(n) <= snum && _ntl_gcompare(n, num_bound) <= 0) break;

      _ntl_gcopy(a, &a_bak);
      _ntl_gcopy(n, &n_bak);
      _ntl_gcopy(w, &w_bak);
      _ntl_gcopy(inv, &inv_bak);

      gotthem = 0;
      sa = SIZE(a);
      sn = SIZE(n);
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = SIZE(a);
         p = DATA(a) + (sa-1);
         num = (double) (*p) * NTL_ZZ_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_ZZ_FRADIX;
         if (sa > 2)
            num += (*(p - 1));

         sn = SIZE(n);
         p = DATA(n) + (sn-1);
         den = (double) (*p) * NTL_ZZ_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_ZZ_FRADIX;
         if (sn > 2)
            den += (*(p - 1));

         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_ZZ_FRADIX;
            lo *= NTL_ZZ_FRADIX;
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
            if (hi >= NTL_SP_BOUND)
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
                     hi = NTL_SP_BOUND;
                  temp = try11;
                  try11 = try21;
                  if ((NTL_WSP_BOUND - temp) / ilo < try21)
                     fast = 0;
                  else
                     try21 = temp + ilo * try21;
                  temp = try12;
                  try12 = try22;
                  if ((NTL_WSP_BOUND - temp) / ilo < try22)
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
         _ntl_gsmul(inv, got11, &x);
         _ntl_gsmul(w, got12, &y);
         _ntl_gsmul(inv, got21, &z);
         _ntl_gsmul(w, got22, &w);
         _ntl_gadd(x, y, &inv);
         _ntl_gadd(z, w, &w);
         _ntl_gsmul(a, got11, &x);
         _ntl_gsmul(n, got12, &y);
         _ntl_gsmul(a, got21, &z);
         _ntl_gsmul(n, got22, &n);
         _ntl_gsub(x, y, &a);
         _ntl_gsub(n, z, &n);
      }
      else
      {
         _ntl_gdiv(a, n, &q, &a);
         _ntl_gmul(q, w, &x);
         _ntl_gadd(inv, x, &inv);
         if (!ZEROP(a))
         {
            _ntl_gdiv(n, a, &q, &n);
            _ntl_gmul(q, inv, &x);
            _ntl_gadd(w, x, &w);
         }
         else
         {
            break;
         }
      }
   }

   _ntl_gcopy(a_bak, &a);
   _ntl_gcopy(n_bak, &n);
   _ntl_gcopy(w_bak, &w);
   _ntl_gcopy(inv_bak, &inv);

   _ntl_gnegate(&w);

   while (1)
   {
      sa = SIZE(w);
      if (sa < 0) SIZE(w) = -sa;
      if (SIZE(w) >= sden && _ntl_gcompare(w, den_bound) > 0) return 0;
      SIZE(w) = sa;

      if (SIZE(n) <= snum && _ntl_gcompare(n, num_bound) <= 0) break;
      
      fast = 0;
      sa = SIZE(a);
      sn = SIZE(n);
      diff = sa - sn;
      if (!diff || diff == 1)
      {
         sa = SIZE(a);
         p = DATA(a) + (sa-1);
         num = (double) (*p) * NTL_ZZ_FRADIX;
         if (sa > 1)
            num += (*(--p));
         num *= NTL_ZZ_FRADIX;
         if (sa > 2)
            num += (*(p - 1));

         sn = SIZE(n);
         p = DATA(n) + (sn-1);
         den = (double) (*p) * NTL_ZZ_FRADIX;
         if (sn > 1)
            den += (*(--p));
         den *= NTL_ZZ_FRADIX;
         if (sn > 2)
            den += (*(p - 1));

         hi = fhi1 * (num + 1.0) / den;
         lo = flo1 * num / (den + 1.0);
         if (diff > 0)
         {
            hi *= NTL_ZZ_FRADIX;
            lo *= NTL_ZZ_FRADIX;
         }

         if (hi < NTL_SP_BOUND)
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
               _ntl_gsub(inv, w, &inv);
               _ntl_gsubpos(a, n, &a);
            }
            else {
               _ntl_gsmul(w, ilo, &x);
               _ntl_gsub(inv, x, &inv);
               _ntl_gsmul(n, ilo, &x);
               _ntl_gsubpos(a, x, &a);
            }
         }
      }
      else {
         _ntl_gdiv(a, n, &q, &a);
         _ntl_gmul(q, w, &x);
         _ntl_gsub(inv, x, &inv);
      }

      _ntl_gswap(&a, &n);
      _ntl_gswap(&inv, &w);
   }

   if (_ntl_gsign(w) < 0) {
      _ntl_gnegate(&w);
      _ntl_gnegate(&n);
   }

   _ntl_gcopy(n, num_out);
   _ntl_gcopy(w, den_out);

   return 1;
}


void
_ntl_gexp(
	_ntl_gbigint a,
	long e,
	_ntl_gbigint *bb
	)
{
	long k;
	long len_a;
	static _ntl_gbigint res = 0;

	if (!e)
	{
		_ntl_gone(bb);
		return;
	}

	if (e < 0)
		ghalt("negative exponent in _ntl_gexp");

	if (_ntl_giszero(a))
	{
		_ntl_gzero(bb);
		return;
	}

	len_a = _ntl_g2log(a);
	if (len_a > (NTL_MAX_LONG-(NTL_ZZ_NBITS-1))/e)
		ghalt("overflow in _ntl_gexp");

	_ntl_gsetlength(&res, (len_a*e+NTL_ZZ_NBITS-1)/NTL_ZZ_NBITS);

	_ntl_gcopy(a, &res);
	k = 1;
	while ((k << 1) <= e)
		k <<= 1;
	while (k >>= 1) {
		_ntl_gsq(res, &res);
		if (e & k)
			_ntl_gmul(a, res, &res);
	}

	_ntl_gcopy(res, bb);
}

void
_ntl_gexps(
	long a,
	long e,
	_ntl_gbigint *bb
	)
{
	long k;
	long len_a;
	static _ntl_gbigint res = 0;

	if (!e)
	{
		_ntl_gone(bb);
		return;
	}

	if (e < 0)
		ghalt("negative exponent in _ntl_zexps");

	if (!a)
	{
		_ntl_gzero(bb);
		return;
	}

	len_a = _ntl_g2logs(a);
	if (len_a > (NTL_MAX_LONG-(NTL_ZZ_NBITS-1))/e)
		ghalt("overflow in _ntl_gexps");

	_ntl_gsetlength(&res, (len_a*e+NTL_ZZ_NBITS-1)/NTL_ZZ_NBITS);

	_ntl_gintoz(a, &res);
	k = 1;
	while ((k << 1) <= e)
		k <<= 1;
	while (k >>= 1) {
		_ntl_gsq(res, &res);
		if (e & k)
			_ntl_gsmul(res, a, &res);
	}

	_ntl_gcopy(res, bb);
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



/* DIRT: will not work with non-empty "nails" */

static
mp_limb_t neg_inv_mod_limb(mp_limb_t m0)
{
   mp_limb_t x; 
   long k;

   x = 1; 
   k = 1; 
   while (k < NTL_ZZ_NBITS) {
      x += x * (1 - x * m0);
      k <<= 1;
   }


   return - x;
}


/* Montgomery reduction:
 * This computes res = T/b^m mod N, where b = 2^{NTL_ZZ_NBITS}.
 * It is assumed that N has n limbs, and that T has at most n + m limbs.
 * Also, inv should be set to -N^{-1} mod b.
 * Finally, it is assumed that T has space allocated for n + m limbs,
 * and that res has space allocated for n limbs.  
 * Note: res should not overlap any inputs, and T is destroyed.
 * Note: res will have at most n limbs, but may not be fully reduced
 * mod N.  In general, we will have res < T/b^m + N.
 */

/* DIRT: this routine may not work with non-empty "nails" */

static
void redc(_ntl_gbigint T, _ntl_gbigint N, long m, mp_limb_t inv, 
          _ntl_gbigint res) 
{
   long n, sT, i;
   mp_limb_t *Ndata, *Tdata, *resdata, q, d, t, c;

   n = SIZE(N);
   Ndata = DATA(N);
   sT = SIZE(T);
   Tdata = DATA(T);
   resdata = DATA(res);

   for (i = sT; i < m+n; i++)
      Tdata[i] = 0;

   c = 0;
   for (i = 0; i < m; i++) {
      q = Tdata[i]*inv;
      d = mpn_addmul_1(Tdata+i, Ndata, n, q);
      t = Tdata[i+n] + d;
      Tdata[i+n] = t + c;
      if (t < d || (c == 1 && t + c  == 0)) 
         c = 1;
      else
         c = 0;
   }

   if (c) {
      mpn_sub_n(resdata, Tdata + m, Ndata, n);
   }
   else {
      for (i = 0; i < n; i++)
         resdata[i] = Tdata[m + i];
   }

   i = n;
   STRIP(i, resdata);

   SIZE(res) = i;
   SIZE(T) = 0;
}



#define REDC_CROSS (32)

void _ntl_gpowermod(_ntl_gbigint g, _ntl_gbigint e, _ntl_gbigint F,
                    _ntl_gbigint *h)

/* h = g^e mod f using "sliding window" algorithm

   remark: the notation (h, g, e, F) is strange, because I
   copied the code from BB.c.
*/

{
   _ntl_gbigint res, gg, *v, t;
   long n, i, k, val, cnt, m;
   long use_redc, sF;
   mp_limb_t inv;

   if (_ntl_gsign(g) < 0 || _ntl_gcompare(g, F) >= 0 || 
       _ntl_gscompare(F, 1) <= 0) 
      ghalt("PowerMod: bad args");

   if (_ntl_gscompare(e, 0) == 0) {
      _ntl_gone(h);
      return;
   }

   if (_ntl_gscompare(e, 1) == 0) {
      _ntl_gcopy(g, h);
      return;
   }

   if (_ntl_gscompare(e, -1) == 0) {
      _ntl_ginvmod(g, F, h);
      return;
   }

   if (_ntl_gscompare(e, 2) == 0) {
      _ntl_gsqmod(g, F, h);
      return;
   }

   if (_ntl_gscompare(e, -2) == 0) {
      res = 0;
      _ntl_gsqmod(g, F, &res);
      _ntl_ginvmod(res, F, h);
      _ntl_gfree(&res);
      return;
   }

   n = _ntl_g2log(e);

   sF = SIZE(F);

   res = 0;
   _ntl_gsetlength(&res, sF*2);

   t = 0;
   _ntl_gsetlength(&t, sF*2);

   use_redc = (DATA(F)[0] & 1) && sF < REDC_CROSS;

   gg = 0;

   if (use_redc) {
      _ntl_glshift(g, sF*NTL_ZZ_NBITS, &res);
      _ntl_gmod(res, F, &gg);

      inv = neg_inv_mod_limb(DATA(F)[0]);
   }
   else
      _ntl_gcopy(g, &gg);


   if (_ntl_gscompare(g, 2) == 0) {
      /* plain square-and-multiply algorithm, optimized for g == 2 */

      _ntl_gbigint F1 = 0;

      if (use_redc) {
         long shamt;

         COUNT_BITS(shamt, DATA(F)[sF-1]);
         shamt = NTL_ZZ_NBITS - shamt;
         _ntl_glshift(F, shamt, &F1);
      }

      _ntl_gcopy(gg, &res);

      for (i = n - 2; i >= 0; i--) {
         _ntl_gsq(res, &t);
         if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);

         if (_ntl_gbit(e, i)) {
            _ntl_gadd(res, res, &res);

            if (use_redc) {
               while (SIZE(res) > sF) {
                  _ntl_gsubpos(res, F1, &res);
               }
            }
            else {
               if (_ntl_gcompare(res, F) >= 0)
                  _ntl_gsubpos(res, F, &res);
            }
         }
      }


      if (use_redc) {
         _ntl_gcopy(res, &t);
         redc(t, F, sF, inv, res);
         if (_ntl_gcompare(res, F) >= 0) {
            _ntl_gsub(res, F, &res);
         }
      }

      if (_ntl_gsign(e) < 0) _ntl_ginvmod(res, F, &res);

      _ntl_gcopy(res, h);
      _ntl_gfree(&res);
      _ntl_gfree(&gg);
      _ntl_gfree(&t);
      _ntl_gfree(&F1);
      return;
   }


   if (n < 16) { 
      /* plain square-and-multiply algorithm */

      _ntl_gcopy(gg, &res);

      for (i = n - 2; i >= 0; i--) {
         _ntl_gsq(res, &t);
         if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);

         if (_ntl_gbit(e, i)) {
            _ntl_gmul(res, gg, &t);
            if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);
         }
      }


      if (use_redc) {
         _ntl_gcopy(res, &t);
         redc(t, F, sF, inv, res);
         if (_ntl_gcompare(res, F) >= 0) {
            _ntl_gsub(res, F, &res);
         }
      }

      if (_ntl_gsign(e) < 0) _ntl_ginvmod(res, F, &res);

      _ntl_gcopy(res, h);
      _ntl_gfree(&res);
      _ntl_gfree(&gg);
      _ntl_gfree(&t);
      return;
   }

   k = OptWinSize(n);

   if (k > 5) k = 5;

   v = (_ntl_gbigint *) NTL_MALLOC((1L << (k-1)), sizeof(_ntl_gbigint), 0);
   if (!v) ghalt("out of memory");
   for (i = 0; i < (1L << (k-1)); i++) {
      v[i] = 0; 
      _ntl_gsetlength(&v[i], sF);
   }

   _ntl_gcopy(gg, &v[0]);
 
   if (k > 1) {
      _ntl_gsq(gg, &t);
      if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);

      for (i = 1; i < (1L << (k-1)); i++) {
         _ntl_gmul(v[i-1], res, &t);
         if (use_redc) redc(t, F, sF, inv, v[i]); else _ntl_gmod(t, F, &v[i]);
      }
   }

   _ntl_gcopy(gg, &res);

   val = 0;
   for (i = n-2; i >= 0; i--) {
      val = (val << 1) | _ntl_gbit(e, i); 
      if (val == 0) {
         _ntl_gsq(res, &t);
         if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);
      }
      else if (val >= (1L << (k-1)) || i == 0) {
         cnt = 0;
         while ((val & 1) == 0) {
            val = val >> 1;
            cnt++;
         }

         m = val;
         while (m > 0) {
            _ntl_gsq(res, &t);
            if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);
            m = m >> 1;
         }

         _ntl_gmul(res, v[val >> 1], &t);
         if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);

         while (cnt > 0) {
            _ntl_gsq(res, &t);
            if (use_redc) redc(t, F, sF, inv, res); else _ntl_gmod(t, F, &res);
            cnt--;
         }

         val = 0;
      }
   }

   if (use_redc) {
      _ntl_gcopy(res, &t);
      redc(t, F, sF, inv, res);
      if (_ntl_gcompare(res, F) >= 0) {
         _ntl_gsub(res, F, &res);
      }
   }

   if (_ntl_gsign(e) < 0) _ntl_ginvmod(res, F, &res);

   _ntl_gcopy(res, h);

   _ntl_gfree(&res);
   _ntl_gfree(&gg);
   _ntl_gfree(&t);
   for (i = 0; i < (1L << (k-1)); i++)
      _ntl_gfree(&v[i]);
   free(v);
}

long _ntl_gsize(_ntl_gbigint rep)
{
   if (!rep)
      return 0;
   else if (SIZE(rep) < 0)
      return -SIZE(rep);
   else
      return SIZE(rep);
}

long _ntl_gisone(_ntl_gbigint rep)
{
   return rep != 0 && SIZE(rep) == 1 && DATA(rep)[0] == 1;
}

long _ntl_gsptest(_ntl_gbigint rep)
{
   return !rep || SIZE(rep) == 0 ||
          ((SIZE(rep) == 1 || SIZE(rep) == -1) && 
           DATA(rep)[0] < ((mp_limb_t) NTL_SP_BOUND));
}

long _ntl_gwsptest(_ntl_gbigint rep)
{
   return !rep || SIZE(rep) == 0 ||
          ((SIZE(rep) == 1 || SIZE(rep) == -1) && 
           DATA(rep)[0] < ((mp_limb_t) NTL_WSP_BOUND));
}



long _ntl_gcrtinrange(_ntl_gbigint g, _ntl_gbigint a)
{
   long sa, sg, i; 
   mp_limb_t carry, u, v;
   mp_limb_t *adata, *gdata;

   if (!a || SIZE(a) <= 0) return 0;

   sa = SIZE(a);

   if (!g) return 1;

   sg = SIZE(g);

   if (sg == 0) return 1;

   if (sg < 0) sg = -sg;

   if (sa-sg > 1) return 1;

   if (sa-sg < 0) return 0;

   adata = DATA(a);
   gdata = DATA(g);

   carry=0;

   if (sa-sg == 1) {
      if (adata[sa-1] > ((mp_limb_t) 1)) return 1;
      carry = 1;
   }

   i = sg-1;
   u = 0;
   v = 0;
   while (i >= 0 && u == v) {
      u = (carry << (NTL_ZZ_NBITS-1)) + (adata[i] >> 1);
      v = gdata[i];
      carry = (adata[i] & 1);
      i--;
   }

   if (u == v) {
      if (carry) return 1;
      return (SIZE(g) > 0);
   }
   else
      return (u > v);
}



/* DIRT: this routine will not work with non-empty "nails" */

void _ntl_gfrombytes(_ntl_gbigint *x, const unsigned char *p, long n)
{
   long BytesPerLimb;
   long lw, r, i, j;
   mp_limb_t *xp, t;

   if (n <= 0) {
      x = 0;
      return;
   }

   BytesPerLimb = NTL_ZZ_NBITS/8;


   lw = n/BytesPerLimb;
   r = n - lw*BytesPerLimb;

   if (r != 0) 
      lw++;
   else
      r = BytesPerLimb;

   _ntl_gsetlength(x, lw); 
   xp = DATA(*x);

   for (i = 0; i < lw-1; i++) {
      t = 0;
      for (j = 0; j < BytesPerLimb; j++) {
         t >>= 8;
         t += (((mp_limb_t)(*p)) & ((mp_limb_t) 255)) << ((BytesPerLimb-1)*8);
         p++;
      }
      xp[i] = t;
   }

   t = 0;
   for (j = 0; j < r; j++) {
      t >>= 8;
      t += (((mp_limb_t)(*p)) & ((mp_limb_t) 255)) << ((BytesPerLimb-1)*8);
      p++;
   }

   t >>= (BytesPerLimb-r)*8;
   xp[lw-1] = t;

   STRIP(lw, xp);
   SIZE(*x) = lw; 
}



/* DIRT: this routine will not work with non-empty "nails" */

void _ntl_gbytesfromz(unsigned char *p, _ntl_gbigint a, long n)
{
   long BytesPerLimb;
   long lbits, lbytes, min_bytes, min_words, r;
   long i, j;
   mp_limb_t *ap, t;

   if (n < 0) n = 0;

   BytesPerLimb = NTL_ZZ_NBITS/8;

   lbits = _ntl_g2log(a);
   lbytes = (lbits+7)/8;

   min_bytes = (lbytes < n) ? lbytes : n;

   min_words = min_bytes/BytesPerLimb;

   r = min_bytes - min_words*BytesPerLimb;
   if (r != 0)
      min_words++;
   else
      r = BytesPerLimb;

   if (a)
      ap = DATA(a);
   else
      ap = 0;


   for (i = 0; i < min_words-1; i++) {
      t = ap[i];
      for (j = 0; j < BytesPerLimb; j++) {
         *p = t & ((mp_limb_t) 255);
         t >>= 8;
         p++;
      }
   }

   if (min_words > 0) {
      t = ap[min_words-1];
      for (j = 0; j < r; j++) {
         *p = t & ((mp_limb_t) 255);
         t >>= 8;
         p++;
      }
   }

   for (j = min_bytes; j < n; j++) {
      *p = 0;
      p++;
   }
}




long _ntl_gblock_construct_alloc(_ntl_gbigint *x, long d, long n)
{
   long d1, sz, AllocAmt, m, j, alloc;
   char *p;
   _ntl_gbigint t;


   /* check n value */

   if (n <= 0)
      ghalt("block construct: n must be positive");



   /* check d value */

   if (d <= 0)
      ghalt("block construct: d must be positive");

   if (NTL_OVERFLOW(d, NTL_ZZ_NBITS, NTL_ZZ_NBITS))
      ghalt("block construct: d too large");

#ifdef NTL_SMALL_MP_SIZE_T
   /* this makes sure that numbers don't get too big for GMP */
   if (d >= (1L << (NTL_BITS_PER_INT-4)))
      ghalt("size too big for GMP");
#endif

   d1 = d + 1;

   if (STORAGE_OVF(d1))
      ghalt("block construct: d too large");



   sz = STORAGE(d1);

   AllocAmt = NTL_MAX_ALLOC_BLOCK/sz;
   if (AllocAmt == 0) AllocAmt = 1;

   if (AllocAmt < n)
      m = AllocAmt;
   else
      m = n;

   p = (char *) NTL_MALLOC(m, sz, 0);
   if (!p) ghalt("out of memory in _ntl_gblock_construct");

   *x = (_ntl_gbigint) p;

   for (j = 0; j < m; j++) {
      t = (_ntl_gbigint) p;
      alloc = (d1 << 2) | 1;
      if (j < m-1) alloc |= 2;
      ALLOC(t) = alloc;
      SIZE(t) = 0;
      p += sz;
   }

   return m;
}


void _ntl_gblock_construct_set(_ntl_gbigint x, _ntl_gbigint *y, long i)
{
   long d1, sz;


   d1 = ALLOC(x) >> 2;
   sz = STORAGE(d1);

   *y = (_ntl_gbigint) (((char *) x) + i*sz);
}


long _ntl_gblock_destroy(_ntl_gbigint x)
{
   long d1, sz, alloc, m;
   char *p;
   _ntl_gbigint t;

   
   d1 = ALLOC(x) >> 2;
   sz = STORAGE(d1);

   p = (char *) x;

   m = 1;

   for (;;) {
      t = (_ntl_gbigint) p;
      alloc = ALLOC(t);
      if ((alloc & 1) == 0) 
         ghalt("corrupted memory detected in _ntl_gblock_destroy");
      if ((alloc & 2) == 0) break;
      m++;
      p += sz;
   }

   free(x);
   return m;
}


long _ntl_gblock_storage(long d)
{
   long d1, sz; 

   d1 = d + 1;
   sz = STORAGE(d1) + sizeof(_ntl_gbigint);

   return sz;
}


/*
 * This is a completely portable MulMod routine.
 */

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




#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))

#define SP_MUL_MOD2(res, a, b, n, bninv) \
do {  \
   long _a = (a);  \
   long _b = (b);  \
   long _n = (n);  \
   double _bninv = (bninv);  \
   long _q, _res;  \
  \
   _q  = (long) (((double) _a) * _bninv);  \
   _res = _a*_b - _q*_n;  \
  \
   _res += (_res >> (NTL_BITS_PER_LONG-1)) & _n;  \
   _res -= _n;  \
   _res += (_res >> (NTL_BITS_PER_LONG-1)) & _n;  \
  \
   res = _res;  \
} while (0)  

#else

/*
 * This is a completely portable MulMod routine.
 */

#define SP_MUL_MOD2(res, a, b, n, bninv) \
do { \
   long _a = (a); \
   long _b = (b); \
   long _n = (n); \
   double _bninv = (bninv); \
   long _q;  \
   unsigned long _res; \
 \
   _q  = (long) (((double) _a) * _bninv); \
   _res = ((unsigned long) _a)*((unsigned long) _b)  - \
          ((unsigned long) _q)*((unsigned long) _n); \
 \
   if (_res >> (NTL_BITS_PER_LONG-1))  \
      _res += _n;  \
   else if (((long) _res) >= _n)  \
      _res -= _n;  \
 \
   res = (long) _res; \
} while (0)

#endif


static
long SpecialPower(long e, long p)
{
   long a;
   long x, y;

   a = (long) ((((mp_limb_t) 1) << (NTL_ZZ_NBITS-2)) % ((mp_limb_t) p));
   SP_MUL_MOD(a, a, 2, p);
   SP_MUL_MOD(a, a, 2, p);

   x = 1;
   y = a;
   while (e) {
      if (e & 1) SP_MUL_MOD(x, x, y, p);
      SP_MUL_MOD(y, y, y, p);
      e = e >> 1;
   }

   return x;
}


static
void sp_ext_eucl(long *dd, long *ss, long *tt, long a, long b)
{
   long  u, v, u0, v0, u1, v1, u2, v2, q, r;

   long aneg = 0, bneg = 0;

   if (a < 0) {
      if (a < -NTL_MAX_LONG) ghalt("integer overflow");
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      if (b < -NTL_MAX_LONG) ghalt("integer overflow");
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
   if (d != 1) ghalt("inverse undefined");
   if (s < 0)
      return s + n;
   else
      return s;
}

/* ------ HERE ------ */



struct crt_body_gmp {
   _ntl_gbigint *v;
   long sbuf;
   long n;
   _ntl_gbigint buf;
};

struct crt_body_gmp1 {
   long n;
   long levels;
   long *primes;
   long *inv_vec;
   long *val_vec;
   long *index_vec;
   _ntl_gbigint *prod_vec;
   _ntl_gbigint *rem_vec;
   _ntl_gbigint *coeff_vec;
   _ntl_gbigint temps[2];
   _ntl_gbigint modulus;
};


struct crt_body {
   long strategy;

   union {
      struct crt_body_gmp G;
      struct crt_body_gmp1 G1;
   } U;
};




void _ntl_gcrt_struct_init(void **crt_struct, long n, _ntl_gbigint p,
                          const long *primes)
{
   struct crt_body *c;

   c = (struct crt_body *) NTL_MALLOC(1, sizeof(struct crt_body), 0);
   if (!c) ghalt("out of memory");

   if (n >= 600) { 
      struct crt_body_gmp1 *C = &c->U.G1;
      long *q;
      long i, j;
      long levels, vec_len;
      long *val_vec, *inv_vec;
      long *index_vec;
      _ntl_gbigint *prod_vec, *rem_vec, *coeff_vec;
      _ntl_gbigint *temps;

      C->modulus = 0;
      _ntl_gcopy(p, &C->modulus);

      temps = &C->temps[0];

      temps[0] = 0;
      temps[1] = 0;
   
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) ghalt("out of memory");

      val_vec = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!val_vec) ghalt("out of memory");

      inv_vec = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!inv_vec) ghalt("out of memory");

      for (i = 0; i < n; i++)
         q[i] = primes[i];

      levels = 0;
      while ((n >> levels) >= 16) levels++;

      vec_len = (1L << levels) - 1;

      index_vec = (long *) NTL_MALLOC((vec_len+1), sizeof(long), 0);
      if (!index_vec) ghalt("out of memory");

      prod_vec = (_ntl_gbigint *) NTL_MALLOC(vec_len, sizeof(_ntl_gbigint), 0);
      if (!prod_vec) ghalt("out of memory");

      rem_vec = (_ntl_gbigint *) NTL_MALLOC(vec_len, sizeof(_ntl_gbigint), 0);
      if (!rem_vec) ghalt("out of memory");

      coeff_vec = (_ntl_gbigint *) NTL_MALLOC(n, sizeof(_ntl_gbigint), 0);
      if (!coeff_vec) ghalt("out of memory");

      for (i = 0; i < vec_len; i++)
         prod_vec[i] = 0;

      for (i = 0; i < vec_len; i++)
         rem_vec[i] = 0;

      for (i = 0; i < n; i++)
         coeff_vec[i] = 0;

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

         _ntl_gone(&prod_vec[i]);
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            _ntl_gsmul(prod_vec[i], q[j], &prod_vec[i]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            _ntl_gsdiv(prod_vec[i], q[j], &coeff_vec[j]);
      }

      for (i = (1L << (levels-1)) - 2; i >= 0; i--)
         _ntl_gmul(prod_vec[2*i+1], prod_vec[2*i+2], &prod_vec[i]);

     /*** new asymptotically fast code to compute inv_vec ***/

      _ntl_gone(&rem_vec[0]);
      for (i = 0; i < (1L << (levels-1)) - 1; i++) {
         _ntl_gmod(rem_vec[i], prod_vec[2*i+1], &temps[0]);
         _ntl_gmul(temps[0], prod_vec[2*i+2], &temps[1]);
         _ntl_gmod(temps[1], prod_vec[2*i+1], &rem_vec[2*i+1]);

         _ntl_gmod(rem_vec[i], prod_vec[2*i+2], &temps[0]);
         _ntl_gmul(temps[0], prod_vec[2*i+1], &temps[1]);
         _ntl_gmod(temps[1], prod_vec[2*i+2], &rem_vec[2*i+2]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         for (j = index_vec[i]; j < index_vec[i+1]; j++) {
            long tt, tt1, tt2;
            _ntl_gsdiv(prod_vec[i], q[j], &temps[0]);
            tt = _ntl_gsmod(temps[0], q[j]);
            tt1 = _ntl_gsmod(rem_vec[i], q[j]);
            SP_MUL_MOD(tt2, tt, tt1, q[j]);
            inv_vec[j] = sp_inv_mod(tt2, q[j]);
         }
      }



#if 0
      /* the following is asymptotically the bottleneck...but it
       * it probably doesn't matter. */

      fprintf(stderr, "checking in lip\n");
      for (i = 0; i < n; i++) {
         long tt;
         _ntl_gsdiv(prod_vec[0], q[i], &temps[0]);
         tt = mpn_mod_1(DATA(temps[0]), SIZE(temps[0]), q[i]);
         if (inv_vec[i] != sp_inv_mod(tt, q[i])) 
            fprintf(stderr, "oops in lip\n");

         /* inv_vec[i] = sp_inv_mod(tt, q[i]); */

      }
#endif

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

   {
      struct crt_body_gmp *C = &c->U.G;
      long i;
      c->strategy = 1;

      C->n = n;
      C->v = (_ntl_gbigint *) NTL_MALLOC(n, sizeof(_ntl_gbigint), 0);
      if (!C->v) ghalt("out of memory");

      for (i = 0; i < n; i++)
         C->v[i] = 0;

      C->sbuf = SIZE(p)+2;

      C->buf = 0;
      _ntl_gsetlength(&C->buf, C->sbuf);

      *crt_struct = (void *) c;
      return;
   }
}

void _ntl_gcrt_struct_insert(void *crt_struct, long i, _ntl_gbigint m)
{
   struct crt_body *c = (struct crt_body *) crt_struct;

   switch (c->strategy) {
   case 1: {
      _ntl_gcopy(m, &c->U.G.v[i]);
      break;
   }

   default:
      ghalt("_ntl_gcrt_struct_insert: inconsistent strategy");

   } /* end switch */
}


void _ntl_gcrt_struct_free(void *crt_struct)
{
   struct crt_body *c = (struct crt_body *) crt_struct;

   switch (c->strategy) {
   case 1: {
      struct crt_body_gmp *C = &c->U.G;
      long i, n;

      n = C->n;

      for (i = 0; i < n; i++)
         _ntl_gfree(&C->v[i]);

      _ntl_gfree(&C->buf);

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
      _ntl_gbigint *prod_vec = C->prod_vec;
      _ntl_gbigint *rem_vec = C->rem_vec;
      _ntl_gbigint *coeff_vec = C->coeff_vec;
      _ntl_gbigint *temps = C->temps;
      _ntl_gbigint modulus = C->modulus;
      long vec_len = (1L << levels) - 1;

      long i;

      for (i = 0; i < vec_len; i++)
         _ntl_gfree(&prod_vec[i]);

      for (i = 0; i < vec_len; i++)
         _ntl_gfree(&rem_vec[i]);

      for (i = 0; i < n; i++)
         _ntl_gfree(&coeff_vec[i]);

      _ntl_gfree(&temps[0]);
      _ntl_gfree(&temps[1]);

      _ntl_gfree(&modulus);

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

   default:

      ghalt("_ntl_gcrt_struct_free: inconsistent strategy");

   } /* end case */
}

static
void gadd_mul_many(_ntl_gbigint *res, _ntl_gbigint *a, long *b, 
                      long n, long sz)

{
   mp_limb_t *xx, *yy; 
   long i, sx;
   long sy;
   mp_limb_t carry;

   sx = sz + 2;
   if (MustAlloc(*res, sx))
      _ntl_gsetlength(res, sx);

   xx = DATA(*res);

   for (i = 0; i < sx; i++)
      xx[i] = 0;

   for (i = 0; i < n; i++) {
      if (!a[i]) continue;

      yy = DATA(a[i]);
      sy = SIZE(a[i]); 

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
   SIZE(*res) = sx;
}

void _ntl_gcrt_struct_eval(void *crt_struct, _ntl_gbigint *x, const long *b)
{
   struct crt_body *c = (struct crt_body *) crt_struct;

   switch (c->strategy) {

   case 1: {
      struct crt_body_gmp *C = &c->U.G;

      mp_limb_t *xx, *yy; 
      _ntl_gbigint *a;
      long i, sx, n;
      long sy;
      mp_limb_t carry;
   
      n = C->n;
      sx = C->sbuf;
   
      xx = DATA(C->buf);

      for (i = 0; i < sx; i++)
         xx[i] = 0;
   
      a = C->v;
   
      for (i = 0; i < n; i++) {
         if (!a[i]) continue;

         yy = DATA(a[i]);
         sy = SIZE(a[i]); 
   
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
      SIZE(C->buf) = sx;
      _ntl_gcopy(C->buf, x);
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
      _ntl_gbigint *prod_vec = C->prod_vec;
      _ntl_gbigint *rem_vec = C->rem_vec;
      _ntl_gbigint *coeff_vec = C->coeff_vec;
      _ntl_gbigint *temps = C->temps;
      long vec_len = (1L << levels) - 1;

      long i, j;

      for (i = 0; i < n; i++) {
         SP_MUL_MOD(val_vec[i], b[i], inv_vec[i], primes[i]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         long j1 = index_vec[i];
         long j2 = index_vec[i+1];
         gadd_mul_many(&rem_vec[i], &coeff_vec[j1], &val_vec[j1], j2-j1, 
                          SIZE(prod_vec[i]));
      }

      for (i = (1L << (levels-1)) - 2; i >= 0; i--) {
         _ntl_gmul(prod_vec[2*i+1], rem_vec[2*i+2], &temps[0]);
         _ntl_gmul(rem_vec[2*i+1], prod_vec[2*i+2], &temps[1]);
         _ntl_gadd(temps[0], temps[1], &rem_vec[i]);
      }

      /* temps[0] = rem_vec[0] mod prod_vec[0] (least absolute residue) */
      _ntl_gmod(rem_vec[0], prod_vec[0], &temps[0]);
      _ntl_gsub(temps[0], prod_vec[0], &temps[1]);
      _ntl_gnegate(&temps[1]);
      if (_ntl_gcompare(temps[0], temps[1]) > 0) {
         _ntl_gnegate(&temps[1]);
         _ntl_gcopy(temps[1], &temps[0]);
      }

      _ntl_gmod(temps[0], C->modulus, &temps[1]);
      _ntl_gcopy(temps[1], x);

      break;
   }

   default:

      ghalt("_crt_gstruct_eval: inconsistent strategy");

   } /* end case */

}


long _ntl_gcrt_struct_special(void *crt_struct)
{
   struct crt_body *c = (struct crt_body *) crt_struct;
   return (c->strategy == 2);
}


struct rem_body_lip {
   long n;
   long *primes;
};

struct rem_body_gmp {
   long n;
   long levels;
   long *primes;
   long *index_vec;
   _ntl_gbigint *prod_vec;
   _ntl_gbigint *rem_vec;
};


struct rem_body_gmp1 {
   long n;
   long levels;
   long *primes;
   long *index_vec;
   long *len_vec;
   mp_limb_t *inv_vec;
   long *corr_vec;
   double *corraux_vec;
   _ntl_gbigint *prod_vec;
   _ntl_gbigint *rem_vec;
};


struct rem_body {
   long strategy;

   union {
      struct rem_body_lip L;
      struct rem_body_gmp G;
      struct rem_body_gmp1 G1;
   } U;
};




void _ntl_grem_struct_init(void **rem_struct, long n, _ntl_gbigint modulus,
                          const long *p)
{
   struct rem_body *r;

   r = (struct rem_body *) NTL_MALLOC(1, sizeof(struct rem_body), 0);
   if (!r) ghalt("out of memory");

   if (n >= 32 && n <= 256) {
      struct rem_body_gmp1 *R = &r->U.G1;

      long *q;
      long i, j;
      long levels, vec_len;
      long *index_vec;
      long *len_vec, *corr_vec;
      double *corraux_vec;
      mp_limb_t *inv_vec;
      _ntl_gbigint *prod_vec, *rem_vec;
   
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) ghalt("out of memory");
   
      for (i = 0; i < n; i++)
         q[i] = p[i];

      levels = 0;
      while ((n >> levels) >= 4) levels++;

      vec_len = (1L << levels) - 1;

      index_vec = (long *) NTL_MALLOC((vec_len+1), sizeof(long), 0);
      if (!index_vec) ghalt("out of memory");

      len_vec = (long *) NTL_MALLOC(vec_len, sizeof(long), 0);
      if (!len_vec) ghalt("out of memory");

      inv_vec = (mp_limb_t *) NTL_MALLOC(vec_len, sizeof(mp_limb_t), 0);
      if (!inv_vec) ghalt("out of memory");

      corr_vec = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!corr_vec) ghalt("out of memory");

      corraux_vec = (double *) NTL_MALLOC(n, sizeof(double), 0);
      if (!corraux_vec) ghalt("out of memory");

      prod_vec = (_ntl_gbigint *) NTL_MALLOC(vec_len, sizeof(_ntl_gbigint), 0);
      if (!prod_vec) ghalt("out of memory");

      rem_vec = (_ntl_gbigint *) NTL_MALLOC(vec_len, sizeof(_ntl_gbigint), 0);
      if (!rem_vec) ghalt("out of memory");

      for (i = 0; i < vec_len; i++)
         prod_vec[i] = 0;

      for (i = 0; i < vec_len; i++)
         rem_vec[i] = 0;

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

         _ntl_gone(&prod_vec[i]);
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            _ntl_gsmul(prod_vec[i], q[j], &prod_vec[i]); 
      }

      for (i = (1L << (levels-1)) - 2; i >= 3; i--)
         _ntl_gmul(prod_vec[2*i+1], prod_vec[2*i+2], &prod_vec[i]);

      
      for (i = 3; i < vec_len; i++)
         len_vec[i] = _ntl_gsize(prod_vec[i]);

      /* Set len_vec[1] = len_vec[2] = 
       *    max(_ntl_gsize(modulus), len_vec[3..6]).
       * This is a bit paranoid, but it makes the code
       * more robust. */

      j = _ntl_gsize(modulus);
      for (i = 3; i <= 6; i++)
         if (len_vec[i] > j) j = len_vec[i];

      len_vec[1] = len_vec[2] = j;

      for (i = 3; i < vec_len; i++)
         inv_vec[i] = neg_inv_mod_limb(DATA(prod_vec[i])[0]);


      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         for (j = index_vec[i]; j < index_vec[i+1]; j++) {
            corr_vec[j] = SpecialPower(len_vec[1] - len_vec[i], q[j]);
            corraux_vec[j] = ((double) corr_vec[j])/((double) q[j]);
         }
      }


      /* allocate length in advance to streamline eval code */

      _ntl_gsetlength(&rem_vec[0], len_vec[1]); /* a special temp */

      for (i = 1; i < vec_len; i++)
         _ntl_gsetlength(&rem_vec[i], len_vec[i]);




      r->strategy = 2;
      R->n = n;
      R->primes = q;
      R->levels = levels;
      R->index_vec = index_vec;
      R->len_vec = len_vec;
      R->inv_vec = inv_vec;
      R->corr_vec = corr_vec;
      R->corraux_vec = corraux_vec;
      R->prod_vec = prod_vec;
      R->rem_vec = rem_vec;

      *rem_struct = (void *) r;
   }
   else if (n >= 32) {
      struct rem_body_gmp *R = &r->U.G;

      long *q;
      long i, j;
      long levels, vec_len;
      long *index_vec;
      _ntl_gbigint *prod_vec, *rem_vec;
   
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) ghalt("out of memory");
   
      for (i = 0; i < n; i++)
         q[i] = p[i];

      levels = 0;
      while ((n >> levels) >= 4) levels++;

      vec_len = (1L << levels) - 1;

      index_vec = (long *) NTL_MALLOC((vec_len+1), sizeof(long), 0);
      if (!index_vec) ghalt("out of memory");

      prod_vec = (_ntl_gbigint *) NTL_MALLOC(vec_len, sizeof(_ntl_gbigint), 0);
      if (!prod_vec) ghalt("out of memory");

      rem_vec = (_ntl_gbigint *) NTL_MALLOC(vec_len, sizeof(_ntl_gbigint), 0);
      if (!rem_vec) ghalt("out of memory");

      for (i = 0; i < vec_len; i++)
         prod_vec[i] = 0;

      for (i = 0; i < vec_len; i++)
         rem_vec[i] = 0;

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

         _ntl_gone(&prod_vec[i]);
         for (j = index_vec[i]; j < index_vec[i+1]; j++)
            _ntl_gsmul(prod_vec[i], q[j], &prod_vec[i]); 
      }

      for (i = (1L << (levels-1)) - 2; i >= 3; i--)
         _ntl_gmul(prod_vec[2*i+1], prod_vec[2*i+2], &prod_vec[i]);


      /* allocate length in advance to streamline eval code */

      _ntl_gsetlength(&rem_vec[1], _ntl_gsize(modulus));
      _ntl_gsetlength(&rem_vec[2], _ntl_gsize(modulus));

      for (i = 1; i < (1L << (levels-1)) - 1; i++) {
         _ntl_gsetlength(&rem_vec[2*i+1], _ntl_gsize(prod_vec[2*i+1]));
         _ntl_gsetlength(&rem_vec[2*i+2], _ntl_gsize(prod_vec[2*i+2]));
      }

      r->strategy = 1;
      R->n = n;
      R->primes = q;
      R->levels = levels;
      R->index_vec = index_vec;
      R->prod_vec = prod_vec;
      R->rem_vec = rem_vec;

      *rem_struct = (void *) r;
   }
   else
   {
      struct rem_body_lip *R = &r->U.L;

      long *q;
      long i;

      r->strategy = 0;
      R->n = n;
      q = (long *) NTL_MALLOC(n, sizeof(long), 0);
      if (!q) ghalt("out of memory");
      R->primes = q;
  
      for (i = 0; i < n; i++)
         q[i] = p[i];
  
      *rem_struct = (void *) r;
   }

}



void _ntl_grem_struct_free(void *rem_struct)
{
   struct rem_body *r = (struct rem_body *) rem_struct;

   switch (r->strategy) {

   case 0: {
      free(r->U.L.primes);
      free(r);
      break;
   }

   case 1: {
      struct rem_body_gmp *R = &r->U.G;

      long levels = R->levels;
      long vec_len = (1L << levels) - 1;
      long i;

      for (i = 0; i < vec_len; i++)
         _ntl_gfree(&R->prod_vec[i]);

      for (i = 0; i < vec_len; i++)
         _ntl_gfree(&R->rem_vec[i]);

      free(R->primes);
      free(R->index_vec);
      free(R->prod_vec);
      free(R->rem_vec);
      free(r);
      break;
   }

   case 2: {
      struct rem_body_gmp1 *R = &r->U.G1;

      long levels = R->levels;
      long vec_len = (1L << levels) - 1;
      long i;

      for (i = 0; i < vec_len; i++)
         _ntl_gfree(&R->prod_vec[i]);

      for (i = 0; i < vec_len; i++)
         _ntl_gfree(&R->rem_vec[i]);

      free(R->primes);
      free(R->index_vec);
      free(R->len_vec);
      free(R->corr_vec);
      free(R->inv_vec);
      free(R->corraux_vec);
      free(R->prod_vec);
      free(R->rem_vec);
      free(r);
      break;
   }


   default:
      ghalt("_ntl_grem_struct_free: inconsistent strategy");

   } /* end switch */
}




void _ntl_grem_struct_eval(void *rem_struct, long *x, _ntl_gbigint a)
{
   struct rem_body *r = (struct rem_body *) rem_struct;

   switch (r->strategy) {

   case 0: {
      struct rem_body_lip *R = &r->U.L;
      long n = R->n;
      long *q = R->primes;

      long j;
      mp_limb_t *adata;
      long sa;

      if (!a) 
         sa = 0;
      else
         sa = SIZE(a);

      if (sa == 0) {
         for (j = 0; j < n; j++)
            x[j] = 0;

         break;
      }

      adata = DATA(a);

      for (j = 0; j < n; j++)
         x[j] = mpn_mod_1(adata, sa, q[j]);

      break;
   }

   case 1: {
      struct rem_body_gmp *R = &r->U.G;

      long n = R->n;
      long levels = R->levels;
      long *q = R->primes;
      long *index_vec = R->index_vec;
      _ntl_gbigint *prod_vec = R->prod_vec;
      _ntl_gbigint *rem_vec = R->rem_vec;
      long vec_len = (1L << levels) - 1;

      long i, j;

      if (ZEROP(a)) {
         for (j = 0; j < n; j++)
            x[j] = 0;

         break;
      }

      _ntl_gcopy(a, &rem_vec[1]);
      _ntl_gcopy(a, &rem_vec[2]);

      for (i = 1; i < (1L << (levels-1)) - 1; i++) {
         gmod_simple(rem_vec[i], prod_vec[2*i+1], &rem_vec[2*i+1]);
         gmod_simple(rem_vec[i], prod_vec[2*i+2], &rem_vec[2*i+2]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         long lo = index_vec[i];
         long hi = index_vec[i+1];
         mp_limb_t *s1p = DATA(rem_vec[i]);
         long s1size = SIZE(rem_vec[i]);
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

   case 2: {
      struct rem_body_gmp1 *R = &r->U.G1;

      long n = R->n;
      long levels = R->levels;
      long *q = R->primes;
      long *index_vec = R->index_vec;
      long *len_vec = R->len_vec;
      long *corr_vec = R->corr_vec;
      double *corraux_vec = R->corraux_vec;
      mp_limb_t *inv_vec = R->inv_vec;
      _ntl_gbigint *prod_vec = R->prod_vec;
      _ntl_gbigint *rem_vec = R->rem_vec;
      long vec_len = (1L << levels) - 1;

      long i, j;

      if (ZEROP(a)) {
         for (j = 0; j < n; j++)
            x[j] = 0;

         break;
      }

      _ntl_gcopy(a, &rem_vec[1]);
      _ntl_gcopy(a, &rem_vec[2]);

      for (i = 1; i < (1L << (levels-1)) - 1; i++) {
         _ntl_gcopy(rem_vec[i], &rem_vec[0]);
         redc(rem_vec[0], prod_vec[2*i+1], len_vec[i]-len_vec[2*i+1],
              inv_vec[2*i+1], rem_vec[2*i+1]);
         redc(rem_vec[i], prod_vec[2*i+2], len_vec[i]-len_vec[2*i+2],
              inv_vec[2*i+2], rem_vec[2*i+2]);
      }

      for (i = (1L << (levels-1)) - 1; i < vec_len; i++) {
         long lo = index_vec[i];
         long hi = index_vec[i+1];
         mp_limb_t *s1p = DATA(rem_vec[i]);
         long s1size = SIZE(rem_vec[i]);
         if (s1size == 0) {
            for (j = lo; j <hi; j++)
               x[j] = 0;
         }
         else {
            for (j = lo; j < hi; j++) {
               long t = mpn_mod_1(s1p, s1size, q[j]);
               SP_MUL_MOD2(x[j], t, corr_vec[j], q[j], corraux_vec[j]);
            }
         }
      }

      break;
   }

   default:
      ghalt("_ntl_grem_struct_eval: inconsistent strategy");


   } /* end switch */


}



/* routines for x += a*b for single-precision b 
 * DIRT: relies crucially on mp_limb_t being at least as
 * wide as a long.
 * Lightly massaged code taken from GMP's mpz routines */



#define _ntl_mpn_com_n(d,s,n)                                \
  do {                                                  \
    mp_limb_t *  __d = (d);                               \
    mp_limb_t *  __s = (s);                               \
    long  __n = (n);                               \
    do                                                  \
      *__d++ = (~ *__s++);              \
    while (--__n);                                      \
  } while (0)


#define _ntl_MPN_MUL_1C(cout, dst, src, size, n, cin)        \
  do {                                                  \
    mp_limb_t __cy;                                     \
    __cy = mpn_mul_1 (dst, src, size, n);               \
    (cout) = __cy + mpn_add_1 (dst, dst, size, cin);    \
  } while (0)

#define _ntl_g_inc(p, n)   \
  do {   \
    mp_limb_t * __p = (p);  \
    long __n = (n);  \
    while (__n > 0) {  \
       (*__p)++;  \
       if (*__p != 0) break;  \
       __p++;  \
       __n--;  \
    }  \
  } while (0);

#define _ntl_g_inc_carry(c, p, n)   \
  do {   \
    mp_limb_t * __p = (p);  \
    long __n = (n);  \
    long __addc = 1; \
    while (__n > 0) {  \
       (*__p)++;  \
       if (*__p != 0) { __addc = 0; break; }  \
       __p++;  \
       __n--;  \
    }  \
    c += __addc; \
  } while (0);

#define _ntl_g_dec(p, n)   \
  do {   \
    mp_limb_t * __p = (p);  \
    mp_limb_t __tmp; \
    long __n = (n);  \
    while (__n > 0) {  \
       __tmp = *__p; \
       (*__p)--;  \
       if (__tmp != 0) break;  \
       __p++;  \
       __n--;  \
    }  \
  } while (0);
  


/* sub==0 means an addmul w += x*y, sub==1 means a submul w -= x*y. */
void
_ntl_gaorsmul_1(_ntl_gbigint x, long yy, long sub, _ntl_gbigint *ww)
{
  long  xsize, wsize, wsize_signed, new_wsize, min_size, dsize;
  _ntl_gbigint w;
  mp_limb_t *xp;
  mp_limb_t *wp;
  mp_limb_t  cy;
  mp_limb_t  y;

  if (ZEROP(x) || yy == 0)
    return;

  if (ZEROP(*ww)) {
    _ntl_gsmul(x, yy, ww);
    if (sub) SIZE(*ww) = -SIZE(*ww);
    return;
  }

  if (yy == 1) {
    if (sub)
      _ntl_gsub(*ww, x, ww);
    else
      _ntl_gadd(*ww, x, ww);
    return;
  }

  if (yy == -1) {
    if (sub)
      _ntl_gadd(*ww, x, ww);
    else
      _ntl_gsub(*ww, x, ww);
    return;
  }

  if (*ww == x) {
    static _ntl_gbigint tmp = 0;
    _ntl_gsmul(x, yy, &tmp);
    if (sub)
       _ntl_gsub(*ww, tmp, ww);
    else
       _ntl_gadd(*ww, tmp, ww);
    return;
  }

  xsize = SIZE(x);
  if (xsize < 0) {
    xsize = -xsize;
    sub = 1-sub;
  }

  if (yy < 0) {
    y = - ((mp_limb_t) yy); /* careful! */
    sub = 1-sub;
  }
  else {
    y = (mp_limb_t) yy;
  }
    

  w = *ww;

  wsize_signed = SIZE(w);
  if (wsize_signed < 0) {
    sub = 1-sub;
    wsize = -wsize_signed;
  }
  else {
    wsize = wsize_signed;
  }


  if (wsize > xsize) {
    new_wsize = wsize;
    min_size = xsize;
  }
  else {
    new_wsize = xsize;
    min_size = wsize;
  }

  if (MustAlloc(w, new_wsize+1)) {
    _ntl_gsetlength(&w, new_wsize+1);
    *ww = w;
  }

  wp = DATA(w);
  xp = DATA(x);

  if (sub == 0)
    {
      /* addmul of absolute values */

      cy = mpn_addmul_1 (wp, xp, min_size, y);
      wp += min_size;
      xp += min_size;

      dsize = xsize - wsize;
      if (dsize != 0)
        {
          mp_limb_t  cy2;
          if (dsize > 0) {
            cy2 = mpn_mul_1 (wp, xp, dsize, y);
          }
          else
            {
              dsize = -dsize;
              cy2 = 0;
            }
          cy = cy2 + mpn_add_1 (wp, wp, dsize, cy);
        }

      wp[dsize] = cy;
      new_wsize += (cy != 0);
    }
  else
    {
      /* submul of absolute values */

      cy = mpn_submul_1 (wp, xp, min_size, y);
      if (wsize >= xsize)
        {
          /* if w bigger than x, then propagate borrow through it */
          if (wsize != xsize) {
            cy = mpn_sub_1 (wp+xsize, wp+xsize, wsize-xsize, cy);
          }

          if (cy != 0)
            {
              /* Borrow out of w, take twos complement negative to get
                 absolute value, flip sign of w.  */
              wp[new_wsize] = ~-cy;  /* extra limb is 0-cy */
              _ntl_mpn_com_n (wp, wp, new_wsize);
              new_wsize++;
              _ntl_g_inc(wp, new_wsize);
              wsize_signed = -wsize_signed;
            }
        }
      else /* wsize < xsize */
        {
          /* x bigger than w, so want x*y-w.  Submul has given w-x*y, so
             take twos complement and use an mpn_mul_1 for the rest.  */

          mp_limb_t  cy2;

          /* -(-cy*b^n + w-x*y) = (cy-1)*b^n + ~(w-x*y) + 1 */
          _ntl_mpn_com_n (wp, wp, wsize);
          _ntl_g_inc_carry(cy, wp, wsize);
          cy -= 1;

          /* If cy-1 == -1 then hold that -1 for latter.  mpn_submul_1 never
             returns cy==MP_LIMB_T_MAX so that value always indicates a -1. */
          cy2 = (cy == ((mp_limb_t) -1));
          cy += cy2;
          _ntl_MPN_MUL_1C (cy, wp+wsize, xp+wsize, xsize-wsize, y, cy);
          wp[new_wsize] = cy;
          new_wsize += (cy != 0);

          /* Apply any -1 from above.  The value at wp+wsize is non-zero
             because y!=0 and the high limb of x will be non-zero.  */
          if (cy2) {
            _ntl_g_dec(wp+wsize, new_wsize-wsize);
          }

          wsize_signed = -wsize_signed;
        }

      /* submul can produce high zero limbs due to cancellation, both when w
         has more limbs or x has more  */
      STRIP(new_wsize, wp);
    }

  SIZE(w) = (wsize_signed >= 0 ? new_wsize : -new_wsize);
}


void
_ntl_gsaddmul(_ntl_gbigint x, long yy,  _ntl_gbigint *ww)
{
  _ntl_gaorsmul_1(x, yy, 0, ww);
}

void
_ntl_gssubmul(_ntl_gbigint x, long yy,  _ntl_gbigint *ww)
{
  _ntl_gaorsmul_1(x, yy, 1, ww);
}


void
_ntl_gaorsmul(_ntl_gbigint x, _ntl_gbigint y, long sub,  _ntl_gbigint *ww)
{
   static _ntl_gbigint tmp = 0;

   _ntl_gmul(x, y, &tmp);
   if (sub)
      _ntl_gsub(*ww, tmp, ww);
   else
      _ntl_gadd(*ww, tmp, ww);
}


void
_ntl_gaddmul(_ntl_gbigint x, _ntl_gbigint y,  _ntl_gbigint *ww)
{
  _ntl_gaorsmul(x, y, 0, ww);
}

void
_ntl_gsubmul(_ntl_gbigint x, _ntl_gbigint y,  _ntl_gbigint *ww)
{
  _ntl_gaorsmul(x, y, 1, ww);
}




