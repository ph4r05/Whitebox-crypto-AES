

typedef long * _ntl_verylong;

#if (defined(NTL_SINGLE_MUL))

#if (defined(NTL_AVOID_FLOAT) || defined(NTL_LONG_LONG))
#error "at most one of NTL_SINGLE_MUL NTL_AVOID_FLOAT NTL_LONG_LONG allowed"
#endif

#elif (defined(NTL_AVOID_FLOAT) && defined(NTL_LONG_LONG))
#error "at most one of NTL_SINGLE_MUL NTL_AVOID_FLOAT NTL_LONG_LONG allowed"
#endif

#if (defined(NTL_SINGLE_MUL))

#if (!NTL_SINGLE_MUL_OK)
#error "NTL_SINGLE_MUL not supported on this platform"
#endif

#if (defined(NTL_CLEAN_INT))
#error "NTL_SINGLE_MUL not allowed with NTL_CLEAN_INT"
#endif


#define NTL_NBITS (26)

#else

#define NTL_NBITS NTL_NBITS_MAX

#endif


#define NTL_RADIX           (1L<<NTL_NBITS)
#define NTL_NBITSH          (NTL_NBITS>>1)
#define NTL_RADIXM          (NTL_RADIX-1)
#define NTL_RADIXROOT       (1L<<NTL_NBITSH)
#define NTL_RADIXROOTM      (NTL_RADIXROOT-1)

#define NTL_FRADIX ((double) NTL_RADIX)
#define NTL_FRADIX_INV  (((double) 1.0)/((double) NTL_RADIX))



#define NTL_ZZ_NBITS NTL_NBITS
#define NTL_ZZ_FRADIX ((double) (1L << NTL_NBITS))

#define NTL_SP_NBITS NTL_NBITS
#define NTL_SP_BOUND (1L << NTL_SP_NBITS)
#define NTL_SP_FBOUND ((double) NTL_SP_BOUND)

#define NTL_WSP_NBITS NTL_ZZ_NBITS
#define NTL_WSP_BOUND (1L << NTL_WSP_NBITS)



#if (defined(NTL_SINGLE_MUL) && !NTL_SINGLE_MUL_OK)
#undef NTL_SINGLE_MUL
#endif

#if (defined(NTL_SINGLE_MUL))


/****************************************************************

The following macros extract the two words of a double,
getting around the type system.
This is only used in the NTL_SINGLE_MUL strategy.

*****************************************************************/

#if (NTL_DOUBLES_LOW_HIGH)
#define NTL_LO_WD 0
#define NTL_HI_WD 1
#else
#define NTL_LO_WD 1
#define NTL_HI_WD 0
#endif


typedef union { double d; unsigned long rep[2]; } _ntl_d_or_rep;

#define NTL_FetchHiLo(hi,lo,x) \
do { \
   _ntl_d_or_rep ll_xx; \
   ll_xx.d = (x); \
   hi = ll_xx.rep[NTL_HI_WD]; \
   lo = ll_xx.rep[NTL_LO_WD]; \
} while (0)


#define NTL_FetchLo(lo,x)  \
do {  \
   _ntl_d_or_rep ll_xx;  \
   ll_xx.d = x;  \
   lo = ll_xx.rep[NTL_LO_WD];  \
} while (0) 

#endif


/**********************************************************************

   A multiprecision integer is represented as a pointer to long.
   Let x be such a pointer.
   x = 0 represents 0.
   Otherwise, let n = abs(x[0]) and s = sign(x[0]);
   the integer represented is then:

      s*(x[1] + x[2]*NTL_RADIX + ... + x[n]*NTL_RADIX^{n-1}),

   where x[n] != 0, unless n = s = 1.
   Notice that the number zero can be represented in precisely 2 ways,
   either as a null pointer, or as x[0] = 1 and x[1] = 0.

   Storage is generally managed via _ntl_zsetlength and _ntl_zfree.
   x[-1] = (k << 1) | b, where k is the maximum value of n allocated,
   and b is a bit that is set is the space is managed by some
   mechanism other than _ntl_zsetlength and _ntl_zfree.

************************************************************************/



#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
extern "C" {
#endif


/***********************************************************************

   Basic Functions

***********************************************************************/
    


    void _ntl_zsadd(_ntl_verylong a, long d, _ntl_verylong *b);
       /* *b = a + d */

    void _ntl_zadd(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
       /*  *c = a + b */

    void _ntl_zsub(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
       /* *c = a - b */

    void _ntl_zsubpos(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
       /* *c = a - b; assumes a >= b >= 0 */

    void _ntl_zsmul(_ntl_verylong a, long d, _ntl_verylong *b);
       /* *b = d * a */

    void _ntl_zmul(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
       /* *c = a * b */

    void _ntl_zsq(_ntl_verylong a, _ntl_verylong *c);
       /* *c = a * a */

    long _ntl_zsdiv(_ntl_verylong a, long b, _ntl_verylong *q);
       /* (*q) = floor(a/b) and a - floor(a/b)*(*q) is returned;
          error is raised if b == 0;
          if b does not divide a, then sign(*q) == sign(b) */

    void _ntl_zdiv(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *q, _ntl_verylong *r);
       /* (*q) = floor(a/b) and (*r) = a - floor(a/b)*(*q);
          error is raised if b == 0;
          if b does not divide a, then sign(*q) == sign(b) */

    void _ntl_zmultirem(_ntl_verylong a, long n, long* dd, long* rr);
    void _ntl_zmultirem2(_ntl_verylong a, long n, long* dd, double **tbl, long* rr);
       /* rr[i] = a % dd[i], i = 0..n-1;
          assumes a >= 0, 0 < dd[i] < NTL_RADIX
          _ntl_zmultirem2 takes an extra argument, tbl, which contains
          pre-computed residues of powers of RADIX */
    void _ntl_zmultirem3(_ntl_verylong a, long n, long* dd, long **tbl, long* rr);
       /* same as above, but tbl has different type */

    long _ntl_zsfastrem(_ntl_verylong a, long d);
       /* return a % d;
          assumes a >= 0, 0 < d < NTL_RADIX */
        

    void _ntl_zmod(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *r);
       /* same as _ntl_zdiv, but only remainder is computed */

    long _ntl_zsmod(_ntl_verylong a, long d);
       /* same as _ntl_zsdiv, but only remainder is computed */

    void _ntl_zquickmod(_ntl_verylong *r, _ntl_verylong b);
       /* *r = *r % b;
	  assumes b > 0 and *r >= 0;
	  The division is performed in place (but may sometimes
          cause *r to grow by one digit) */

    void _ntl_zsaddmul(_ntl_verylong x, long y,  _ntl_verylong *ww);
      /* *ww += x*y */

    void _ntl_zaddmul(_ntl_verylong x, _ntl_verylong y,  _ntl_verylong *ww);
      /* *ww += x*y */

    void _ntl_zssubmul(_ntl_verylong x, long y,  _ntl_verylong *ww);
      /* *ww -= x*y */

    void _ntl_zsubmul(_ntl_verylong x, _ntl_verylong y,  _ntl_verylong *ww);
      /* *ww -= x*y */


/********************************************************************

   Shifting and bit manipulation

*********************************************************************/

    void _ntl_z2mul(_ntl_verylong n, _ntl_verylong *a);
       /* *a = 2 * n */

    long _ntl_z2div(_ntl_verylong n, _ntl_verylong *a);
       /* *a = sign(n) * (|n|/2) */ 

    void _ntl_zlshift(_ntl_verylong n, long k, _ntl_verylong *a);
       /* *a = sign(n) * (|n| << k);
          shift is in reverse direction for negative k */

    void _ntl_zrshift(_ntl_verylong n, long k, _ntl_verylong *a);
       /* *a = sign(n) * (|n| >> k);
          shift is in reverse direction for negative k */
    
    long _ntl_zmakeodd(_ntl_verylong *n);
       /*
          if (n != 0)
              *n = m;
              return (k such that n == 2 ^ k * m with m odd);
          else
              return (0); 
        */

    long _ntl_znumtwos(_ntl_verylong n);
        /* return largest e such that 2^e divides n, or zero if n is zero */



    long _ntl_zodd(_ntl_verylong a);
       /* returns 1 if n is odd and 0 if it is even */

    long _ntl_zbit(_ntl_verylong a, long p);
       /* returns p-th bit of a, where the low order bit is indexed by 0;
          p out of range returns 0 */

    long _ntl_zsetbit(_ntl_verylong *a, long p);
       /* returns original value of p-th bit of |a|, and replaces
          p-th bit of a by 1 if it was zero;
          error if p < 0 */

    long _ntl_zswitchbit(_ntl_verylong *a, long p);
       /* returns original value of p-th bit of |a|, and switches
          the value of p-th bit of a;
          p starts counting at 0;
          error if p < 0 */


     void _ntl_zlowbits(_ntl_verylong a, long k, _ntl_verylong *b);
        /* places k low order bits of |a| in b */ 

     long _ntl_zslowbits(_ntl_verylong a, long k);
        /* returns k low order bits of |a| */

    long _ntl_zweights(long a);
        /* returns Hamming weight of |a| */

    long _ntl_zweight(_ntl_verylong a);
        /* returns Hamming weight of |a| */

    void _ntl_zand(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
        /* c gets bit pattern `bits of |a|` and `bits of |b|` */

    void _ntl_zor(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
        /* c gets bit pattern `bits of |a|` inclusive or `bits of |b|` */

    void _ntl_zxor(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
        /* c gets bit pattern `bits of |a|` exclusive or `bits of |b|` */




/************************************************************************

   Comparison

*************************************************************************/

    long _ntl_zcompare(_ntl_verylong a, _ntl_verylong b);
       /*
          if (a > b)
              return (1);
          if (a == b)
              return (0);
          if (a < b)
              return (-1);
         */

    long _ntl_zscompare(_ntl_verylong a, long b);
       /* single-precision version of the above */

    long _ntl_ziszero (_ntl_verylong a);
       /* test for 0 */


    long _ntl_zsign(_ntl_verylong a);
       /* 
          if (a > 0)
              return (1);
          if (a == 0)
              return (0);
          if (a < 0)
              return (-1);
        */

    void _ntl_zabs(_ntl_verylong *a);
       /* *a = |a| */

    void _ntl_znegate(_ntl_verylong *a);
       /* *a = -a */

    void _ntl_zcopy(_ntl_verylong a, _ntl_verylong *b);
       /* *b = a;  space is allocated  */

    void _ntl_zcopy1(_ntl_verylong a, _ntl_verylong *b);
       /* *b = a;  space not necessarily allocated  */

    void _ntl_zswap(_ntl_verylong *a, _ntl_verylong *b);
       /* swap a and b (by swaping pointers) */

    long _ntl_z2log(_ntl_verylong a);
       /* number of bits in |a|; returns 0 if a = 0 */

    long _ntl_z2logs(long a);
        /* single-precision version of the above */


/********************************************************************

   Conversion

*********************************************************************/
        
    void _ntl_zzero(_ntl_verylong *a);
       /* *a = 0;  space is allocated */

    void _ntl_zzero1(_ntl_verylong *a);
       /* *a = 0;  space not necessarily allocated */

    void _ntl_zone(_ntl_verylong *a);
       /* *a = 1 */

    void _ntl_zintoz(long d, _ntl_verylong *a);
       /* *a = d;  space is allocated  */

    void _ntl_zintoz1(long d, _ntl_verylong *a);
       /* *a = d;  space not necessarily allocated  */

    void _ntl_zuintoz(unsigned long d, _ntl_verylong *a);
       /* *a = d;  space is allocated  */

    long _ntl_ztoint(_ntl_verylong a);
       /* converts a to a long;  overflow results in value
          mod 2^{NTL_BITS_PER_LONG}. */

    unsigned long _ntl_ztouint(_ntl_verylong a);
       /* converts a to a long;  overflow results in value
          mod 2^{NTL_BITS_PER_LONG}. */


    double _ntl_zdoub(_ntl_verylong n);
       /* converts a to a double;  no overflow check */

    long _ntl_zround_correction(_ntl_verylong a, long k, long residual);
       /* k >= 1, |a| >= 2^k, and residual is 0, 1, or -1.
          The result is what we should add to (a >> k) to round
          x = a/2^k to the nearest integer using IEEE-like rounding rules
          (i.e., round to nearest, and round to even to break ties).
          The result is either 0 or sign(a).
          If residual is not zero, it is as if x were replaced by
          x' = x + residual*2^{-(k+1)}.
          This can be used to break ties when x is exactly
          half way between two integers. */

    double _ntl_zlog(_ntl_verylong a);
       /* computes log(a), protecting against overflow */
      

    void _ntl_zdoubtoz(double a, _ntl_verylong *x);
       /* x = floor(a); */  
    



/************************************************************************

   Square roots

*************************************************************************/


    long _ntl_zsqrts(long n);
       /* return floor(sqrt(n));  error raised in n < 0 */

    void _ntl_zsqrt(_ntl_verylong n, _ntl_verylong *r);
       /* *r =  floor(sqrt(n));  error raised in n < 0 */

/*********************************************************************
 
    Exponentiation
 
**********************************************************************/

   void _ntl_zexp(_ntl_verylong a, long e, _ntl_verylong *b);
       /* *b = a^e;  error raised if e < 0 */

   void _ntl_zexps(long a, long e, _ntl_verylong *b);
       /* *b = a^e;  error raised if e < 0 */
       

/*********************************************************************

   Modular Arithmetic

   Addition, subtraction, multiplication, squaring division, inversion,
   and exponentiation modulo a positive modulus n, where all operands
   (except for the exponent in exponentiation) and results are in the
   range [0, n-1].   

***********************************************************************/

    void _ntl_zaddmod(_ntl_verylong a, _ntl_verylong b, _ntl_verylong n, _ntl_verylong *c);
       /* *c = (a + b) % n */

    void _ntl_zsubmod(_ntl_verylong a, _ntl_verylong b, _ntl_verylong n, _ntl_verylong *c);
       /* *c = (a - b) % n */

    void _ntl_zsmulmod(_ntl_verylong a, long b, _ntl_verylong n, _ntl_verylong *c);
       /* *c = (a * b) % n */

    void _ntl_zmulmod(_ntl_verylong a, _ntl_verylong b, _ntl_verylong n, _ntl_verylong *c);
       /* *c = (a * b) % n */

    void _ntl_zsqmod(_ntl_verylong a, _ntl_verylong n, _ntl_verylong *c);
       /* *c = (a ^ 2) % n */

    void _ntl_zinvmod(_ntl_verylong a, _ntl_verylong n, _ntl_verylong *c);
       /* *c = (1 / a) % n; error raised if gcd(b, n) != 1 */

    void _ntl_zpowermod(_ntl_verylong g, _ntl_verylong e, _ntl_verylong F,
                        _ntl_verylong *h);

       /* *b = (a ^ e) % n;  */



/**************************************************************************

   Euclidean Algorithms

***************************************************************************/
    void _ntl_zgcd(_ntl_verylong m1, _ntl_verylong m2, _ntl_verylong *r);
       /* *r = greatest common divisor of m1 and m2; 
          uses binary gcd algorithm */


    void _ntl_zexteucl(_ntl_verylong a, _ntl_verylong *xa,
                 _ntl_verylong b, _ntl_verylong *xb,
                 _ntl_verylong *d);
       /*
          *d = a * *xa + b * *xb = gcd(a, b);
          sets *d, *xa and *xb given a and b;
          uses Lehmer`s trick
        */


    long _ntl_zinv(_ntl_verylong a, _ntl_verylong b, _ntl_verylong *c);
       /*
          if (a and b coprime)
          {
              *c = inv; 
              return(0);
          }
          else
          {
              *c = gcd(a, b);
              return(1);
          }
          
          where inv is such that (inv * a)  == 1 mod b;
          error raised if b <= 1 or a < 0 or a >= b
        */

     long _ntl_zxxratrecon(_ntl_verylong x, _ntl_verylong m,  
                      _ntl_verylong a_bound, _ntl_verylong b_bound,
                      _ntl_verylong *a, _ntl_verylong *b);

        /* rational reconstruction: see doc in ZZ.txt */


        
/**********************************************************************

    Storage Allocation

    These routines use malloc and free.

***********************************************************************/


    long _ntl_zmaxalloc(_ntl_verylong x); 
      /* max allocation request, possibly rounded up a bit */


    void _ntl_zsetlength(_ntl_verylong *v, long len);
       /* Allocates enough space to hold a len-digit number,
          where each digit has NTL_NBITS bits.
          If space must be allocated, space for one extra digit
          is always allocated. */

    void _ntl_zfree(_ntl_verylong *x);
       /* Free's space held by x, and sets x back to 0. */


/*******************************************************************

    Special routines

********************************************************************/



long _ntl_zsize(_ntl_verylong n);
long _ntl_zisone(_ntl_verylong n);
long _ntl_zdigit(_ntl_verylong a, long i);

long _ntl_zsptest(_ntl_verylong a);
long _ntl_zwsptest(_ntl_verylong a);

long _ntl_zcrtinrange(_ntl_verylong g, _ntl_verylong a);

void _ntl_zfrombytes(_ntl_verylong *x, const unsigned char *p, long n);
void _ntl_zbytesfromz(unsigned char *p, _ntl_verylong a, long nn);

long _ntl_zblock_construct_alloc(_ntl_verylong *x, long d, long n);
void _ntl_zblock_construct_set(_ntl_verylong x, _ntl_verylong *y, long i);
long _ntl_zblock_destroy(_ntl_verylong x);
long _ntl_zblock_storage(long d);



void _ntl_crt_struct_init(void **crt_struct, long n, _ntl_verylong p,
                          const long *primes);
void _ntl_crt_struct_insert(void *crt_struct, long i, _ntl_verylong m);
void _ntl_crt_struct_free(void *crt_struct);
void _ntl_crt_struct_eval(void *crt_struct, _ntl_verylong *t, const long *a);
long _ntl_crt_struct_special(void *crt_struct);

void _ntl_rem_struct_init(void **rem_struct, long n, _ntl_verylong p, 
                          const long *primes);
void _ntl_rem_struct_free(void *rem_struct);
void _ntl_rem_struct_eval(void *rem_struct, long *x, _ntl_verylong a);



#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
}
#endif


extern int _ntl_gmp_hack;

#define NTL_crt_struct_eval _ntl_crt_struct_eval
#define NTL_crt_struct_free _ntl_crt_struct_free
#define NTL_crt_struct_init _ntl_crt_struct_init
#define NTL_crt_struct_insert _ntl_crt_struct_insert
#define NTL_crt_struct_special _ntl_crt_struct_special
#define NTL_rem_struct_eval _ntl_rem_struct_eval
#define NTL_rem_struct_free _ntl_rem_struct_free
#define NTL_rem_struct_init _ntl_rem_struct_init
#define NTL_verylong _ntl_verylong
#define NTL_z2log _ntl_z2log
#define NTL_zabs _ntl_zabs
#define NTL_zadd _ntl_zadd
#define NTL_zaddmod _ntl_zaddmod
#define NTL_zand _ntl_zand
#define NTL_zbit _ntl_zbit
#define NTL_zblock_construct_alloc _ntl_zblock_construct_alloc
#define NTL_zblock_construct_set _ntl_zblock_construct_set
#define NTL_zblock_destroy _ntl_zblock_destroy
#define NTL_zblock_storage _ntl_zblock_storage
#define NTL_zbytesfromz _ntl_zbytesfromz
#define NTL_zcompare _ntl_zcompare
#define NTL_zcopy _ntl_zcopy1
#define NTL_zcrtinrange _ntl_zcrtinrange
#define NTL_zdigit _ntl_zdigit
#define NTL_zdiv _ntl_zdiv
#define NTL_zdoub _ntl_zdoub
#define NTL_zdoubtoz _ntl_zdoubtoz
#define NTL_zexp _ntl_zexp
#define NTL_zexps _ntl_zexps
#define NTL_zexteucl _ntl_zexteucl
#define NTL_zfree _ntl_zfree
#define NTL_zfrombytes _ntl_zfrombytes
#define NTL_zgcd _ntl_zgcd
#define NTL_zintoz _ntl_zintoz1
#define NTL_zinv _ntl_zinv
#define NTL_zinvmod _ntl_zinvmod
#define NTL_zisone _ntl_zisone
#define NTL_ziszero _ntl_ziszero
#define NTL_zlog _ntl_zlog
#define NTL_zlowbits _ntl_zlowbits
#define NTL_zlshift _ntl_zlshift
#define NTL_zmakeodd _ntl_zmakeodd
#define NTL_zmod _ntl_zmod
#define NTL_zmul _ntl_zmul
#define NTL_zmulmod _ntl_zmulmod
#define NTL_znegate _ntl_znegate
#define NTL_znumtwos _ntl_znumtwos
#define NTL_zodd _ntl_zodd
#define NTL_zone _ntl_zone
#define NTL_zor _ntl_zor
#define NTL_zpowermod _ntl_zpowermod
#define NTL_zquickmod _ntl_zquickmod
#define NTL_zround_correction _ntl_zround_correction
#define NTL_zrshift _ntl_zrshift
#define NTL_zsadd _ntl_zsadd
#define NTL_zscompare _ntl_zscompare
#define NTL_zsdiv _ntl_zsdiv
#define NTL_zsetbit _ntl_zsetbit
#define NTL_zmaxalloc _ntl_zmaxalloc
#define NTL_zsetlength _ntl_zsetlength
#define NTL_zsign _ntl_zsign
#define NTL_zsize _ntl_zsize
#define NTL_zslowbits _ntl_zslowbits
#define NTL_zsmod _ntl_zsmod
#define NTL_zsmul _ntl_zsmul
#define NTL_zsmulmod _ntl_zsmulmod
#define NTL_zsptest _ntl_zsptest
#define NTL_zsq _ntl_zsq
#define NTL_zsqmod _ntl_zsqmod
#define NTL_zsqrt _ntl_zsqrt
#define NTL_zsqrts _ntl_zsqrts
#define NTL_zsub _ntl_zsub
#define NTL_zsubmod _ntl_zsubmod
#define NTL_zsubpos _ntl_zsubpos
#define NTL_zswap _ntl_zswap
#define NTL_zswitchbit _ntl_zswitchbit
#define NTL_ztoint _ntl_ztoint
#define NTL_ztouint _ntl_ztouint
#define NTL_zuintoz _ntl_zuintoz
#define NTL_zweight _ntl_zweight
#define NTL_zweights _ntl_zweights
#define NTL_zwsptest _ntl_zwsptest
#define NTL_zxor _ntl_zxor
#define NTL_zxxratrecon _ntl_zxxratrecon
#define NTL_zzero _ntl_zzero1

#define NTL_zsaddmul _ntl_zsaddmul
#define NTL_zaddmul _ntl_zaddmul
#define NTL_zssubmul _ntl_zssubmul
#define NTL_zsubmul _ntl_zsubmul

