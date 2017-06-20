

#ifdef NTL_SINGLE_MUL
#error "do not set NTL_SINGLE_MUL when NTL_GMP_LIP is set"
#endif

#if 1

typedef void *_ntl_gbigint;

#else

/*
 * This way of defining the bigint handle type is a bit non-standard,
 * but better for debugging.
 */

struct _ntl_gbigint_is_opaque { int _x_; };
typedef struct _ntl_gbigint_is_opaque * _ntl_gbigint;

#endif

#define NTL_SP_NBITS NTL_NBITS_MAX
#define NTL_SP_BOUND (1L << NTL_SP_NBITS)
#define NTL_SP_FBOUND ((double) NTL_SP_BOUND)

#define NTL_WSP_NBITS (NTL_BITS_PER_LONG-2)
#define NTL_WSP_BOUND (1L << NTL_WSP_NBITS)

/* define the following so an error is raised */

#define NTL_RADIX ......
#define NTL_NBITSH ......
#define NTL_RADIXM ......
#define NTL_RADIXROOT ......
#define NTL_RADIXROOTM ......
#define NTL_FRADIX_INV ......




#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
extern "C" {
#endif


/***********************************************************************

   Basic Functions

***********************************************************************/
    


    void _ntl_gsadd(_ntl_gbigint a, long d, _ntl_gbigint *b);
       /* *b = a + d */

    void _ntl_gadd(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
       /*  *c = a + b */

    void _ntl_gsub(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
       /* *c = a - b */

    void _ntl_gsubpos(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
       /* *c = a - b; assumes a >= b >= 0 */

    void _ntl_gsmul(_ntl_gbigint a, long d, _ntl_gbigint *b);
       /* *b = d * a */

    void _ntl_gmul(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
       /* *c = a * b */

    void _ntl_gsq(_ntl_gbigint a, _ntl_gbigint *c);
       /* *c = a * a */

    long _ntl_gsdiv(_ntl_gbigint a, long b, _ntl_gbigint *q);
       /* (*q) = floor(a/b) and a - floor(a/b)*(*q) is returned;
          error is raised if b == 0;
          if b does not divide a, then sign(*q) == sign(b) */

    void _ntl_gdiv(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *q, _ntl_gbigint *r);
       /* (*q) = floor(a/b) and (*r) = a - floor(a/b)*(*q);
          error is raised if b == 0;
          if b does not divide a, then sign(*q) == sign(b) */

    void _ntl_gmod(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *r);
       /* same as _ntl_gdiv, but only remainder is computed */

    long _ntl_gsmod(_ntl_gbigint a, long d);
       /* same as _ntl_gsdiv, but only remainder is computed */

    void _ntl_gquickmod(_ntl_gbigint *r, _ntl_gbigint b);
       /* *r = *r % b; 
	  The division is performed in place (but may sometimes
	  assumes b > 0 and *r >= 0;
          cause *r to grow by one digit) */

    void _ntl_gsaddmul(_ntl_gbigint x, long y,  _ntl_gbigint *ww);
      /* *ww += x*y */

    void _ntl_gaddmul(_ntl_gbigint x, _ntl_gbigint y,  _ntl_gbigint *ww);
      /* *ww += x*y */

    void _ntl_gssubmul(_ntl_gbigint x, long y,  _ntl_gbigint *ww);
      /* *ww -= x*y */

    void _ntl_gsubmul(_ntl_gbigint x, _ntl_gbigint y,  _ntl_gbigint *ww);
      /* *ww -= x*y */


/********************************************************************

   Shifting and bit manipulation

*********************************************************************/


    void _ntl_glshift(_ntl_gbigint n, long k, _ntl_gbigint *a);
       /* *a = sign(n) * (|n| << k);
          shift is in reverse direction for negative k */

    void _ntl_grshift(_ntl_gbigint n, long k, _ntl_gbigint *a);
       /* *a = sign(n) * (|n| >> k);
          shift is in reverse direction for negative k */
    
    long _ntl_gmakeodd(_ntl_gbigint *n);
       /*
          if (n != 0)
              *n = m;
              return (k such that n == 2 ^ k * m with m odd);
          else
              return (0); 
        */

    long _ntl_gnumtwos(_ntl_gbigint n);
        /* return largest e such that 2^e divides n, or zero if n is zero */

    long _ntl_godd(_ntl_gbigint a);
       /* returns 1 if n is odd and 0 if it is even */

    long _ntl_gbit(_ntl_gbigint a, long p);
       /* returns p-th bit of a, where the low order bit is indexed by 0;
          p out of range returns 0 */

    long _ntl_gsetbit(_ntl_gbigint *a, long p);
       /* returns original value of p-th bit of |a|, and replaces
          p-th bit of a by 1 if it was zero;
          error if p < 0 */

    long _ntl_gswitchbit(_ntl_gbigint *a, long p);
       /* returns original value of p-th bit of |a|, and switches
          the value of p-th bit of a;
          p starts counting at 0;
          error if p < 0 */


     void _ntl_glowbits(_ntl_gbigint a, long k, _ntl_gbigint *b);
        /* places k low order bits of |a| in b */ 

     long _ntl_gslowbits(_ntl_gbigint a, long k);
        /* returns k low order bits of |a| */

    long _ntl_gweights(long a);
        /* returns Hamming weight of |a| */

    long _ntl_gweight(_ntl_gbigint a);
        /* returns Hamming weight of |a| */

    void _ntl_gand(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
        /* c gets bit pattern `bits of |a|` and `bits of |b|` */

    void _ntl_gor(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
        /* c gets bit pattern `bits of |a|` inclusive or `bits of |b|` */

    void _ntl_gxor(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
        /* c gets bit pattern `bits of |a|` exclusive or `bits of |b|` */




/************************************************************************

   Comparison

*************************************************************************/

    long _ntl_gcompare(_ntl_gbigint a, _ntl_gbigint b);
       /*
          if (a > b)
              return (1);
          if (a == b)
              return (0);
          if (a < b)
              return (-1);
         */

    long _ntl_gscompare(_ntl_gbigint a, long b);
       /* single-precision version of the above */

    long _ntl_giszero (_ntl_gbigint a);
       /* test for 0 */


    long _ntl_gsign(_ntl_gbigint a);
       /* 
          if (a > 0)
              return (1);
          if (a == 0)
              return (0);
          if (a < 0)
              return (-1);
        */

    void _ntl_gabs(_ntl_gbigint *a);
       /* *a = |a| */

    void _ntl_gnegate(_ntl_gbigint *a);
       /* *a = -a */

    void _ntl_gcopy(_ntl_gbigint a, _ntl_gbigint *b);
       /* *b = a;  */

    void _ntl_gswap(_ntl_gbigint *a, _ntl_gbigint *b);
       /* swap a and b (by swaping pointers) */

    long _ntl_g2log(_ntl_gbigint a);
       /* number of bits in |a|; returns 0 if a = 0 */

    long _ntl_g2logs(long a);
        /* single-precision version of the above */


/********************************************************************

   Conversion

*********************************************************************/
        
    void _ntl_gzero(_ntl_gbigint *a);
       /* *a = 0;  */

    void _ntl_gone(_ntl_gbigint *a);
       /* *a = 1 */

    void _ntl_gintoz(long d, _ntl_gbigint *a);
       /* *a = d;  */


    void _ntl_guintoz(unsigned long d, _ntl_gbigint *a);
       /* *a = d;  space is allocated  */

    long _ntl_gtoint(_ntl_gbigint a);
       /* converts a to a long;  overflow results in value
          mod 2^{NTL_BITS_PER_LONG}. */

    unsigned long _ntl_gtouint(_ntl_gbigint a);
       /* converts a to a long;  overflow results in value
          mod 2^{NTL_BITS_PER_LONG}. */

   


    double _ntl_gdoub(_ntl_gbigint n);
       /* converts a to a double;  no overflow check */

    long _ntl_ground_correction(_ntl_gbigint a, long k, long residual);
       /* k >= 1, |a| >= 2^k, and residual is 0, 1, or -1.
          The result is what we should add to (a >> k) to round
          x = a/2^k to the nearest integer using IEEE-like rounding rules
          (i.e., round to nearest, and round to even to break ties).
          The result is either 0 or sign(a).
          If residual is not zero, it is as if x were replaced by
          x' = x + residual*2^{-(k+1)}.
          This can be used to break ties when x is exactly
          half way between two integers. */

    double _ntl_glog(_ntl_gbigint a);
       /* computes log(a), protecting against overflow */

    void _ntl_gdoubtoz(double a, _ntl_gbigint *x);
       /* x = floor(a);  */
    



/************************************************************************

   Square roots

*************************************************************************/


    long _ntl_gsqrts(long n);
       /* return floor(sqrt(n));  error raised in n < 0 */

    void _ntl_gsqrt(_ntl_gbigint n, _ntl_gbigint *r);
       /* *r =  floor(sqrt(n));  error raised in n < 0 */

/*********************************************************************
 
    Exponentiation
 
**********************************************************************/

   void _ntl_gexp(_ntl_gbigint a, long e, _ntl_gbigint *b);
       /* *b = a^e;  error raised if e < 0 */

   void _ntl_gexps(long a, long e, _ntl_gbigint *b);
       /* *b = a^e;  error raised if e < 0 */
       

/*********************************************************************

   Modular Arithmetic

   Addition, subtraction, multiplication, squaring division, inversion,
   and exponentiation modulo a positive modulus n, where all operands
   (except for the exponent in exponentiation) and results are in the
   range [0, n-1].   

   ALIAS RESTRICTION:  output parameters should not alias n

***********************************************************************/

    void _ntl_gaddmod(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint n, _ntl_gbigint *c);
       /* *c = (a + b) % n */

    void _ntl_gsubmod(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint n, _ntl_gbigint *c);
       /* *c = (a - b) % n */

    void _ntl_gsmulmod(_ntl_gbigint a, long b, _ntl_gbigint n, _ntl_gbigint *c);
       /* *c = (a * b) % n */

    void _ntl_gmulmod(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint n, _ntl_gbigint *c);
       /* *c = (a * b) % n */

    void _ntl_gsqmod(_ntl_gbigint a, _ntl_gbigint n, _ntl_gbigint *c);
       /* *c = (a ^ 2) % n */

    void _ntl_ginvmod(_ntl_gbigint a, _ntl_gbigint n, _ntl_gbigint *c);
       /* *c = (1 / a) % n; error raised if gcd(b, n) != 1 */

    void _ntl_gpowermod(_ntl_gbigint g, _ntl_gbigint e, _ntl_gbigint F,
                        _ntl_gbigint *h);

       /* *b = (a ^ e) % n; */




/**************************************************************************

   Euclidean Algorithms

***************************************************************************/
    void _ntl_ggcd(_ntl_gbigint m1, _ntl_gbigint m2, _ntl_gbigint *r);
       /* *r = greatest common divisor of m1 and m2; 
          uses binary gcd algorithm */


    void _ntl_gexteucl(_ntl_gbigint a, _ntl_gbigint *xa,
                 _ntl_gbigint b, _ntl_gbigint *xb,
                 _ntl_gbigint *d);
       /*
          *d = a * *xa + b * *xb = gcd(a, b);
          sets *d, *xa and *xb given a and b;
          uses Lehmer`s trick
        */


    long _ntl_ginv(_ntl_gbigint a, _ntl_gbigint b, _ntl_gbigint *c);
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
          error raised if a < 0 or b <= 0
        */

     long _ntl_gxxratrecon(_ntl_gbigint x, _ntl_gbigint m,  
                      _ntl_gbigint a_bound, _ntl_gbigint b_bound,
                      _ntl_gbigint *a, _ntl_gbigint *b);

        /* rational reconstruction: see doc in ZZ.txt */


        
/**********************************************************************

    Storage Allocation

    These routines use malloc and free.

***********************************************************************/

    long _ntl_gmaxalloc(_ntl_gbigint x);
       /* max allocation request, possibly rounded up a bit */

    void _ntl_gsetlength(_ntl_gbigint *v, long len);
       /* Allocates enough space to hold a len-digit number,
          where each digit has NTL_NBITS bits.
          If space must be allocated, space for one extra digit
          is always allocated. */

    void _ntl_gfree(_ntl_gbigint *x);
       /* Free's space held by x, and sets x back to 0. */


/*******************************************************************

    Special routines

********************************************************************/

long _ntl_gsize(_ntl_gbigint n);
long _ntl_gisone(_ntl_gbigint n);

long _ntl_gsptest(_ntl_gbigint a);
long _ntl_gwsptest(_ntl_gbigint a);
long _ntl_gcrtinrange(_ntl_gbigint g, _ntl_gbigint a);

void _ntl_gfrombytes(_ntl_gbigint *x, const unsigned char *p, long n);
void _ntl_gbytesfromz(unsigned char *p, _ntl_gbigint a, long nn);


long _ntl_gblock_construct_alloc(_ntl_gbigint *x, long d, long n);
void _ntl_gblock_construct_set(_ntl_gbigint x, _ntl_gbigint *y, long i);
long _ntl_gblock_destroy(_ntl_gbigint x);
long _ntl_gblock_storage(long d);


void _ntl_gcrt_struct_init(void **crt_struct, long n, _ntl_gbigint p,
                          const long *primes);
void _ntl_gcrt_struct_insert(void *crt_struct, long i, _ntl_gbigint m);
void _ntl_gcrt_struct_free(void *crt_struct);
void _ntl_gcrt_struct_eval(void *crt_struct, _ntl_gbigint *t, const long *a);
long _ntl_gcrt_struct_special(void *crt_struct);

void _ntl_grem_struct_init(void **rem_struct, long n, _ntl_gbigint p,
                          const long *primes);
void _ntl_grem_struct_free(void *rem_struct);
void _ntl_grem_struct_eval(void *rem_struct, long *x, _ntl_gbigint a);




#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
}
#endif


extern int _ntl_gmp_hack;

#define NTL_crt_struct_eval _ntl_gcrt_struct_eval
#define NTL_crt_struct_free _ntl_gcrt_struct_free
#define NTL_crt_struct_init _ntl_gcrt_struct_init
#define NTL_crt_struct_insert _ntl_gcrt_struct_insert
#define NTL_crt_struct_special _ntl_gcrt_struct_special
#define NTL_rem_struct_eval _ntl_grem_struct_eval
#define NTL_rem_struct_free _ntl_grem_struct_free
#define NTL_rem_struct_init _ntl_grem_struct_init
#define NTL_verylong _ntl_gbigint
#define NTL_z2log _ntl_g2log
#define NTL_zabs _ntl_gabs
#define NTL_zadd _ntl_gadd
#define NTL_zaddmod _ntl_gaddmod
#define NTL_zand _ntl_gand
#define NTL_zbit _ntl_gbit
#define NTL_zblock_construct_alloc _ntl_gblock_construct_alloc
#define NTL_zblock_construct_set _ntl_gblock_construct_set
#define NTL_zblock_destroy _ntl_gblock_destroy
#define NTL_zblock_storage _ntl_gblock_storage
#define NTL_zbytesfromz _ntl_gbytesfromz
#define NTL_zcompare _ntl_gcompare
#define NTL_zcopy _ntl_gcopy
#define NTL_zcrtinrange _ntl_gcrtinrange
#define NTL_zdiv _ntl_gdiv
#define NTL_zdoub _ntl_gdoub
#define NTL_zdoubtoz _ntl_gdoubtoz
#define NTL_zexp _ntl_gexp
#define NTL_zexps _ntl_gexps
#define NTL_zexteucl _ntl_gexteucl
#define NTL_zfree _ntl_gfree
#define NTL_zfrombytes _ntl_gfrombytes
#define NTL_zgcd _ntl_ggcd
#define NTL_zintoz _ntl_gintoz
#define NTL_zinv _ntl_ginv
#define NTL_zinvmod _ntl_ginvmod
#define NTL_zisone _ntl_gisone
#define NTL_ziszero _ntl_giszero
#define NTL_zlog _ntl_glog
#define NTL_zlowbits _ntl_glowbits
#define NTL_zlshift _ntl_glshift
#define NTL_zmakeodd _ntl_gmakeodd
#define NTL_zmod _ntl_gmod
#define NTL_zmul _ntl_gmul
#define NTL_zmulmod _ntl_gmulmod
#define NTL_znegate _ntl_gnegate
#define NTL_znumtwos _ntl_gnumtwos
#define NTL_zodd _ntl_godd
#define NTL_zone _ntl_gone
#define NTL_zor _ntl_gor
#define NTL_zpowermod _ntl_gpowermod
#define NTL_zquickmod _ntl_gquickmod
#define NTL_zround_correction _ntl_ground_correction
#define NTL_zrshift _ntl_grshift
#define NTL_zsadd _ntl_gsadd
#define NTL_zscompare _ntl_gscompare
#define NTL_zsdiv _ntl_gsdiv
#define NTL_zsetbit _ntl_gsetbit
#define NTL_zmaxalloc _ntl_gmaxalloc
#define NTL_zsetlength _ntl_gsetlength
#define NTL_zsign _ntl_gsign
#define NTL_zsize _ntl_gsize
#define NTL_zslowbits _ntl_gslowbits
#define NTL_zsmod _ntl_gsmod
#define NTL_zsmul _ntl_gsmul
#define NTL_zsmulmod _ntl_gsmulmod
#define NTL_zsptest _ntl_gsptest
#define NTL_zsq _ntl_gsq
#define NTL_zsqmod _ntl_gsqmod
#define NTL_zsqrt _ntl_gsqrt
#define NTL_zsqrts _ntl_gsqrts
#define NTL_zsub _ntl_gsub
#define NTL_zsubmod _ntl_gsubmod
#define NTL_zsubpos _ntl_gsubpos
#define NTL_zswap _ntl_gswap
#define NTL_zswitchbit _ntl_gswitchbit
#define NTL_ztoint _ntl_gtoint
#define NTL_ztouint _ntl_gtouint
#define NTL_zuintoz _ntl_guintoz
#define NTL_zweight _ntl_gweight
#define NTL_zweights _ntl_gweights
#define NTL_zwsptest _ntl_gwsptest
#define NTL_zxor _ntl_gxor
#define NTL_zxxratrecon _ntl_gxxratrecon
#define NTL_zzero _ntl_gzero

#define NTL_zsaddmul _ntl_gsaddmul
#define NTL_zaddmul _ntl_gaddmul
#define NTL_zssubmul _ntl_gssubmul
#define NTL_zsubmul _ntl_gsubmul

#define NTL_GMP_LIP

