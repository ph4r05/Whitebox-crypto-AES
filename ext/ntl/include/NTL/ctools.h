
#ifndef NTL_ctools__H
#define NTL_ctools__H

#include <NTL/config.h>
#include <NTL/mach_desc.h>


/*
 * Resolve double-word integer types.
 *
 * Unfortunately, there is no "standard" way to do this.
 * On 32-bit machines, 'long long' usually works (but not
 * on MSVC++ or BORLAND), and on 64-bit machines, there is
 * no standard.  However, most compilers do offer *some*
 * non-standard double-word type.  
 *
 * Note that C99 creates a standard header <stdint.h>,
 * but it is not clear how widely this is implemented yet,
 * and for example, GCC does not provide a type int128_t 
 * in <stdint.h> on 64-bit machines.
 */


#if (defined(NTL_LONG_LONG_TYPE))

#define NTL_LL_TYPE NTL_LONG_LONG_TYPE

#elif (NTL_BITS_PER_LONG == 64 && defined(__GNUC__))

#define NTL_LL_TYPE __int128_t

#elif (NTL_BITS_PER_LONG == 32 && (defined(_MSC_VER) || defined(__BORLANDC__)))

#define NTL_LL_TYPE __int64

#elif (NTL_BITS_PER_LONG == 64 && (defined(_MSC_VER) || defined(__BORLANDC__)))

#define NTL_LL_TYPE __int128

#endif

#if (!defined(NTL_LL_TYPE))

#define NTL_LL_TYPE long long

#endif



#if (defined(NTL_UNSIGNED_LONG_LONG_TYPE))

#define NTL_ULL_TYPE NTL_UNSIGNED_LONG_LONG_TYPE

#elif (NTL_BITS_PER_LONG == 64 && defined(__GNUC__))

#define NTL_ULL_TYPE __uint128_t

#elif (NTL_BITS_PER_LONG == 32 && (defined(_MSC_VER) || defined(__BORLANDC__)))

#define NTL_ULL_TYPE unsigned __int64

#elif (NTL_BITS_PER_LONG == 64 && (defined(_MSC_VER) || defined(__BORLANDC__)))

#define NTL_ULL_TYPE unsigned __int128

#endif

#if (!defined(NTL_ULL_TYPE))

#define NTL_ULL_TYPE unsigned long long

#endif


/********************************************************/







#define NTL_OVFBND (1L << (NTL_BITS_PER_LONG-4))

/*
 * NTL_OVFBND is the general bound used throughout NTL to keep various
 * integer values comfortably bounded away from an integer overflow
 * condition.  Do not change this value!
 */





#if ((NTL_BITS_PER_SIZE_T-1) < (NTL_BITS_PER_LONG-4))
#define NTL_OVFBND1 (1L << (NTL_BITS_PER_SIZE_T-1))
#else
#define NTL_OVFBND1 NTL_OVFBND
#endif

/*
 * NTL_OVFBND1 is a smaller bound than NTL_OVF when size_t is
 * narrower than long.  This prevents overflow on calls to malloc
 * and realloc.
 */






#define NTL_OVERFLOW(n, a, b) \
   (((b) >= NTL_OVFBND) || (((long) (n)) > 0 && (((a) >= NTL_OVFBND) || \
    (((long) (n)) >= (NTL_OVFBND-((long)(b))+((long)(a))-1)/((long)(a))))))

/*
 * NTL_OVERFLOW(n, a, b) returns 1 if n*a + b >= NTL_OVFBND,
 * and returns 0 otherwise.  The value n is effectively treated as type long,
 * while the values a and b may be *any* integral type.  It is assumed that
 * n >= 0, a > 0, and b >= 0.  Care is taken to ensure that overflow does
 * not occur. If a and b are constants, and n has no side effects,
 * a good optimizing compiler will * translate this into a single test 
 * of the form n >= c, where c is a constant.
 */






#define NTL_OVERFLOW1(n, a, b) \
   (((b) >= NTL_OVFBND1) || (((long) (n)) > 0 && (((a) >= NTL_OVFBND1) || \
    (((long) (n)) >= (NTL_OVFBND1-((long)(b))+((long)(a))-1)/((long)(a))))))

/*
 * NTL_OVERFLOW1 is the same as NTL_OVERFLOW, except that it uses the
 * bound NTL_OVFBND1 instead of NTL_OVFBND.
 */





#define NTL_MALLOC(n, a, b) \
   (NTL_OVERFLOW1(n, a, b) ? ((void *) 0) : \
    ((void *) malloc(((long)(n))*((long)(a)) + ((long)(b)))))

/*
 * NTL_MALLOC(n, a, b) returns 0 if a*n + b >= NTL_OVFBND1, and otherwise
 * returns malloc(n*a + b). 
 * The programmer must ensure that the name "malloc" is visible
 * at the point in the source code where this macro is expanded.
 */






#define NTL_SNS_MALLOC(n, a, b) \
   (NTL_OVERFLOW1(n, a, b) ? ((void *) 0) : \
    ((void *) NTL_SNS malloc(((long)(n))*((long)(a)) + ((long)(b)))))

/*
 * NTL_SNS_MALLOC is the same as NTL_MALLOC, except that the call
 * to malloc is prefixed by NTL_SNS.
 */








#define NTL_REALLOC(p, n, a, b) \
   (NTL_OVERFLOW1(n, a, b) ? ((void *) 0) : \
    ((void *) realloc((p), ((long)(n))*((long)(a)) + ((long)(b)))))

/*
 * NTL_REALLOC(n, a, b) returns 0 if a*n + b >= NTL_OVFBND1, and otherwise
 * returns realloc(p, n*a + b).
 * The programmer must ensure that the name "realloc" is visible
 * at the point in the source code where this macro is expanded.
 */






#define NTL_SNS_REALLOC(p, n, a, b) \
   (NTL_OVERFLOW1(n, a, b) ? ((void *) 0) : \
    ((void *) NTL_SNS realloc((p), ((long)(n))*((long)(a)) + ((long)(b)))))

/*
 * NTL_SNS_REALLOC is the same as NTL_REALLOC, except that the call
 * to realloc is prefixed by NTL_SNS.
 */





#define NTL_MAX_ALLOC_BLOCK (40000)

/*
 * NTL_MAX_ALLOC_BLOCK is the number of bytes that are allocated in
 * a single block in a number of places throughout NTL (for
 * vec_ZZ_p, ZZVec, vec_GF2X, and GF2XVec).
 */


#define NTL_ULONG_TO_LONG(a) \
   ((((unsigned long) a) >> (NTL_BITS_PER_LONG-1)) ? \
    (((long) (((unsigned long) a) - ((unsigned long) NTL_MIN_LONG))) + \
       NTL_MIN_LONG) : \
    ((long) a))

/* 
 * This macro converts from unsigned long to signed long.  It is portable
 * among platforms for which a long has a 2's complement representation
 * of the same width as an unsigned long.  While it avoids assumptions
 * about the behavior of non-standard conversions,  a good optimizing
 * compiler should turn it into the identity function.
 */


#define NTL_UINT_TO_INT(a) \
   ((((unsigned int) a) >> (NTL_BITS_PER_INT-1)) ? \
    (((int) (((unsigned int) a) - ((unsigned int) NTL_MIN_INT))) + \
       NTL_MIN_INT) : \
    ((int) a))

/* 
 * This macro converts from unsigned int to signed int.  It is portable
 * among platforms for which an int has a 2's complement representation
 * of the same width as an unsigned int.  While it avoids assumptions
 * about the behavior of non-standard conversions,  a good optimizing
 * compiler should turn it into the identity function.
 */


#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
extern "C" {
#endif

long _ntl_IsFinite(double *p);
/* This forces a double into memory, and tests if it is "normal";
   that means, not NaN, not +/- infinity, not denormalized, etc.
   Forcing into memory is sometimes necessary on machines 
   with "extended" double precision registers (e.g., Intel x86s)
   to force the standard IEEE format. */

void _ntl_ForceToMem(double *p);
/* This is do-nothing routine that has the effect of forcing
   a double into memory (see comment above). */

double _ntl_ldexp(double x, long e);

void _ntl_abort(void);
/* This is the routine called by NTL to abort a program in case of error. */

void _ntl_abort_cxx_callback(void);
/* This is a C++ function (implemented in tools.c) that is
   used to implement the callback mechanism.  The issue here
   is that I don't want a C function to call a C++ function
   via a function pointer.  This could potentially be problematic. */

   
   
#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
}
#endif

#endif

