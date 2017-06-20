
/*************************************************************


   Assembly code support for computing the high-order
   word of a word * word product (unsigned).

   Note that these typically only make a significant difference
   on some 64-bit machines, as on 32-bit machines, the "long long"
   solution is usually just as good.

   These code sequences were extracted from a recent version of
   the file longlong.h from gmp.  Copyright notice follows:

      Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003
      Free Software Foundation, Inc.

      This file is free software; you can redistribute it and/or modify
      it under the terms of the GNU Lesser General Public License as published by
      the Free Software Foundation; either version 2.1 of the License, or (at your
      option) any later version.

      This file is distributed in the hope that it will be useful, but
      WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
      or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
      License for more details.

      You should have received a copy of the GNU Lesser General Public License
      along with this file; see the file COPYING.LIB.  If not, write to
      the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
      MA 02111-1307, USA.
   

   
*************************************************************/






#if (defined(__GNUC__) && (__GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 7)))

//  To simplify things, we require gcc v2.7 or higher.



// ------ POWERPC ------

#if defined (_ARCH_PPC) || defined (__powerpc__) || defined (__POWERPC__) \
  || defined (__ppc__) || defined(__ppc64__)                                                   \
  || (defined (PPC) && ! defined (CPU_FAMILY)) /* gcc 2.7.x GNU&SysV */       \
  || (defined (PPC) && defined (CPU_FAMILY)    /* VxWorks */                  \
         && CPU_FAMILY == PPC)

#if (NTL_BITS_PER_LONG == 32)

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi;
    __asm__ ("mulhwu %0,%1,%2" : "=r" (hi) : "%r" (a), "r" (b));
   return hi;
} 

#elif (NTL_BITS_PER_LONG == 64)


static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi;
    __asm__ ("mulhdu %0,%1,%2" : "=r" (hi) : "%r" (a), "r" (b));
   return hi;
} 

#endif

#endif




// ------ ALPHA ------

#if (defined (__alpha) && NTL_BITS_PER_LONG == 64)

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi;
   __asm__ ("umulh %r1,%2,%0" : "=r" (hi) : "%rJ" (a), "rI" (b));
   return hi;
} 

#endif


// ------ IA64 ------

#if (defined (__ia64) && NTL_BITS_PER_LONG == 64)

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi;
   __asm__ ("xma.hu %0 = %1, %2, f0" : "=f" (hi) : "f" (a), "f" (b));
   return hi;
} 


#endif


//  ------ x86 ------

#if ((defined (__i386__) || defined (__i486__)) && NTL_BITS_PER_LONG == 32)

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi, lo;
   __asm__ ("mull %3" : "=a" (lo), "=d" (hi) : "%0" (a), "rm" (b));
   
   return hi;
} 

#endif


// ------ x86-64 ------

#if (defined (__x86_64__) && NTL_BITS_PER_LONG == 64)

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi, lo;
   __asm__ ("mulq %3" : "=a" (lo), "=d" (hi) : "%0" (a), "rm" (b));
   
   return hi;
} 


#endif


// ------ MIPS ------

#if (defined (__mips))

#if (NTL_BITS_PER_LONG == 32)

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi, lo;
   __asm__ ("multu %2,%3" : "=l" (lo), "=h" (hi) : "d" (a), "d" (b));
   return hi;
} 



#elif (NTL_BITS_PER_LONG == 64)


static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi, lo;
   __asm__ ("dmultu %2,%3" : "=l" (lo), "=h" (hi) : "d" (a), "d" (b));
   return hi;
}


#endif

#endif


//  -------- SPARC --------


#if (defined (__sparc__) && NTL_BITS_PER_LONG == 32)

#if (defined (__sparc_v9__) || defined (__sparcv9) || \
     defined (__sparc_v8__) || defined (__sparcv8) || defined (__sparclite__))

static inline unsigned long MulHiUL(unsigned long a, unsigned long b)
{
   unsigned long hi, lo;
   __asm__ ("umul %2,%3,%1;rd %%y,%0" : "=r" (hi), "=r" (lo) : "r" (a), "r" (b));
   return hi;
}


#endif





#endif



#endif // __GNUC__ 
