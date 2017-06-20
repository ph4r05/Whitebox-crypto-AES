
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <NTL/config.h>


#ifndef NTL_GMP_LIP


int main()
{
   fprintf(stderr, "NTL_GMP_LIP flag not set\n");

   return 0;
}



#else


#include <gmp.h>
#include <NTL/mach_desc.h>

void print2k(FILE *f, long k, long bpl)
{
   long m, l;
   long first;

   if (k <= 0) {
      fprintf(f, "((double) 1.0)");
      return;
   }

   m = bpl - 2;
   first = 1;

   fprintf(f, "(");

   while (k > 0) {
      if (k > m)
         l = m;
      else
         l = k;

      k = k - l;


      if (first)
         first = 0;
      else
         fprintf(f, "*");

      fprintf(f, "((double)(1L<<%ld))", l);
   }

   fprintf(f, ")");
}

void Error(const char *s)
{
   fprintf(stderr, "%s\n", s);
   abort();
}


int main()
{
   long bpl;
   long ntl_zz_nbits, ntl_wsp_nbits, ntl_sp_nbits;

   fprintf(stderr, "NTL_GMP_LIP flag set\n");

   bpl = NTL_BITS_PER_LONG;


   /*
    * We require that the number of bits per limb quantity correspond to the
    * number of bits of a long, or possibly a "long long" that is twice as wide
    * as a long.  These restrictions will almost certainly be satisfied, unless
    * GMP is installed using the newly proposed "nail" option.
    */

   ntl_zz_nbits = 0;

   if (sizeof(mp_limb_t) == sizeof(long) && mp_bits_per_limb == bpl)
      ntl_zz_nbits = bpl;
   else if (sizeof(mp_limb_t) == 2*sizeof(long) && mp_bits_per_limb == 2*bpl)
      ntl_zz_nbits = 2*bpl;
   else
      Error("sorry...this is a funny gmp");

   if (sizeof(mp_size_t) != sizeof(long) &&
       sizeof(mp_size_t) != sizeof(int))

   Error("sorry...this is a funny gmp");

   ntl_wsp_nbits = bpl - 2;

   ntl_sp_nbits = NTL_NBITS_MAX;

   if (sizeof(mp_size_t) < sizeof(long)) {
      printf("#define NTL_SMALL_MP_SIZE_T\n");
      fprintf(stderr, "setting NTL_SMALL_MP_SIZE_T\n");
   }

   fprintf(stderr, "NTL_ZZ_NBITS = %ld\n", ntl_zz_nbits);
   fprintf(stderr, "NTL_WSP_NBITS = %ld\n", ntl_wsp_nbits);
   fprintf(stderr, "NTL_SP_NBITS = %ld\n", ntl_sp_nbits);

   printf("#define NTL_ZZ_NBITS (%ld)\n",  ntl_zz_nbits);

   printf("#define NTL_ZZ_FRADIX ");
   print2k(stdout, ntl_zz_nbits, bpl);
   printf("\n");

   return 0;
}

#endif
