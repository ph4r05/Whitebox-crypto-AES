

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <NTL/config.h>

#ifdef NTL_GMP_HACK

#include <gmp.h>
#include <NTL/mach_desc.h>


int gcd(int a, int b)
{
   int u, v, t, x;

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

int pow2(int a)
{
   int m;
   int k;

   m = 1;
   k = 0;

   while (m < a) {
      m = 2*m;
      k++;
   }

   if (m == a) 
      return k;
   else
      return -1;
}
  
   

const char *accum[2] = { "", "yy | " };

void lip_to_gmp(int A, int B)
{
   int d, na, nb;
   int i, j, r;
   int xx, yy;
   int shamt;

   int iter, na2, nb2;

   d = gcd(A, B);
   na = B/d;
   nb = A/d;

   na2 = pow2(na);
   nb2 = pow2(nb);

   printf("static\n");
   printf("void lip_to_gmp(const long *x, mp_limb_t *y, long n)\n");
   printf("{\n");
   printf("   long r, q, xx;\n");
   printf("   mp_limb_t yy;\n\n");

   if (na2 != -1) {
      printf("   r = n & %d;\n", na-1);
      printf("   q = n >> %d;\n", na2);
   }
   else {
      printf("   r = n % %d;\n", na);
      printf("   q = n / %d;\n", na);
   }

   printf("\n");

   printf("   if (q > 0) {\n");

   if (na2 != -1) 
      printf("      x += (q << %d);\n", na2);
   else 
      printf("      x += (q * %d);\n", na);

   if (nb2 != -1) 
      printf("      y += (q << %d);\n", nb2);
   else if (pow2(nb+1) != -1)
      printf("      y += (q << %d) - q;\n", pow2(nb+1));
   else
      printf("      y += (q * %d);\n", nb);

   printf("   }\n");

   printf("\n");
   printf("   if (r > 0) {\n");
   printf("      yy = 0;\n");
   printf("\n");

   printf("      switch (r-1) {\n");
      

   xx = 0;
   yy = 0;

   for (i = 1; ; i++) {
      j = (i*A)/B;
      r = (i*A)%B;


      if (xx) {
         printf("         ");
         printf("yy = %s (((mp_limb_t)(xx)) << %d);\n", accum[yy], B-r);
         yy = 1;
      }

      printf("      ");
      printf("case %d:\n", na-1-i);

      shamt = A-B+r;

      printf("         ");

      if (shamt >= 0) {
         printf("y[%d] = ", nb-1-j);
      }
      else {
         printf("yy = ");
      }

      printf("%s", accum[yy]);

      if (shamt == 0) {
         printf("x[0];\n");
         break;
      }
      else if (shamt < 0) {
         printf("(((mp_limb_t)(x[%d])) << %d);\n", na-1-i, -shamt);
         xx = 0;
         yy = 1;
      }
      else {
         printf("((xx=x[%d]) >> %d);\n", na-1-i, shamt);
         xx = 1;
         yy = 0;
      }
   }

   printf("      }\n\n");
   printf("      n -= r;\n");
   printf("   }\n\n");
   printf("   while (n > 0) {\n");
   printf("      x -= %d;\n", na);
   printf("      n -= %d;\n", na);
   printf("      y -= %d;\n", nb);
   printf("\n");


   xx = 0;
   yy = 0;

   for (i = 0; ; i++) {
      j = (i*A)/B;
      r = (i*A)%B;

      
      if (xx) {
         printf("      ");
         printf("yy = %s (((mp_limb_t)(xx)) << %d);\n", accum[yy], B-r);
         yy = 1;
      }

      shamt = A-B+r;

      printf("      ");

      if (shamt >= 0) {
         printf("y[%d] = ", nb-1-j);
      }
      else {
         printf("yy = ");
      }

      printf("%s", accum[yy]);

      if (shamt == 0) {
         printf("x[0];\n");
         break;
      }
      else if (shamt < 0) {
         printf("(((mp_limb_t)(x[%d])) << %d);\n", na-1-i, -shamt);
         xx = 0;
         yy = 1;
      }
      else {
         printf("((xx=x[%d]) >> %d);\n", na-1-i, shamt);
         xx = 1;
         yy = 0;
      }
   }

   printf("   }\n");
   printf("}\n");
}


void gmp_to_lip(int A, int B, int alt)
{
   int d, na, nb;
   int i, j, r;
   int shamt;

   int na2, nb2;

   d = gcd(A, B);
   na = B/d;
   nb = A/d;

   na2 = pow2(na);
   nb2 = pow2(nb);

   printf("static\n");
   if (!alt)
      printf("void gmp_to_lip(long *x, const mp_limb_t *y, long n)\n");
   else
      printf("void gmp_to_lip1(long *x, const mp_limb_t *y, long n)\n");
   printf("{\n");
   printf("   long r, q; unsigned long xx;\n");
   printf("   mp_limb_t yy;\n\n");

   if (na2 != -1) {
      printf("   r = n & %d;\n", na-1);
      printf("   q = n >> %d;\n", na2);
   }
   else {
      printf("   r = n % %d;\n", na);
      printf("   q = n / %d;\n", na);
   }

   printf("\n");

   printf("   if (q > 0) {\n");

   if (na2 != -1)
      printf("      x += (q << %d);\n", na2);
   else
      printf("      x += (q * %d);\n", na);

   if (nb2 != -1)
      printf("      y += (q << %d);\n", nb2);
   else if (pow2(nb+1) != -1)
      printf("      y += (q << %d) - q;\n", pow2(nb+1));
   else
      printf("      y += (q * %d);\n", nb);

   printf("   }\n");

   printf("\n");


   printf("   if (r > 0) {\n");

   if (!alt) {
      if (na - nb == 1)
         printf("      yy = y[r-1];\n");
      else if (na2 != -1) 
         printf("      yy = y[((r*%d + %d) >> %d) - 1];\n", nb, na-1, na2);
      else
         printf("      yy = y[((r*%d + %d) / %d) - 1];\n", nb, na-1, na);
   
   }
   else
      printf("      yy = 0;\n");

   printf("\n");

   printf("      switch (r-1) {\n");
      


   for (i = 1; ; i++) {
      j = (i*A)/B;
      r = (i*A)%B;


      printf("      ");
      printf("case %d:\n", na-1-i);

      shamt = A-B+r;


      if (shamt > 0) {
         printf("         ");
         printf("xx = ((unsigned long)(yy)) << %d;\n", shamt);
         printf("         ");
         printf("x[%d] = (xx | ((unsigned long)((yy = y[%d]) >> %d))) & NTL_RADIXM;\n", 
                na-1-i, nb-2-j, 2*B-r-A);
      }
      else if (shamt == 0) {
         printf("         ");
         printf("x[0] = ((unsigned long)(yy)) & NTL_RADIXM;\n");
         break;
      }
      else {
         printf("         ");
         printf("x[%d] = ((unsigned long)(yy >> %d)) & NTL_RADIXM;\n", na-1-i, -shamt);
      }
   }

   printf("      }\n\n");
   printf("      n -= r;\n");
   printf("   }\n\n");
   printf("   while (n > 0) {\n");
   printf("      x -= %d;\n", na);
   printf("      n -= %d;\n", na);
   printf("      y -= %d;\n", nb);
   printf("\n");

   printf("      yy = y[%d];\n", nb-1);

   for (i = 0; ; i++) {
      j = (i*A)/B;
      r = (i*A)%B;


      shamt = A-B+r;

      if (shamt > 0) {
         printf("      ");
         printf("xx = ((unsigned long)(yy)) << %d;\n", shamt);
         printf("      ");
         printf("x[%d] = (xx | ((unsigned long)((yy = y[%d]) >> %d))) & NTL_RADIXM;\n", 
                na-1-i, nb-2-j, 2*B-r-A);
      }
      else if (shamt == 0) {
         printf("      ");
         printf("x[0] = ((unsigned long)(yy)) & NTL_RADIXM;\n");
         break;
      }
      else {
         printf("      ");
         printf("x[%d] = ((unsigned long)(yy >> %d)) & NTL_RADIXM;\n", na-1-i, -shamt);
      }
   }


   printf("   }\n");
   printf("}\n");

}

void rdup(int a, int b)
{
   printf("((");
   if (pow2(a) != -1) {
      printf("(x << %d)", pow2(a)); 
   }
   else if (pow2(a+1) != -1) {
      printf("((x << %d) - x)", pow2(a+1));
   }
   else {
      printf("(x*%d)", a);
   }
   
   printf(" + %d)", b-1);

   if (pow2(b) != -1) {
      printf(" >> %d)", pow2(b));
   }
   else {
      printf(" / %d)", b);
   }
}


void G_TO_L(int A, int B)
{
   int d, na, nb;

   d = gcd(A, B);
   na = B/d;
   nb = A/d;

   printf("#define G_TO_L(x) ");
   
   rdup(na, nb);

   printf("\n");
}

void L_TO_G(int A, int B)
{
   int d, na, nb;

   d = gcd(A, B);
   na = B/d;
   nb = A/d;

   printf("#define L_TO_G(x) ");
   
   rdup(nb, na);

   printf("\n");
}


void L_TO_G_CHECK(int BPL, int BPI)
{
   if (BPL != BPI) {
      printf("#define L_TO_G_CHECK_LEN\n");
   }
      
}


void Error(const char *s)
{
   fprintf(stderr, "%s\n", s);
   abort();
}


int main()
{
   mpz_t tt;
   int A, B, BPL, BPI;

   A = NTL_NBITS_MAX;

   /*
    * We compute B as the number of bits of a gmp limb.
    * We require that this quantity correspond to the number of bits
    * of a long, or possibly a "long long" that is twice as
    * wide as a long.  These restrictions may not be entirely 
    * necessary, but they are satisfied on all platforms that I know of.
    */

   if (sizeof(mp_limb_t)==sizeof(long) &&
            mp_bits_per_limb == NTL_BITS_PER_LONG)

      B = NTL_BITS_PER_LONG;

   else if (sizeof(mp_limb_t) == 2*sizeof(long) &&
            mp_bits_per_limb == 2*NTL_BITS_PER_LONG)

      B = 2*NTL_BITS_PER_LONG;

   else
      Error("sorry...this is a funny gmp");

   /*
    * The following test is a bit redundant, but it doesn't hurt.
    */

   if (A >= B) Error("sorry...this is a funny gmp");



   /*
    * Next, we check if either the _mp_size field of an mpz struct
    * or the type mp_size_t is narrower than type "long".
    * This is done to enable some overflow checks.
    * For simplicity, we require that the sizeof
    * these types is that of an int or a long, and we also make
    * the somewhat DIRTY assumption that this sizeof value implies
    * a corresponding bit count.  Since this assumption is true
    * on all platforms that I know of, and since this assumption
    * only affects some overflow tests, it seems reasonable.
    */

   if (sizeof(tt->_mp_size) == sizeof(int))
      BPI = NTL_BITS_PER_INT;
   else if (sizeof(tt->_mp_size) == sizeof(long))
      BPI = NTL_BITS_PER_LONG;
   else
      Error("sorry...this is a funny gmp");

   if (sizeof(mp_size_t) != sizeof(int) && sizeof(mp_size_t) != sizeof(long))
      Error("sorry...this is a funny gmp");

   if (sizeof(mp_size_t) < sizeof(tt->_mp_size))
      BPI = NTL_BITS_PER_INT;


   BPL = NTL_BITS_PER_LONG;

   fprintf(stderr, "\ngmp looks OK...%d, %d, %d, %d.\n", A, B, BPI, BPL);
   fprintf(stderr, "generating file lip_gmp_aux.c.\n");


   lip_to_gmp(A, B);
   printf("\n");

   gmp_to_lip(A, B, 0);
   printf("\n");
   gmp_to_lip(A, B, 1);
   printf("\n");

   G_TO_L(A, B);
   L_TO_G(A, B);
   L_TO_G_CHECK(BPL, BPI);

   printf("#define HAVE_LIP_GMP_AUX\n");

   return 0;
}

#else

int main()
{
   fprintf(stderr, "NTL_GMP_HACK flag not set.\n");
   return 0;
}

#endif
