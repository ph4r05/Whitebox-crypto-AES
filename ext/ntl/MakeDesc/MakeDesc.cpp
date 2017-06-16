

#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>


#include <NTL/version.h>


#if (defined(__GNUC__) && (defined(__i386__) || defined(__i486__) || defined(__i586__)))


#define AutoFix (1)

#else

#define AutoFix (0)

#endif


int val_int(int x);
unsigned int val_uint(unsigned int x);

long val_long(long x);
unsigned long val_ulong(unsigned long x);

size_t val_size_t(size_t x);

double val_double(double x);

void touch_int(int* x);
void touch_uint(unsigned int* x);

void touch_long(long* x);
void touch_ulong(unsigned long* x);

void touch_size_t(size_t* x);

void touch_double(double* x);




double power2(long k)
{
   long i;
   double res;

   res = 1;

   for (i = 1; i <= k; i++)
      res = res * 2;

   return res;
}


long DoubleRounding(long dp)
{
   double a = power2(dp-1) + 1;
   double b = (power2(dp)-1)/power2(dp+1);
   register double x = a + b;
   double y = x;

   touch_double(&y);

   if (y != power2(dp-1) + 1)
      return 1;
   else 
      return 0; 
}



long DoublePrecision(void)
{
   double eps, one, res;
   long k;

   one = val_double(1.0);
   eps = val_double(1.0);

   k = 0;

   do {
      double tmp;

      k++;
      eps *= 1.0/2.0;
      tmp = 1.0 + eps;
      touch_double(&tmp);
      res = tmp - one;
   } while (res == eps);

   return k;
}

long DoublePrecision1(void)
{
   double eps, one, res;
   long k;

   one = val_double(1.0);
   eps = val_double(1.0);

   k = 0;

   do {
      register double tmp;

      k++;
      eps *= 1.0/2.0;
      tmp = 1.0 + eps;
      res = tmp - one;
   } while (res == eps);

   return k;
}


union d_or_rep {
   double d;
   unsigned long rep[2];
};

long RepTest(void)
{
   union d_or_rep v;

   if (sizeof(double) != 2*sizeof(long))
      return 0;

   v.rep[0] = v.rep[1] = 0;

   v.d = 565656565656.0;

   if (v.rep[0] == 0x42607678 && v.rep[1] == 0x46f30000)
      return 1;
   else if (v.rep[1] == 0x42607678 && v.rep[0] == 0x46f30000)
      return -1;
   else
      return 0;
}

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






void print_mul_body(FILE *f, long n1, long k, long fn,
                    long half_flag, long short_flag)
{

  long n, i, chop, r;
  unsigned long mask, mask2;

  if (half_flag)
    n = n1/2;
  else
    n = n1;

  chop = n % k; /* first block */
  if (chop == 0)
    chop = k;
  r = n - k;
  mask = (1UL << k) - 1UL;

  fprintf(f, "\n\n#define NTL_");
  if (half_flag) fprintf(f, "HALF_");
  if (short_flag) fprintf(f, "SHORT_");

  fprintf(f, "BB_MUL_CODE%ld \\\n", fn);

  
  if (fn > 0) /* Mul1/AddMul1 */
    {
      fprintf(f, "   long i;\\\n");
      fprintf(f, "   _ntl_ulong carry = 0, b;\\\n");
    }
  fprintf(f, "   _ntl_ulong hi, lo, t;\\\n");
  fprintf(f, "   _ntl_ulong A[%ld];\\\n", 1L << k);
  fprintf(f, "   A[0] = 0;\\\n");

  fprintf(f, "   A[1] = a;\\\n");

  for (i = 2; i < (1L << k); i++)
    {
      if (i % 2 == 0)
        fprintf(f, "   A[%ld] = A[%ld] << 1;\\\n", i, i / 2);
      else
        fprintf(f, "   A[%ld] = A[%ld] ^ A[1];\\\n", i, i - 1);
    }

  if (fn > 0)
    {
      fprintf(f, "   for (i = 0; i < sb; i++) {\\\n");
      fprintf(f, "      b = bp[i];\\\n");
      fprintf(f, "   ");
    }

  fprintf(f, "   lo = A[b & %lu]; ", mask);
  fprintf(f, "t = A[(b >> %ld) & %lu]; ", k, mask);
  fprintf(f, "hi = t >> %ld; lo ^= t << %ld;\\\n", n1-k, k);

  for (i = 2*k; i < n - chop; i += k) {
    if (fn > 0) fprintf(f, "   ");
    fprintf(f, "   t = A[(b >> %ld) & %lu]; ", i, mask);
    fprintf(f, "hi ^= t >> %ld; lo ^= t << %ld;\\\n", n1-i, i);
  }

  if (fn > 0) fprintf(f, "   ");
  fprintf(f, "   t = A[b >> %ld]; ", n-chop);
  fprintf(f, "hi ^= t >> %ld; lo ^= t << %ld;\\\n", n1-i, i);

  mask = 0;
  for (i = 0; i < n; i += k)
    mask |= 1UL << i;
  mask = ~mask;
  if (half_flag) mask &= (1UL << n) - 1UL;
  mask2 = mask;
  
  if (!short_flag) {
    for (i = 1; i < k; i++)
      {
	if (fn > 0) fprintf(f, "   ");

        if (i == 1)
	   fprintf(f, "   if (a >> %ld) ", n1-i);
        else
           fprintf(f, "   if ((a >> %ld) & 1) ", n1-i);

	/* bit n1-i from a was not considered in blocks of
           k bits from b for index j >= i */
	fprintf(f, "hi ^= ((b & 0x%lxUL) >> %ld);\\\n", mask2, i);
	mask2 = (mask2 << 1) & mask;
      }
   }

  if (fn > 0) fprintf(f, "   ");

  if (fn == 0)
    {
      fprintf(f, "   c[0] = lo; ");
      fprintf(f, "   c[1] = hi;\\\n");
    }
  else if (fn == 1 || fn == 3)
    {
      fprintf(f, "   cp[i] = carry ^ lo; ");
      fprintf(f, "   carry = hi;\\\n");
    }
  else if (fn == 2)
    {
      fprintf(f, "   cp[i] ^= (carry ^ lo); ");
      fprintf(f, "   carry = hi;\\\n");
    }
  if (fn > 0)
    {

      fprintf(f, "   }\\\n");

      if (fn == 1 || fn == 3) 
        fprintf(f, "   cp[sb] = carry;\\\n");
      else
        fprintf(f, "   cp[sb] ^= carry;\\\n");

    }
  fprintf(f, "\n\n\n");


}





/*
 * This generates anternative code that runs significantly faster 
 * on some machines, like a PowerPC (and probably other RISC machines).
 * It makes it easier for the compiler to schedule instrucyions better,
 * and it avoids branches.  It seems like this does not help
 * on x86 machines (and can even make things worse).
 */


void print_alt_mul_body(FILE *f, long n1, long k, long fn, 
                        long half_flag, long short_flag)
{

  long n, i, chop, r;
  unsigned long mask, mask2;

  if (half_flag)
    n = n1/2;
  else
    n = n1;

  chop = n % k; /* first block */
  if (chop == 0)
    chop = k;
  r = n - k;
  mask = (1UL << k) - 1UL;

  fprintf(f, "\n\n#define NTL_ALT_");
  if (half_flag) fprintf(f, "HALF_");
  if (short_flag) fprintf(f, "SHORT_");

  fprintf(f, "BB_MUL_CODE%ld \\\n", fn);
  
  if (fn > 0) /* Mul1/AddMul1 */
    {
      fprintf(f, "   long i;\\\n");
      fprintf(f, "   _ntl_ulong carry = 0;\\\n");
    }
  fprintf(f, "   _ntl_ulong A[%ld];\\\n", 1L << k);
  fprintf(f, "   A[0] = 0;\\\n");

  fprintf(f, "   A[1] = a;\\\n");

  for (i = 2; i < (1L << k); i++)
    {
      if (i % 2 == 0)
        fprintf(f, "   A[%ld] = A[%ld] << 1;\\\n", i, i / 2);
      else
        fprintf(f, "   A[%ld] = A[%ld] ^ A[1];\\\n", i, i - 1);
    }

  if (fn > 0)
    {
      fprintf(f, "   for (i = 0; i < sb; i++) {\\\n");
      fprintf(f, "      const _ntl_ulong b = bp[i];\\\n");
    }

  for (i = k; i < n - chop; i += k) {
    if (fn > 0) fprintf(f, "   ");
    fprintf(f, "   const _ntl_ulong t%ld = A[(b >> %ld) & %lu]; \\\n", i, i, mask);

  }
  if (fn > 0) fprintf(f, "   ");
  fprintf(f, "   const _ntl_ulong t%ld = A[b >> %ld]; \\\n", n-chop, n-chop);

  if (fn > 0) fprintf(f, "   ");
  fprintf(f, "   const _ntl_ulong lo = A[b & %lu] \\\n", mask);

  for (i = k; i < n; i += k) {
    if (fn > 0) fprintf(f, "   ");
    fprintf(f, "      ^ (t%ld << %ld)", i, i);
    if (i == n - chop)
      fprintf(f, ";\\\n");
    else
      fprintf(f, "\\\n");
  }

  for (i = k; i < n; i += k) {
    if (fn > 0) fprintf(f, "   ");
    if (i == k)
      fprintf(f, "   const _ntl_ulong hi = ");
    else
      fprintf(f, "      ^ ");
    fprintf(f, "(t%ld >> %ld)", i, n1-i);
    if (i == n - chop && short_flag)
      fprintf(f, ";\\\n");
    else
      fprintf(f, "\\\n");
      
      
  }
 
  mask = 0;
  for (i = 0; i < n; i += k)
    mask |= 1UL << i;
  mask = ~mask;
  if (half_flag) mask &= (1UL << n) - 1UL;
  mask2 = mask;
  
  if (!short_flag) {
    for (i = 1; i < k; i++)
      {

        /* bit n1-i from a was not considered in blocks of
           k bits from b for index j >= i */

        if (fn > 0) fprintf(f, "   ");

	
        if (i == 1)
          fprintf(f, 
             "      ^ (((b & 0x%lxUL) >> %ld) & (-(a >> %ld)))", 
        	mask2, i, n1-1);
        else {
          fprintf(f, 
             "      ^ (((b & 0x%lxUL) >> %ld) & (-((a >> %ld) & 1UL)))",  
                  mask2, i, n1-i);
        }
        if (i == k-1) 
           fprintf(f, ";\\\n");
        else
            fprintf(f, "\\\n");


	
	mask2 = (mask2 << 1) & mask;
      }
   }

  if (fn > 0) fprintf(f, "   ");

  if (fn == 0)
    {
      fprintf(f, "   c[0] = lo; ");
      fprintf(f, "   c[1] = hi;\\\n");
    }
  else if (fn == 1)
    {
      fprintf(f, "   cp[i] = carry ^ lo; ");
      fprintf(f, "   carry = hi;\\\n");
    }
  else if (fn == 2)
    {
      fprintf(f, "   cp[i] ^= (carry ^ lo); ");
      fprintf(f, "   carry = hi;\\\n");
    }
  if (fn > 0)
    {

      fprintf(f, "   }\\\n");

      if (fn == 1 || fn == 3) 
        fprintf(f, "   cp[sb] = carry;\\\n");
      else
        fprintf(f, "   cp[sb] ^= carry;\\\n");

    }
  fprintf(f, "\n\n\n");


}




void print_alt1_mul_body(FILE *f, long n1, long k, long fn,
                    long half_flag, long short_flag)
{

  long n, i, chop, r;
  unsigned long mask, mask2;

  if (half_flag)
    n = n1/2;
  else
    n = n1;

  chop = n % k; /* first block */
  if (chop == 0)
    chop = k;
  r = n - k;
  mask = (1UL << k) - 1UL;

  fprintf(f, "\n\n#define NTL_ALT1_");
  if (half_flag) fprintf(f, "HALF_");
  if (short_flag) fprintf(f, "SHORT_");

  fprintf(f, "BB_MUL_CODE%ld \\\n", fn);

  
  if (fn > 0) /* Mul1/AddMul1 */
    {
      fprintf(f, "   long i;\\\n");
      fprintf(f, "   _ntl_ulong carry = 0, b;\\\n");
    }
  fprintf(f, "   _ntl_ulong hi, lo, t;\\\n");
  fprintf(f, "   _ntl_ulong A[%ld];\\\n", 1L << k);
  fprintf(f, "   A[0] = 0;\\\n");

  fprintf(f, "   A[1] = a;\\\n");

  for (i = 2; i < (1L << k); i++)
    {
      if (i % 2 == 0)
        fprintf(f, "   A[%ld] = A[%ld] << 1;\\\n", i, i / 2);
      else
        fprintf(f, "   A[%ld] = A[%ld] ^ A[1];\\\n", i, i - 1);
    }

  if (fn > 0)
    {
      fprintf(f, "   for (i = 0; i < sb; i++) {\\\n");
      fprintf(f, "      b = bp[i];\\\n");
      fprintf(f, "   ");
    }

  fprintf(f, "   lo = A[b & %lu]; ", mask);
  fprintf(f, "t = A[(b >> %ld) & %lu]; ", k, mask);
  fprintf(f, "hi = t >> %ld; lo ^= t << %ld;\\\n", n1-k, k);

  for (i = 2*k; i < n - chop; i += k) {
    if (fn > 0) fprintf(f, "   ");
    fprintf(f, "   t = A[(b >> %ld) & %lu]; ", i, mask);
    fprintf(f, "hi ^= t >> %ld; lo ^= t << %ld;\\\n", n1-i, i);
  }

  if (fn > 0) fprintf(f, "   ");
  fprintf(f, "   t = A[b >> %ld]; ", n-chop);
  fprintf(f, "hi ^= t >> %ld; lo ^= t << %ld;\\\n", n1-i, i);

  mask = 0;
  for (i = 0; i < n; i += k)
    mask |= 1UL << i;
  mask = ~mask;
  if (half_flag) mask &= (1UL << n) - 1UL;
  mask2 = mask;

  
  if (!short_flag) {
    for (i = 1; i < k; i++)
      {

        /* bit n1-i from a was not considered in blocks of
           k bits from b for index j >= i */

        if (fn > 0) fprintf(f, "   ");

	
        if (i == 1)
          fprintf(f, 
             "   hi ^= (((b & 0x%lxUL) >> %ld) & (-(a >> %ld)))", 
        	mask2, i, n1-1);
        else {
          fprintf(f, 
             "      ^ (((b & 0x%lxUL) >> %ld) & (-((a >> %ld) & 1UL)))",  
                  mask2, i, n1-i);
        }
        if (i == k-1) 
           fprintf(f, ";\\\n");
        else
            fprintf(f, "\\\n");


	
	mask2 = (mask2 << 1) & mask;
      }
   }

  


  if (fn > 0) fprintf(f, "   ");

  if (fn == 0)
    {
      fprintf(f, "   c[0] = lo; ");
      fprintf(f, "   c[1] = hi;\\\n");
    }
  else if (fn == 1 || fn == 3)
    {
      fprintf(f, "   cp[i] = carry ^ lo; ");
      fprintf(f, "   carry = hi;\\\n");
    }
  else if (fn == 2)
    {
      fprintf(f, "   cp[i] ^= (carry ^ lo); ");
      fprintf(f, "   carry = hi;\\\n");
    }
  if (fn > 0)
    {

      fprintf(f, "   }\\\n");

      if (fn == 1 || fn == 3) 
        fprintf(f, "   cp[sb] = carry;\\\n");
      else
        fprintf(f, "   cp[sb] ^= carry;\\\n");

    }
  fprintf(f, "\n\n\n");


}







void print_BB_mul_code(FILE *f, long n)
{
   long k;

   if (n >= 64)
      k = 4;
   else
      k = 3;


   print_mul_body(f, n, k, 0, 0, 0);
   print_mul_body(f, n, 4, 1, 0, 0);
   print_mul_body(f, n, 4, 2, 0, 0);
   print_mul_body(f, n, 4, 1, 0, 1);
   print_mul_body(f, n, 2, 0, 1, 0);



   print_alt_mul_body(f, n, k, 0, 0, 0);
   print_alt_mul_body(f, n, 4, 1, 0, 0);
   print_alt_mul_body(f, n, 4, 2, 0, 0);
   print_alt_mul_body(f, n, 4, 1, 0, 1);
   print_alt_mul_body(f, n, 2, 0, 1, 0);



   print_alt1_mul_body(f, n, k, 0, 0, 0);
   print_alt1_mul_body(f, n, 4, 1, 0, 0);
   print_alt1_mul_body(f, n, 4, 2, 0, 0);
   print_alt1_mul_body(f, n, 4, 1, 0, 1);
   print_alt1_mul_body(f, n, 2, 0, 1, 0);

   fprintf(f, "#define NTL_BB_MUL1_BITS (4)\n\n");
}





void print_BB_sqr_code(FILE *f, long n)
{
   long i, pos;

   fprintf(f, "\n\n");
   fprintf(f, "#define NTL_BB_SQR_CODE \\\n");
   fprintf(f, "lo=sqrtab[a&255];\\\n");
   pos = 16;

   for (i = 8; i < n; i += 8) {
      if (2*(i+8) <= n) {
         fprintf(f, "lo=lo|(sqrtab[(a>>%ld)&255]<<%ld);\\\n", i, pos);
         pos += 16;
      }
      else if (2*i == n) {
         fprintf(f, "hi=sqrtab[(a>>%ld)&255];\\\n", i);
         pos = 16;
      }
      else if (2*i > n) {
         fprintf(f, "hi=hi|(sqrtab[(a>>%ld)&255]<<%ld);\\\n", i, pos);
         pos += 16;
      }
      else { /* only applies if word size is not a multiple of 16 */
         fprintf(f, "_ntl_ulong t=sqrtab[(a>>%ld)&255];\\\n", i);
         fprintf(f, "lo=lo|(t<<%ld);\\\n", pos);
         fprintf(f, "hi=t>>%ld;\\\n", n-8);
         pos = 8;
      }
   }

   fprintf(f, "\n\n");
}


void print_BB_rev_code(FILE *f, long n)
{
   long i;

   fprintf(f, "\n\n");
   fprintf(f, "#define NTL_BB_REV_CODE ");

   for (i = 0; i < n; i += 8) {
      if (i != 0) fprintf(f, "\\\n|");
      fprintf(f, "(revtab[(a>>%ld)&255]<<%ld)", i, n-8-i);
   }

   fprintf(f, "\n\n");
}
   



char *yn_vec[2] = { "no", "yes" }; 



int main()
{
   long bpl, bpi, bpt, rs_arith, nbits, single_mul_ok;
   long dp, dp1, dr;
   FILE *f;
   long warnings = 0;

   unsigned long ulval;
   unsigned int uival;
   size_t tval;
   long slval;

   fprintf(stderr, "This is NTL version %s\n\n", NTL_VERSION);




   /*
    * compute bpl =  bits per long 
    */

   ulval = val_ulong(1);
   bpl = 0;

   while (ulval) {
      ulval <<= 1;
      touch_ulong(&ulval); 
      bpl++;
   }




   /*
    * compute  bpi = bits per int 
    */

   uival = val_uint(1);
   bpi = 0;

   while (uival) {
      uival <<= 1;
      touch_uint(&uival);
      bpi++;
   }



   /*
    * compute bpt = bits per size_t
    */

   tval = val_size_t(1);
   bpt = 0;

   while (tval) {
      tval <<= 1;
      touch_size_t(&tval);
      bpt++;
   }


   /*
    * check if bpl and bpi are not too small --- any standard conforming
    * platform should pass this test.
    */

   if (bpi < 16) {
      fprintf(stderr, "BAD NEWS: int type too short.\n");
      return 1;
   }

   if (bpl < 32) {
      fprintf(stderr, "BAD NEWS: long type too short.\n");
      return 1;
   }

   /*
    * check that ints are bigger than chars.
    */

   if (bpi <= CHAR_BIT) {
      fprintf(stderr, "BAD NEWS: int type must be longer than char type.\n");
      return 1;
   }



   /*
    * check that bpl is a multiple of 8.
    */

   if (bpl % 8 != 0) {
      fprintf(stderr, "BAD NEWS: word size must be multiple of 8 bits.\n");
      return 1;
   }




   /*
    * check if width of signed versions of int and long agree with that of
    * the unsigned versions, and that negative numbers are represented
    * using 2's compliment.
    *
    * The C99 standard, at least, is very precise about the possible
    * representations of unsigned and signed integer types, and so if
    * the following tests pass, we can be sure that the desired
    * properties hold.  NTL relies implicitly and crucially on 
    * these properties.
    *
    * I know of no machines for which these properties do not hold.
    */

   if (((unsigned int) val_int(INT_MIN)) != val_uint(1U << (bpi-1))) {
      fprintf(stderr, "BAD NEWS: machine must be 2's compliment.\n");
      return 1;
   }

   if (((unsigned int) val_int(INT_MAX)) != val_uint((1U << (bpi-1)) - 1U)) {
      fprintf(stderr, "BAD NEWS: machine must be 2's compliment.\n");
      return 1;
   }

   if (((unsigned long) val_long(LONG_MIN)) != val_ulong(1UL << (bpl-1))) {
      fprintf(stderr, "BAD NEWS: machine must be 2's compliment.\n");
      return 1;
   }

   if (((unsigned long) val_long(LONG_MAX)) != val_ulong((1UL<<(bpl-1))-1UL)) {
      fprintf(stderr, "BAD NEWS: machine must be 2's compliment.\n");
      return 1;
   }



   /*
    * check that floating point to integer conversions truncates toward zero
    * --- any standard conforming platform should pass this test.
    */

   if (((long) val_double(1.75)) != 1L) {
      fprintf(stderr, 
         "BAD NEWS: machine must truncate floating point toward zero.\n");
      return 1;
   }

   if (((long) val_double(-1.75)) != -1L) {
      fprintf(stderr, 
         "BAD NEWS: machine must truncate floating point toward zero.\n");
      return 1;
   }



   /*
    * Test if right shift is arithemtic or not.  According to the
    * standards, the result of right-shifting a negative number is
    * "implementation defined", which almost surely means the right shift
    * is *always* arithmetic or *always* logical.  However, this cannot
    * be guaranteed, and so this test is *not* 100% portable --- but I
    * know of no machine for which this test does not correctly
    * predict the general behavior.  One should set the NTL_CLEAN_INT
    * flag if one wants to avoid such machine dependencies.
    */

   slval = val_long(-1);
   if ((slval >> 1) == slval)
      rs_arith = 1;
   else
      rs_arith = 0;



   /*
    * Next, we check some properties of floating point arithmetic.
    * An implementation should conform to the IEEE floating
    * point standard --- essentially all modern platforms do,
    * except for a few very old Cray's.  There is no easy way
    * to check this, so we simply make a few simple sanity checks,
    * calculate the precision, and if the platform performs
    * double precision arithemtic in extended double precision registers.
    * The last property is one that the IEE standard allows, and which
    * some important platforms (like x86) have --- this is quite
    * unfortunate, as it really makes many of the other properties
    * of the IEEE standard unusable.
    */

   /*
    * First, we simply check that we are using a machine with radix 2.
    */

   if (FLT_RADIX != 2) {
      fprintf(stderr, "BAD NEWS: machine must use IEEE floating point.\n");
      return 1;
   }

   /*
    * Next, we calculate the precision of "in memory" doubles,
    * and check that it is at least 53.
    */

   dp = DoublePrecision();

   if (dp < 53) {
      fprintf(stderr, "BAD NEWS: machine must use IEEE floating point (*).\n");
      return 1;
   }

   /*
    * Next, we check that the *range* of doubles is sufficiently large.
    * Specifically, we require that DBL_MAX > 2^{7*max(bpl, dp)} 
    * and 1/DBL_MIN > 2^{7*max(bpl, dp)}.  
    * On IEEE floating point compliant machines, this
    * will hold, and the following test will pass, if bpl is at most 128, which
    * should be true for the foreseeable future.
    */

   if (log(DBL_MAX)/log(2.0) < 7.01*bpl ||  log(DBL_MAX)/log(2.0) < 7.01*dp ||
      -log(DBL_MIN)/log(2.0) < 7.01*bpl || -log(DBL_MIN)/log(2.0) < 7.01*dp) {
      fprintf(stderr, "BAD NEWS: range of doubles too small.\n");
      return 1;
   }

   

   /*
    * Next, we check if the machine has wider "in-register" doubles or not.
    * This test almost always yields the correct result --- if not,
    * you will have to set the NTL_EXT_DOUBLE in "mach_desc.h"
    * by hand.
    * 
    * The test effectively proves that in-register doubles are wide
    * if dp1 > dp || dr.
    */


   dp1 = DoublePrecision1();
   dr = DoubleRounding(dp);


   /*
    * Set nbits --- the default radix size for NTL's "built in"
    * long integer arithmetic.
    *
    *  Given the minimum size of blp and dp, the smallest possible
    *  value of nbits is 30.
    */


   if (bpl-2 < dp-3)
      nbits = bpl-2;
   else
      nbits = dp-3;

   if (nbits % 2 != 0) nbits--;


   /* 
    * We next test if the NTL_SINGLE_MUL option is valid.  This test is
    * inherently DIRTY (i.e., the behavior of the test itself is not well
    * defined according to the standard), but in practice should not cause any
    * bad behavior, especially if the NTL_SINGLE_MUL option is never used.
    * This option is anyway considered fairly obsolete, and if desired, one may
    * change the following "if 1" to "if 0" and avoid performing this test
    * altogether.
    */

#if 1
   single_mul_ok = RepTest();
#else
   single_mul_ok = 0;
#endif


   /*
    * That's it!  All tests have passed.
    */

   fprintf(stderr, "GOOD NEWS: compatible machine.\n");
   fprintf(stderr, "summary of machine characteristics:\n");
   fprintf(stderr, "bits per long = %ld\n", bpl);
   fprintf(stderr, "bits per int = %ld\n", bpi);
   fprintf(stderr, "bits per size_t = %ld\n", bpt);
   fprintf(stderr, "arith right shift = %s\n", yn_vec[rs_arith]);
   fprintf(stderr, "double precision = %ld\n", dp);
   fprintf(stderr, "NBITS (maximum) = %ld\n", nbits);
   fprintf(stderr, "single mul ok = %s\n", yn_vec[single_mul_ok != 0]);
   fprintf(stderr, "register double precision = %ld\n", dp1);
   fprintf(stderr, "double rounding detected = %s\n", yn_vec[dr]);  

   if (((dp1 > dp) || dr) && AutoFix)
      fprintf(stderr, "-- auto x86 fix\n");

   if (dp != 53) {
      warnings = 1;

      fprintf(stderr, "\n\nWARNING:\n\n"); 
      fprintf(stderr, "Nonstandard floating point precision.\n");
      fprintf(stderr, "IEEE standard is 53 bits.\n"); 
   }

#if (defined(__sparc__) && !defined(__sparc_v8__) && \
    !defined(__sparcv8) && !defined(__sparc_v9__) && !defined(__sparcv9))

   warnings = 1;

   fprintf(stderr, "\n\nWARNING:\n\n");
   fprintf(stderr, "If this Sparc is a Sparc-10 or later (so it has\n");
   fprintf(stderr, "a hardware integer multiply instruction) you\n");
   fprintf(stderr, "should specify the -mv8 option in the makefile\n"); 
   fprintf(stderr, "to obtain more efficient code.\n");

#endif

   if (((dp1 > dp) || dr) && !AutoFix) {
      warnings = 1;
      fprintf(stderr, "\n\nWARNING:\n\n");
      fprintf(stderr, "This platform has extended double precision registers.\n");
      fprintf(stderr, "While that may sound like a good thing, it actually is not.\n");
      fprintf(stderr, "If this is a Pentium or other x86 and your compiler\n");
      fprintf(stderr, "is g++ or supports GNU 'asm' constructs, it is recommended\n");
      fprintf(stderr, "to compile NTL with the NTL_X86_FIX flag to get true IEEE floating point.\n");
      fprintf(stderr, "Set this flag by editing the file config.h.\n");
      fprintf(stderr, "The code should still work even if you don't set\n");
      fprintf(stderr, "this flag.  See quad_float.txt for details.\n\n");
   }


#if 0

   /* better not to be interactive */

   if (warnings) {
      int c;

      fprintf(stderr, "Do you want to continue anyway[y/n]? ");
      c = getchar();
      if (c == 'n' || c == 'N') {
         fprintf(stderr, "Make the necessary changes to the makefile and/or config.h,\n");
         fprintf(stderr, "then type 'make clobber' and then 'make'.\n\n\n");
         return 1;
      }
   }

#endif

   f = fopen("mach_desc.h", "w");
   if (!f) {
      fprintf(stderr, "can't open mach_desc.h for writing\n");
      return 1;
   }

   fprintf(f, "#ifndef NTL_mach_desc__H\n");
   fprintf(f, "#define NTL_mach_desc__H\n\n\n");
   fprintf(f, "#define NTL_BITS_PER_LONG (%ld)\n", bpl);
   fprintf(f, "#define NTL_MAX_LONG (%ldL)\n", ((long) ((1UL<<(bpl-1))-1UL)));
   fprintf(f, "#define NTL_MAX_INT (%ld)\n", ((long) ((1UL<<(bpi-1))-1UL)));
   fprintf(f, "#define NTL_BITS_PER_INT (%ld)\n", bpi);
   fprintf(f, "#define NTL_BITS_PER_SIZE_T (%ld)\n", bpt);
   fprintf(f, "#define NTL_ARITH_RIGHT_SHIFT (%ld)\n", rs_arith);
   fprintf(f, "#define NTL_NBITS_MAX (%ld)\n", nbits);
   fprintf(f, "#define NTL_DOUBLE_PRECISION (%ld)\n", dp);
   fprintf(f, "#define NTL_FDOUBLE_PRECISION ");
   print2k(f, dp-1, bpl);
   fprintf(f, "\n");
   fprintf(f, "#define NTL_QUAD_FLOAT_SPLIT (");
   print2k(f, dp - (dp/2), bpl);
   fprintf(f, "+1.0)\n");
   fprintf(f, "#define NTL_EXT_DOUBLE (%d)\n", ((dp1 > dp) || dr));
   fprintf(f, "#define NTL_SINGLE_MUL_OK (%d)\n", single_mul_ok != 0);
   fprintf(f, "#define NTL_DOUBLES_LOW_HIGH (%d)\n\n\n", single_mul_ok < 0);
   print_BB_mul_code(f, bpl);
   print_BB_sqr_code(f, bpl);
   print_BB_rev_code(f, bpl);

   fprintf(f, "#define NTL_MIN_LONG (-NTL_MAX_LONG - 1L)\n");
   fprintf(f, "#define NTL_MIN_INT  (-NTL_MAX_INT - 1)\n");

   fprintf(f, "#endif\n\n");

   fclose(f);

   fprintf(stderr, "\n\n");
   
   return 0;
}
