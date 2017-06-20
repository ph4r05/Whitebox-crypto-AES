

#include <NTL/FFT.h>

#include <NTL/new.h>

NTL_START_IMPL


FFTPrimeInfo *FFTTables = 0;
Vec<FFTPrimeInfo> FFTTables_store;

long *FFTPrime = 0;
Vec<long> FFTPrime_store;

double *FFTPrimeInv = 0;
Vec<double> FFTPrimeInv_store;

// We separate the pointer from the Vec, to ensure
// portability: global initialization of C++ objects
// can be problematic.



long NumFFTPrimes = 0;




static
long IsFFTPrime(long n, long& w)
{
   long  m, x, y, z;
   long j, k;

   if (n % 3 == 0) return 0;

   if (n % 5 == 0) return 0;

   if (n % 7 == 0) return 0;
   
   m = n - 1;
   k = 0;
   while ((m & 1) == 0) {
      m = m >> 1;
      k++;
   }

   for (;;) {
      x = RandomBnd(n);

      if (x == 0) continue;
      z = PowerMod(x, m, n);
      if (z == 1) continue;

      x = z;
      j = 0;
      do {
         y = z;
         z = MulMod(y, y, n);
         j++;
      } while (j != k && z != 1);

      if (z != 1 || y !=  n-1) return 0;

      if (j == k) 
         break;
   }

   /* x^{2^k} = 1 mod n, x^{2^{k-1}} = -1 mod n */

   long TrialBound;

   TrialBound = m >> k;
   if (TrialBound > 0) {
      if (!ProbPrime(n, 5)) return 0;
   
      /* we have to do trial division by special numbers */
   
      TrialBound = SqrRoot(TrialBound);
   
      long a, b;
   
      for (a = 1; a <= TrialBound; a++) {
         b = (a << k) + 1;
         if (n % b == 0) return 0; 
      }
   }

   /* n is an FFT prime */

   for (j = NTL_FFTMaxRoot; j < k; j++)
      x = MulMod(x, x, n);

   w = x;
   return 1;
}


static
void NextFFTPrime(long& q, long& w)
{
   static long m = NTL_FFTMaxRootBnd + 1;
   static long k = 0;

   long t, cand;

   for (;;) {
      if (k == 0) {
         m--;
         if (m < 5) Error("ran out of FFT primes");
         k = 1L << (NTL_SP_NBITS-m-2);
      }

      k--;

      cand = (1L << (NTL_SP_NBITS-1)) + (k << (m+1)) + (1L << m) + 1;

      if (!IsFFTPrime(cand, t)) continue;
      q = cand;
      w = t;
      return;
   }
}


long CalcMaxRoot(long p)
{
   p = p-1;
   long k = 0;
   while ((p & 1) == 0) {
      p = p >> 1;
      k++;
   }

   if (k > NTL_FFTMaxRoot)
      return NTL_FFTMaxRoot;
   else
      return k; 
}



void UseFFTPrime(long index)
{
   long numprimes = FFTTables_store.length();

   if (index < 0 || index > numprimes)
      Error("invalid FFT prime index");

   if (index < numprimes) return;

   // index == numprimes

   long q, w;

   NextFFTPrime(q, w);

   double qinv = 1/((double) q);

   long mr = CalcMaxRoot(q);

   FFTTables_store.SetLength(numprimes+1);
   FFTTables = FFTTables_store.elts();

   FFTPrimeInfo& info = FFTTables[numprimes];

   info.q = q;
   info.qinv = qinv;

   info.RootTable.SetLength(mr+1);
   info.RootInvTable.SetLength(mr+1);
   info.TwoInvTable.SetLength(mr+1);
   info.TwoInvPreconTable.SetLength(mr+1);

   long *rt = &info.RootTable[0];
   long *rit = &info.RootInvTable[0];
   long *tit = &info.TwoInvTable[0];
   mulmod_precon_t *tipt = &info.TwoInvPreconTable[0];

   long j;
   long t;

   rt[mr] = w;
   for (j = mr-1; j >= 0; j--)
      rt[j] = MulMod(rt[j+1], rt[j+1], q);

   rit[mr] = InvMod(w, q);
   for (j = mr-1; j >= 0; j--)
      rit[j] = MulMod(rit[j+1], rit[j+1], q);

   t = InvMod(2, q);
   tit[0] = 1;
   for (j = 1; j <= mr; j++)
      tit[j] = MulMod(tit[j-1], t, q);

   for (j = 0; j <= mr; j++)
      tipt[j] = PrepMulModPrecon(tit[j], q, qinv);


   // initialize data structures for the legacy inteface

   NumFFTPrimes = FFTTables_store.length();
   
   FFTPrime_store.SetLength(NumFFTPrimes);
   FFTPrime = FFTPrime_store.elts();
   FFTPrime[NumFFTPrimes-1] = q;

   FFTPrimeInv_store.SetLength(NumFFTPrimes);
   FFTPrimeInv = FFTPrimeInv_store.elts();
   FFTPrimeInv[NumFFTPrimes-1] = qinv;
}



/*
 * Our FFT is based on the routine in Cormen, Leiserson, Rivest, and Stein.
 * For very large inputs, it should be relatively cache friendly.
 * The inner loop has been unrolled and pipelined, to exploit any
 * low-level parallelism in the machine.
 * 
 * This version now allows input to alias output.
 */




static
long RevInc(long a, long k)
{
   long j, m;

   j = k; 
   m = 1L << (k-1);

   while (j && (m & a)) {
      a ^= m;
      m >>= 1;
      j--;
   }
   if (j) a ^= m;
   return a;
}



static Vec<long> brc_mem[NTL_FFTMaxRoot+1];

static
void BitReverseCopy(long *A, const long *a, long k)
{
   long n = 1L << k;
   long* rev;
   long i, j;

   rev = brc_mem[k].elts();
   if (!rev) {
      brc_mem[k].SetLength(n);
      rev = brc_mem[k].elts();
      for (i = 0, j = 0; i < n; i++, j = RevInc(j, k))
         rev[i] = j;
   }

   for (i = 0; i < n; i++)
      A[rev[i]] = a[i];
}

static
void BitReverseCopy(unsigned long *A, const long *a, long k)
{
   long n = 1L << k;
   long* rev;
   long i, j;

   rev = brc_mem[k].elts();
   if (!rev) {
      brc_mem[k].SetLength(n);
      rev = brc_mem[k].elts();
      for (i = 0, j = 0; i < n; i++, j = RevInc(j, k))
         rev[i] = j;
   }

   for (i = 0; i < n; i++)
      A[rev[i]] = a[i];
}



void FFT(long* A, const long* a, long k, long q, const long* root)
// performs a 2^k-point convolution modulo q

{
   if (k <= 1) {
      if (k == 0) {
	 A[0] = a[0];
	 return;
      }
      if (k == 1) {
	 long a0 = AddMod(a[0], a[1], q);
	 long a1 = SubMod(a[0], a[1], q);
         A[0] = a0;
         A[1] = a1;
	 return;
      }
   }

   // assume k > 1

   

   static Vec<long> wtab_store;
   static Vec<mulmod_precon_t> wqinvtab_store;
   static Vec<long> AA_store;

   wtab_store.SetLength(1L << (k-2));
   wqinvtab_store.SetLength(1L << (k-2));
   AA_store.SetLength(1L << k);

   long * NTL_RESTRICT wtab = wtab_store.elts();
   mulmod_precon_t * NTL_RESTRICT wqinvtab = wqinvtab_store.elts();
   long *AA = AA_store.elts();

   double qinv = 1/((double) q);

   wtab[0] = 1;
   wqinvtab[0] = PrepMulModPrecon(1, q, qinv);


   BitReverseCopy(AA, a, k);

   long n = 1L << k;

   long s, m, m_half, m_fourth, i, j, t, u, t1, u1, tt, tt1;

   long w;
   mulmod_precon_t wqinv;

   // s = 1

   for (i = 0; i < n; i += 2) {
      t = AA[i + 1];
      u = AA[i];
      AA[i] = AddMod(u, t, q);
      AA[i+1] = SubMod(u, t, q);
   }

   
  
   for (s = 2; s < k; s++) {
      m = 1L << s;
      m_half = 1L << (s-1);
      m_fourth = 1L << (s-2);

      w = root[s];
      wqinv = PrepMulModPrecon(w, q, qinv);

      // prepare wtab...

      if (s == 2) {
         wtab[1] = MulModPrecon(wtab[0], w, q, wqinv);
         wqinvtab[1] = PrepMulModPrecon(wtab[1], q, qinv);
      }
      else {
         // some software pipelining

         i = m_half-1; j = m_fourth-1;
         wtab[i-1] = wtab[j];
         wqinvtab[i-1] = wqinvtab[j];
         wtab[i] = MulModPrecon(wtab[i-1], w, q, wqinv);

         i -= 2; j --;

         for (; i >= 0; i -= 2, j --) {
            long wp2 = wtab[i+2];
            long wm1 = wtab[j];
            wqinvtab[i+2] = PrepMulModPrecon(wp2, q, qinv);
            wtab[i-1] = wm1;
            wqinvtab[i-1] = wqinvtab[j];
            wtab[i] = MulModPrecon(wm1, w, q, wqinv);
         }

         wqinvtab[1] = PrepMulModPrecon(wtab[1], q, qinv);
      }

      for (i = 0; i < n; i+= m) {

         long * NTL_RESTRICT AA0 = &AA[i];
         long * NTL_RESTRICT AA1 = &AA[i + m_half];
          
         t = AA1[0];
         u = AA0[0];
         t1 = MulModPrecon(AA1[1], w, q, wqinv);
         u1 = AA0[1];

         for (j = 0; j < m_half-2; j += 2) {
            long a02 = AA0[j+2];
            long a03 = AA0[j+3];
            long a12 = AA1[j+2];
            long a13 = AA1[j+3];
            long w2 = wtab[j+2];
            long w3 = wtab[j+3];
            mulmod_precon_t wqi2 = wqinvtab[j+2];
            mulmod_precon_t wqi3 = wqinvtab[j+3];

            tt = MulModPrecon(a12, w2, q, wqi2);
            long b00 = AddMod(u, t, q);
            long b10 = SubMod(u, t, q);
            t = tt;
            u = a02;

            tt1 = MulModPrecon(a13, w3, q, wqi3);
            long b01 = AddMod(u1, t1, q);
            long b11 = SubMod(u1, t1, q);
            t1 = tt1;
            u1 = a03;

            AA0[j] = b00;
            AA1[j] = b10;
            AA0[j+1] = b01;
            AA1[j+1] = b11;
         }


         AA0[j] = AddMod(u, t, q);
         AA1[j] = SubMod(u, t, q);
         AA0[j + 1] = AddMod(u1, t1, q);
         AA1[j + 1] = SubMod(u1, t1, q);
      }
   }


   // s == k...special case

   m = 1L << s;
   m_half = 1L << (s-1);
   m_fourth = 1L << (s-2);


   w = root[s];
   wqinv = PrepMulModPrecon(w, q, qinv);

   // j = 0, 1

   t = AA[m_half];
   u = AA[0];
   t1 = MulModPrecon(AA[1+ m_half], w, q, wqinv);
   u1 = AA[1];

   A[0] = AddMod(u, t, q);
   A[m_half] = SubMod(u, t, q);
   A[1] = AddMod(u1, t1, q);
   A[1 + m_half] = SubMod(u1, t1, q);

   for (j = 2; j < m_half; j += 2) {
      t = MulModPrecon(AA[j + m_half], wtab[j >> 1], q, wqinvtab[j >> 1]);
      u = AA[j];
      t1 = MulModPrecon(AA[j + 1+ m_half], wtab[j >> 1], q, 
                        wqinvtab[j >> 1]);
      t1 = MulModPrecon(t1, w, q, wqinv);
      u1 = AA[j + 1];

      A[j] = AddMod(u, t, q);
      A[j + m_half] = SubMod(u, t, q);
      A[j + 1] = AddMod(u1, t1, q);
      A[j + 1 + m_half] = SubMod(u1, t1, q);
     
   }
}



#if (!defined(NTL_FFT_LAZYMUL) || defined(NTL_SINGLE_MUL) || \
     (!defined(NTL_SPMM_ULL) && !defined(NTL_SPMM_ASM)))

// FFT with precomputed tables 

#define NTL_PIPELINE (1)

static
void PrecompFFTMultipliers(long k, long q, const long *root, FFTMultipliers& tab)
{
   if (k < 1) Error("PrecompFFTMultipliers: bad input");

   if (k <= tab.MaxK) return;

   tab.wtab_precomp.SetLength(k+1);
   tab.wqinvtab_precomp.SetLength(k+1);

   double qinv = 1/((double) q);

   if (tab.MaxK == -1) {
      tab.wtab_precomp[1].SetLength(1);
      tab.wqinvtab_precomp[1].SetLength(1);
      tab.wtab_precomp[1][0] = 1;
      tab.wqinvtab_precomp[1][0] = PrepMulModPrecon(1, q, qinv);
      tab.MaxK = 1;
   }

   for (long s = tab.MaxK+1; s <= k; s++) {
      tab.wtab_precomp[s].SetLength(1L << (s-1));
      tab.wqinvtab_precomp[s].SetLength(1L << (s-1));

      long m = 1L << s;
      long m_half = 1L << (s-1);
      long m_fourth = 1L << (s-2);

      long *wtab_last = tab.wtab_precomp[s-1].elts();
      mulmod_precon_t *wqinvtab_last = tab.wqinvtab_precomp[s-1].elts();

      long *wtab = tab.wtab_precomp[s].elts();
      mulmod_precon_t *wqinvtab = tab.wqinvtab_precomp[s].elts();

      for (long i = 0; i < m_fourth; i++) {
         wtab[i] = wtab_last[i];
         wqinvtab[i] = wqinvtab_last[i];
      } 

      long w = root[s];
      mulmod_precon_t wqinv = PrepMulModPrecon(w, q, qinv);

      // prepare wtab...

      if (s == 2) {
         wtab[1] = MulModPrecon(wtab[0], w, q, wqinv);
         wqinvtab[1] = PrepMulModPrecon(wtab[1], q, qinv);
      }
      else {
         // some software pipelining
         long i, j;

         i = m_half-1; j = m_fourth-1;
         wtab[i-1] = wtab[j];
         wqinvtab[i-1] = wqinvtab[j];
         wtab[i] = MulModPrecon(wtab[i-1], w, q, wqinv);

         i -= 2; j --;

         for (; i >= 0; i -= 2, j --) {
            long wp2 = wtab[i+2];
            long wm1 = wtab[j];
            wqinvtab[i+2] = PrepMulModPrecon(wp2, q, qinv);
            wtab[i-1] = wm1;
            wqinvtab[i-1] = wqinvtab[j];
            wtab[i] = MulModPrecon(wm1, w, q, wqinv);
         }

         wqinvtab[1] = PrepMulModPrecon(wtab[1], q, qinv);
      }
   }

   tab.MaxK = k;
}



void FFT(long* A, const long* a, long k, long q, const long* root, FFTMultipliers& tab)
// performs a 2^k-point convolution modulo q

{
   if (k <= 1) {
      if (k == 0) {
	 A[0] = a[0];
	 return;
      }
      if (k == 1) {
	 long a0 = AddMod(a[0], a[1], q);
	 long a1 = SubMod(a[0], a[1], q);
         A[0] = a0;
         A[1] = a1;
	 return;
      }
   }



   // assume k > 1

   if (k > tab.MaxK) PrecompFFTMultipliers(k, q, root, tab);

   static Vec<long> AA_store;
   AA_store.SetLength(1L << k);
   long *AA = AA_store.elts();

   BitReverseCopy(AA, a, k);

   long n = 1L << k;

   long s, m, m_half, m_fourth, i, j, t, u, t1, u1, tt, tt1;

   // s = 1

   for (i = 0; i < n; i += 2) {
      t = AA[i + 1];
      u = AA[i];
      AA[i] = AddMod(u, t, q);
      AA[i+1] = SubMod(u, t, q);
   }
   
  
   for (s = 2; s < k; s++) {
      m = 1L << s;
      m_half = 1L << (s-1);
      m_fourth = 1L << (s-2);

      const long* wtab = tab.wtab_precomp[s].elts();
      const mulmod_precon_t *wqinvtab = tab.wqinvtab_precomp[s].elts();

      for (i = 0; i < n; i+= m) {

         long *AA0 = &AA[i];
         long *AA1 = &AA[i + m_half];

#if (NTL_PIPELINE)

// pipelining: seems to be faster
          
         t = AA1[0];
         u = AA0[0];
         t1 = MulModPrecon(AA1[1], wtab[1], q, wqinvtab[1]);
         u1 = AA0[1];

         for (j = 0; j < m_half-2; j += 2) {
            long a02 = AA0[j+2];
            long a03 = AA0[j+3];
            long a12 = AA1[j+2];
            long a13 = AA1[j+3];
            long w2 = wtab[j+2];
            long w3 = wtab[j+3];
            mulmod_precon_t wqi2 = wqinvtab[j+2];
            mulmod_precon_t wqi3 = wqinvtab[j+3];

            tt = MulModPrecon(a12, w2, q, wqi2);
            long b00 = AddMod(u, t, q);
            long b10 = SubMod(u, t, q);

            tt1 = MulModPrecon(a13, w3, q, wqi3);
            long b01 = AddMod(u1, t1, q);
            long b11 = SubMod(u1, t1, q);

            AA0[j] = b00;
            AA1[j] = b10;
            AA0[j+1] = b01;
            AA1[j+1] = b11;


            t = tt;
            u = a02;
            t1 = tt1;
            u1 = a03;
         }


         AA0[j] = AddMod(u, t, q);
         AA1[j] = SubMod(u, t, q);
         AA0[j + 1] = AddMod(u1, t1, q);
         AA1[j + 1] = SubMod(u1, t1, q);
      }
#else
         for (j = 0; j < m_half; j += 2) {
            const long a00 = AA0[j];
            const long a01 = AA0[j+1];
            const long a10 = AA1[j];
            const long a11 = AA1[j+1];

            const long w0 = wtab[j];
            const long w1 = wtab[j+1];
            const mulmod_precon_t wqi0 = wqinvtab[j];
            const mulmod_precon_t wqi1 = wqinvtab[j+1];

            const long tt = MulModPrecon(a10, w0, q, wqi0);
            const long uu = a00;
            const long b00 = AddMod(uu, tt, q); 
            const long b10 = SubMod(uu, tt, q);

            const long tt1 = MulModPrecon(a11, w1, q, wqi1);
            const long uu1 = a01;
            const long b01 = AddMod(uu1, tt1, q); 
            const long b11 = SubMod(uu1, tt1, q);

            AA0[j] = b00;
            AA0[j+1] = b01;
            AA1[j] = b10;
            AA1[j+1] = b11;
         }
      }
#endif
   }


   // s == k, special case
   {
      m = 1L << s;
      m_half = 1L << (s-1);
      m_fourth = 1L << (s-2);

      const long* wtab = tab.wtab_precomp[s].elts();
      const mulmod_precon_t *wqinvtab = tab.wqinvtab_precomp[s].elts();

      for (i = 0; i < n; i+= m) {

         long *AA0 = &AA[i];
         long *AA1 = &AA[i + m_half];
         long *A0 = &A[i];
         long *A1 = &A[i + m_half];

#if (NTL_PIPELINE)

// pipelining: seems to be faster
          
         t = AA1[0];
         u = AA0[0];
         t1 = MulModPrecon(AA1[1], wtab[1], q, wqinvtab[1]);
         u1 = AA0[1];

         for (j = 0; j < m_half-2; j += 2) {
            long a02 = AA0[j+2];
            long a03 = AA0[j+3];
            long a12 = AA1[j+2];
            long a13 = AA1[j+3];
            long w2 = wtab[j+2];
            long w3 = wtab[j+3];
            mulmod_precon_t wqi2 = wqinvtab[j+2];
            mulmod_precon_t wqi3 = wqinvtab[j+3];

            tt = MulModPrecon(a12, w2, q, wqi2);
            long b00 = AddMod(u, t, q);
            long b10 = SubMod(u, t, q);

            tt1 = MulModPrecon(a13, w3, q, wqi3);
            long b01 = AddMod(u1, t1, q);
            long b11 = SubMod(u1, t1, q);

            A0[j] = b00;
            A1[j] = b10;
            A0[j+1] = b01;
            A1[j+1] = b11;


            t = tt;
            u = a02;
            t1 = tt1;
            u1 = a03;
         }


         A0[j] = AddMod(u, t, q);
         A1[j] = SubMod(u, t, q);
         A0[j + 1] = AddMod(u1, t1, q);
         A1[j + 1] = SubMod(u1, t1, q);
      }
#else
         for (j = 0; j < m_half; j += 2) {
            const long a00 = AA0[j];
            const long a01 = AA0[j+1];
            const long a10 = AA1[j];
            const long a11 = AA1[j+1];

            const long w0 = wtab[j];
            const long w1 = wtab[j+1];
            const mulmod_precon_t wqi0 = wqinvtab[j];
            const mulmod_precon_t wqi1 = wqinvtab[j+1];

            const long tt = MulModPrecon(a10, w0, q, wqi0);
            const long uu = a00;
            const long b00 = AddMod(uu, tt, q); 
            const long b10 = SubMod(uu, tt, q);

            const long tt1 = MulModPrecon(a11, w1, q, wqi1);
            const long uu1 = a01;
            const long b01 = AddMod(uu1, tt1, q); 
            const long b11 = SubMod(uu1, tt1, q);

            A0[j] = b00;
            A0[j+1] = b01;
            A1[j] = b10;
            A1[j+1] = b11;
         }
      }
#endif
   }

}


#else

// FFT with precomputed tables and David Harvey's lazy multiplication
// strategy.


static inline 
unsigned long LazyPrepMulModPrecon(long b, long n, double ninv)
{
   unsigned long q, r;

   q = (long) ( (((double) b) * NTL_SP_FBOUND) * ninv ); 
   r = (((unsigned long) b) << NTL_SP_NBITS ) - q * ((unsigned long) n);

   if (r >> (NTL_BITS_PER_LONG-1)) {
      q--;
      r += n;
   }
   else if (((long) r) >= n) {
      q++;
      r -=n;
   }

   unsigned long res = q << (NTL_BITS_PER_LONG - NTL_SP_NBITS);
   long qq, rr;

   rr = MulDivRem(qq, (long) r, 4, n, 4*ninv);

   res = res + (qq << (NTL_BITS_PER_LONG - NTL_SP_NBITS-2));

   return res;
}

static inline 
unsigned long LazyMulModPrecon(unsigned long a, unsigned long b, 
                               unsigned long n, unsigned long bninv)
{
   unsigned long q = MulHiUL(a, bninv);
   unsigned long res = a*b - q*n;
   return res;
}


static inline 
unsigned long LazyReduce(unsigned long a, unsigned long q)
{
  unsigned long res;
#if (NTL_ARITH_RIGHT_SHIFT && defined(NTL_AVOID_BRANCHING) && !defined(NTL_CLEAN_INT))
  res = a - q;
  res  += (((long) res) >> (NTL_BITS_PER_LONG-1)) & q; 
#elif (defined(NTL_AVOID_BRANCHING))
  res = a - q;
  res  += (-(res >> (NTL_BITS_PER_LONG-1))) & q; 
#else
  if (a >= q)
    res = a - q;
  else
    res = a;
#endif

  return res;
}


static
void LazyPrecompFFTMultipliers(long k, long q, const long *root, FFTMultipliers& tab)
{
   if (k < 1) Error("LazyPrecompFFTMultipliers: bad input");

   if (k <= tab.MaxK) return;

   tab.wtab_precomp.SetLength(k+1);
   tab.wqinvtab_precomp.SetLength(k+1);

   double qinv = 1/((double) q);

   if (tab.MaxK == -1) {
      tab.wtab_precomp[1].SetLength(1);
      tab.wqinvtab_precomp[1].SetLength(1);
      tab.wtab_precomp[1][0] = 1;
      tab.wqinvtab_precomp[1][0] = LazyPrepMulModPrecon(1, q, qinv);
      tab.MaxK = 1;
   }

   for (long s = tab.MaxK+1; s <= k; s++) {
      tab.wtab_precomp[s].SetLength(1L << (s-1));
      tab.wqinvtab_precomp[s].SetLength(1L << (s-1));

      long m = 1L << s;
      long m_half = 1L << (s-1);
      long m_fourth = 1L << (s-2);

      long *wtab_last = tab.wtab_precomp[s-1].elts();
      mulmod_precon_t *wqinvtab_last = tab.wqinvtab_precomp[s-1].elts();

      long *wtab = tab.wtab_precomp[s].elts();
      mulmod_precon_t *wqinvtab = tab.wqinvtab_precomp[s].elts();

      for (long i = 0; i < m_fourth; i++) {
         wtab[i] = wtab_last[i];
         wqinvtab[i] = wqinvtab_last[i];
      } 

      long w = root[s];
      mulmod_precon_t wqinv = LazyPrepMulModPrecon(w, q, qinv);

      // prepare wtab...

      if (s == 2) {
         wtab[1] = LazyReduce(LazyMulModPrecon(wtab[0], w, q, wqinv), q);
         wqinvtab[1] = LazyPrepMulModPrecon(wtab[1], q, qinv);
      }
      else {
         // some software pipelining
         long i, j;

         i = m_half-1; j = m_fourth-1;
         wtab[i-1] = wtab[j];
         wqinvtab[i-1] = wqinvtab[j];
         wtab[i] = LazyReduce(LazyMulModPrecon(wtab[i-1], w, q, wqinv), q);

         i -= 2; j --;

         for (; i >= 0; i -= 2, j --) {
            long wp2 = wtab[i+2];
            long wm1 = wtab[j];
            wqinvtab[i+2] = LazyPrepMulModPrecon(wp2, q, qinv);
            wtab[i-1] = wm1;
            wqinvtab[i-1] = wqinvtab[j];
            wtab[i] = LazyReduce(LazyMulModPrecon(wm1, w, q, wqinv), q);
         }

         wqinvtab[1] = LazyPrepMulModPrecon(wtab[1], q, qinv);
      }
   }

   tab.MaxK = k;
}




void FFT(long* A, const long* a, long k, long q, const long* root, FFTMultipliers& tab)

// performs a 2^k-point convolution modulo q

{
   if (k <= 1) {
      if (k == 0) {
	 A[0] = a[0];
	 return;
      }
      if (k == 1) {
	 long a0 = AddMod(a[0], a[1], q);
	 long a1 = SubMod(a[0], a[1], q);
         A[0] = a0;
         A[1] = a1;
	 return;
      }
   }

   // assume k > 1

   if (k > tab.MaxK) LazyPrecompFFTMultipliers(k, q, root, tab);

   static Vec<unsigned long> AA_store;
   AA_store.SetLength(1L << k);
   unsigned long *AA = AA_store.elts();



   BitReverseCopy(AA, a, k);

   long n = 1L << k;


   /* we work with redundant representations, in the range [0, 4q) */



   long s, m, m_half, m_fourth, i, j; 
   unsigned long t, u, t1, u1;

   long two_q = 2 * q; 

   // s = 1

   for (i = 0; i < n; i += 2) {
      t = AA[i + 1];
      u = AA[i];
      AA[i] = u + t;
      AA[i+1] = u - t + q;
   }


   // s = 2

   {
      const long * NTL_RESTRICT wtab = tab.wtab_precomp[2].elts();
      const mulmod_precon_t * NTL_RESTRICT wqinvtab = tab.wqinvtab_precomp[2].elts();

      for (i = 0; i < n; i += 4) {

         unsigned long * NTL_RESTRICT AA0 = &AA[i];
         unsigned long * NTL_RESTRICT AA1 = &AA[i + 2];

         {
            const long w1 = wtab[0];
            const mulmod_precon_t wqi1 = wqinvtab[0];
            const unsigned long a11 = AA1[0];
            const unsigned long a01 = AA0[0];

            const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
            const unsigned long uu1 = LazyReduce(a01, two_q);
            const unsigned long b01 = uu1 + tt1; 
            const unsigned long b11 = uu1 - tt1 + two_q;

            AA0[0] = b01;
            AA1[0] = b11;
         }
         {
            const long w1 = wtab[1];
            const mulmod_precon_t wqi1 = wqinvtab[1];
            const unsigned long a11 = AA1[1];
            const unsigned long a01 = AA0[1];

            const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
            const unsigned long uu1 = LazyReduce(a01, two_q);
            const unsigned long b01 = uu1 + tt1; 
            const unsigned long b11 = uu1 - tt1 + two_q;

            AA0[1] = b01;
            AA1[1] = b11;
         }
      }
   }


   //  s = 3..k

   for (s = 3; s <= k; s++) {
      m = 1L << s;
      m_half = 1L << (s-1);
      m_fourth = 1L << (s-2);

      const long* NTL_RESTRICT wtab = tab.wtab_precomp[s].elts();
      const mulmod_precon_t * NTL_RESTRICT wqinvtab = tab.wqinvtab_precomp[s].elts();

      for (i = 0; i < n; i += m) {

         unsigned long * NTL_RESTRICT AA0 = &AA[i];
         unsigned long * NTL_RESTRICT AA1 = &AA[i + m_half];

         for (j = 0; j < m_half; j += 4) {
            {
               const long w1 = wtab[j+0];
               const mulmod_precon_t wqi1 = wqinvtab[j+0];
               const unsigned long a11 = AA1[j+0];
               const unsigned long a01 = AA0[j+0];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce(a01, two_q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + two_q;

               AA0[j+0] = b01;
               AA1[j+0] = b11;
            }
            {
               const long w1 = wtab[j+1];
               const mulmod_precon_t wqi1 = wqinvtab[j+1];
               const unsigned long a11 = AA1[j+1];
               const unsigned long a01 = AA0[j+1];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce(a01, two_q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + two_q;

               AA0[j+1] = b01;
               AA1[j+1] = b11;
            }
            {
               const long w1 = wtab[j+2];
               const mulmod_precon_t wqi1 = wqinvtab[j+2];
               const unsigned long a11 = AA1[j+2];
               const unsigned long a01 = AA0[j+2];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce(a01, two_q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + two_q;

               AA0[j+2] = b01;
               AA1[j+2] = b11;
            }
            {
               const long w1 = wtab[j+3];
               const mulmod_precon_t wqi1 = wqinvtab[j+3];
               const unsigned long a11 = AA1[j+3];
               const unsigned long a01 = AA0[j+3];

               const unsigned long tt1 = LazyMulModPrecon(a11, w1, q, wqi1);
               const unsigned long uu1 = LazyReduce(a01, two_q);
               const unsigned long b01 = uu1 + tt1; 
               const unsigned long b11 = uu1 - tt1 + two_q;

               AA0[j+3] = b01;
               AA1[j+3] = b11;
            }
         }
      }
   }

   /* need to reduce redundant representations */

   for (i = 0; i < n; i++) {
      unsigned long tmp = LazyReduce(AA[i], two_q);
      A[i] = LazyReduce(tmp, q);
   }
}


#endif





NTL_END_IMPL
