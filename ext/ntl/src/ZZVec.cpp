
#include <NTL/ZZVec.h>

#include <NTL/new.h>

NTL_START_IMPL

void ZZVec::SetSize(long n, long d)
{
   if (n < 0 || d <= 0) Error("bad args to ZZVec::SetSize()");

   if (v)
      Error("illegal ZZVec initialization");

   len = n;
   bsize = d;

   if (n == 0) return;

   v = (ZZ*) NTL_MALLOC(n, sizeof(ZZ), 0);
   if (!v) Error("out of memory in ZZVec::SetSize()");

   long i = 0;
   long m;
   long j;

   while (i < n) {
      m = ZZ_BlockConstructAlloc(v[i], d, n-i);
      for (j = 1; j < m; j++)
         ZZ_BlockConstructSet(v[i], v[i+j], j);
      i += m;
   }
}

void ZZVec::kill()
{
   long n = len;

   len = 0; bsize = 0;

   if (n == 0) return;

   long i = 0;
   long m;

   while (i < n) {
      m = ZZ_BlockDestroy(v[i]);
      i += m;
   }

   free(v);
   v = 0; 
}


ZZVec& ZZVec::operator=(const ZZVec& a) 
{
   if (this == &a)
      return *this;

   kill();
   SetSize(a.len, a.bsize);

   long i;
   for (i = 0; i < a.len; i++)
      v[i] = (a.v)[i];

  return *this;
}
   
ZZVec::ZZVec(const ZZVec& a)
{
   v = 0; len = 0; bsize = 0;

   SetSize(a.len, a.bsize);

   long i;
   for (i = 0; i < a.len; i++)
      v[i] = (a.v)[i];
}

void ZZVec::swap_impl(ZZVec& x, ZZVec& y)
{
   ZZ* t1;
   long t2;

   t1 = x.v;
   x.v = y.v;
   y.v = t1;

   t2 = x.len;
   x.len = y.len;
   y.len = t2;

   t2 = x.bsize;
   x.bsize = y.bsize;
   y.bsize = t2;
}

NTL_END_IMPL
