
#include <NTL/GF2XVec.h>

#include <NTL/new.h>

NTL_START_IMPL


void GF2XVec::SetSize(long n, long d)
{
   if (n < 0 || d <= 0) Error("bad args to GF2XVec::SetSize()");

   if (v)
      Error("illegal GF2XVec initialization");

   len = n;
   bsize = d;

   if (n == 0) return;

   v = (GF2X*) NTL_MALLOC(n, sizeof(GF2X), 0);
   if (!v) Error("out of memory in GF2XVec::SetSize()");

   long i = 0;
   long m;
   long j;

   while (i < n) {
      m = WV_BlockConstructAlloc(v[i].xrep, d, n-i);
      for (j = 1; j < m; j++)
         WV_BlockConstructSet(v[i].xrep, v[i+j].xrep, j);
      i += m;
   }
}


void GF2XVec::kill()
{
   long n = len;

   len = 0; bsize = 0;

   if (n == 0) return;

   long i = 0;
   long m;

   while (i < n) {
      m = WV_BlockDestroy(v[i].xrep);
      i += m;
   }

   free(v);
   v = 0; 
}


GF2XVec& GF2XVec::operator=(const GF2XVec& a) 
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
   
GF2XVec::GF2XVec(const GF2XVec& a)
{
   v = 0; len = 0; bsize = 0;

   SetSize(a.len, a.bsize);

   long i;
   for (i = 0; i < a.len; i++)
      v[i] = (a.v)[i];
}

void GF2XVec::swap_impl(GF2XVec& x, GF2XVec& y)
{
   GF2X* t1;
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
