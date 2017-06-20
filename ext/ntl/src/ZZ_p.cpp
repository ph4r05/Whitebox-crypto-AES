

#include <NTL/ZZ_p.h>
#include <NTL/FFT.h>

#include <NTL/new.h>


NTL_START_IMPL


ZZ_pInfoT::ZZ_pInfoT(const ZZ& NewP)
{
   if (NewP <= 1) Error("ZZ_pContext: p must be > 1");

   ref_count = 1;
   p = NewP;
   size = p.size();

   ExtendedModulusSize = 2*size + 
                 (NTL_BITS_PER_LONG + NTL_ZZ_NBITS - 1)/NTL_ZZ_NBITS;

   initialized = 0;
   x = 0;
   u = 0;
   tbl = 0;
   tbl1 = 0;

   long i;
   for (i = 0; i < MAX_ZZ_p_TEMPS; i++)
      temps[i] = 0;

   temps_top = 0;
}



void ZZ_pInfoT::init()
{
   ZZ B, M, M1, M2, M3;
   long n, i;
   long q, t;

   initialized = 1;

   sqr(B, p);

   LeftShift(B, B, NTL_FFTMaxRoot+NTL_FFTFudge);

   set(M);
   n = 0;
   while (M <= B) {
      UseFFTPrime(n);
      q = FFTPrime[n];
      n++;
      mul(M, M, q);
   }

   NumPrimes = n;
   MaxRoot = CalcMaxRoot(q);


   double fn = double(n);

   if (8.0*fn*(fn+32) > NTL_FDOUBLE_PRECISION)
      Error("modulus too big");


   if (8.0*fn*(fn+32) > NTL_FDOUBLE_PRECISION/double(NTL_SP_BOUND))
      QuickCRT = 0;
   else
      QuickCRT = 1;


   if (!(x = (double *) NTL_MALLOC(n, sizeof(double), 0)))
      Error("out of space");

   if (!(u = (long *) NTL_MALLOC(n,  sizeof(long), 0)))
      Error("out of space");

   ZZ_p_rem_struct_init(&rem_struct, n, p, FFTPrime);

   ZZ_p_crt_struct_init(&crt_struct, n, p, FFTPrime);

   if (ZZ_p_crt_struct_special(crt_struct)) return;

   ZZ qq, rr;

   DivRem(qq, rr, M, p);

   NegateMod(MinusMModP, rr, p);

   for (i = 0; i < n; i++) {
      q = FFTPrime[i];

      long tt = rem(qq, q);

      mul(M2, p, tt);
      add(M2, M2, rr); 
      div(M2, M2, q);  // = (M/q) rem p
      

      div(M1, M, q);
      t = rem(M1, q);
      t = InvMod(t, q);

      mul(M3, M2, t);
      rem(M3, M3, p);

      ZZ_p_crt_struct_insert(crt_struct, i, M3);


      x[i] = ((double) t)/((double) q);
      u[i] = t;
   }
}



ZZ_pInfoT::~ZZ_pInfoT()
{
   long i;

   for (i = 0; i < MAX_ZZ_p_TEMPS; i++)
      if (temps[i]) delete temps[i];

   if (initialized) {
      ZZ_p_rem_struct_free(rem_struct);
      ZZ_p_crt_struct_free(crt_struct);

      free(x);
      free(u);
   }
}


ZZ_pInfoT *ZZ_pInfo = 0; 

typedef ZZ_pInfoT *ZZ_pInfoPtr;


static 
void CopyPointer(ZZ_pInfoPtr& dst, ZZ_pInfoPtr src)
{
   if (src == dst) return;

   if (dst) {
      dst->ref_count--;

      if (dst->ref_count < 0) 
         Error("internal error: negative ZZ_pContext ref_count");

      if (dst->ref_count == 0) delete dst;
   }

   if (src) {
      if (src->ref_count == NTL_MAX_LONG)
         Error("internal error: ZZ_pContext ref_count overflow");

      src->ref_count++;
   }

   dst = src;
}
   


void ZZ_p::init(const ZZ& p)
{
   ZZ_pContext c(p);
   c.restore();
}


ZZ_pContext::ZZ_pContext(const ZZ& p)
{
   ptr = NTL_NEW_OP ZZ_pInfoT(p);
}

ZZ_pContext::ZZ_pContext(const ZZ_pContext& a)
{
   ptr = 0;
   CopyPointer(ptr, a.ptr);
}

ZZ_pContext& ZZ_pContext::operator=(const ZZ_pContext& a)
{
   CopyPointer(ptr, a.ptr);
   return *this;
}


ZZ_pContext::~ZZ_pContext()
{
   CopyPointer(ptr, 0);
}

void ZZ_pContext::save()
{
   CopyPointer(ptr, ZZ_pInfo);
}

void ZZ_pContext::restore() const
{
   CopyPointer(ZZ_pInfo, ptr);
}



ZZ_pBak::~ZZ_pBak()
{
   if (MustRestore)
      CopyPointer(ZZ_pInfo, ptr);

   CopyPointer(ptr, 0);
}

void ZZ_pBak::save()
{
   MustRestore = 1;
   CopyPointer(ptr, ZZ_pInfo);
}


void ZZ_pBak::restore()
{
   MustRestore = 0;
   CopyPointer(ZZ_pInfo, ptr);
}


ZZ_pTemp::ZZ_pTemp()
{
   if (ZZ_pInfo->temps_top == MAX_ZZ_p_TEMPS)
      Error("ZZ_p temporary: out of temps");

   pos = ZZ_pInfo->temps_top;
   ZZ_pInfo->temps_top++;
}

ZZ_pTemp::~ZZ_pTemp()
{
   ZZ_pInfo->temps_top--;
}

ZZ_p& ZZ_pTemp::val() const
{
   if (!ZZ_pInfo->temps[pos]) 
      ZZ_pInfo->temps[pos] = NTL_NEW_OP ZZ_p;

   return *(ZZ_pInfo->temps[pos]);
}




const ZZ_p& ZZ_p::zero()
{
   static ZZ_p z(ZZ_p_NoAlloc);
   return z;
}

ZZ_p::DivHandlerPtr ZZ_p::DivHandler = 0;

ZZ_p::ZZ_p()
{
   _ZZ_p__rep.SetSize(ModulusSize());
}
   

ZZ_p::ZZ_p(INIT_VAL_TYPE, const ZZ& a) 
{
   _ZZ_p__rep.SetSize(ModulusSize());
   conv(*this, a);
} 

ZZ_p::ZZ_p(INIT_VAL_TYPE, long a)
{
   _ZZ_p__rep.SetSize(ModulusSize());
   conv(*this, a);
}


void conv(ZZ_p& x, long a)
{
   if (a == 0)
      clear(x);
   else if (a == 1)
      set(x);
   else {
      static ZZ y;

      conv(y, a);
      conv(x, y);
   }
}

istream& operator>>(istream& s, ZZ_p& x)
{
   static ZZ y;

   s >> y;
   conv(x, y);

   return s;
}

void div(ZZ_p& x, const ZZ_p& a, const ZZ_p& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val(); 

   inv(T, b);
   mul(x, a, T);
}

void inv(ZZ_p& x, const ZZ_p& a)
{
   if (InvModStatus(x._ZZ_p__rep, a._ZZ_p__rep, ZZ_p::modulus())) {
      if (IsZero(a._ZZ_p__rep))
         Error("ZZ_p: division by zero");
      else if (ZZ_p::DivHandler)
         (*ZZ_p::DivHandler)(a);
      else
         Error("ZZ_p: division by non-invertible element");
   }
}

long operator==(const ZZ_p& a, long b)
{
   if (b == 0)
      return IsZero(a);

   if (b == 1)
      return IsOne(a);

   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   return a == T;
}



void add(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   add(x, a, T);
}

void sub(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   sub(x, a, T);
}

void sub(ZZ_p& x, long a, const ZZ_p& b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, a);
   sub(x, T, b);
}

void mul(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   mul(x, a, T);
}

void div(ZZ_p& x, const ZZ_p& a, long b)
{
   ZZ_pTemp TT; ZZ_p& T = TT.val();
   conv(T, b);
   div(x, a, T);
}

void div(ZZ_p& x, long a, const ZZ_p& b)
{
   if (a == 1) {
      inv(x, b);
   }
   else {
      ZZ_pTemp TT; ZZ_p& T = TT.val();
      conv(T, a);
      div(x, T, b);
   }
}

NTL_END_IMPL
