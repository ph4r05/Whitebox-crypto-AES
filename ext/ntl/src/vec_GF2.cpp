
#include <NTL/vec_GF2.h>

#include <NTL/new.h>
#include <stdio.h>

NTL_START_IMPL

void vec_GF2::SetLength(long n)
{
   long len = length();

   if (n == len) return;

   if (n < 0) Error("negative length in vec_GF2::SetLength");

   if (NTL_OVERFLOW(n, 1, 0))
      Error("vec_GF2::SetLength: excessive length");

   if (fixed()) Error("SetLength: can't change this vector's length");

   long wdlen = (n+NTL_BITS_PER_LONG-1)/NTL_BITS_PER_LONG;

   if (n < len) {
      // have to clear bits n..len-1

      long q = n/NTL_BITS_PER_LONG;
      long p = n - q*NTL_BITS_PER_LONG;

      _ntl_ulong *x = rep.elts();

      x[q] &= (1UL << p) - 1UL;

      long q1 = (len-1)/NTL_BITS_PER_LONG;
      long i;

      for (i = q+1; i <= q1; i++)
         x[i] = 0;

      _len = n;

      rep.QuickSetLength(wdlen);

      return;
   }

   long maxlen = MaxLength();

   if (n <= maxlen) {
      _len = n;
      rep.QuickSetLength(wdlen);
      return;
   }

   long alloc = rep.MaxLength();

   if (wdlen <= alloc) {
      _len = n;
      _maxlen = (n << 1);
      rep.QuickSetLength(wdlen);
      return;
   }

   // have to grow vector and initialize to zero

   rep.SetLength(wdlen);

   wdlen = rep.MaxLength(); // careful! rep.MaxLength() may exceed the
                            // old value of wdlen...this is due to 
                            // the awkward semantics of WordVector.

   _ntl_ulong *x = rep.elts();

   long i;
   for (i = alloc; i < wdlen; i++)
      x[i] = 0;

   _len = n;
   _maxlen = (n << 1);

}


vec_GF2& vec_GF2::operator=(const vec_GF2& a)
{
   if (this == &a) return *this;

   long n = a.length();

   SetLength(n);

   long wdlen = (n+NTL_BITS_PER_LONG-1)/NTL_BITS_PER_LONG;

   _ntl_ulong *x = rep.elts();
   const _ntl_ulong *y = a.rep.elts();

   long i;
   for (i = 0; i < wdlen; i++)
      x[i] = y[i];

   return *this;
}

void vec_GF2::kill()
{
   if (fixed()) Error("can't kill this vec_GF2");
   rep.kill();
   _len = _maxlen = 0;
}


void vec_GF2::SetMaxLength(long n)
{
   long oldlen = length();
   if (n > oldlen) {
      SetLength(n);
      SetLength(oldlen);
   }
}

void vec_GF2::FixLength(long n)
{
   if (MaxLength() > 0 || fixed()) Error("can't fix this vector");

   SetLength(n);
   _maxlen |= 1;
}

const GF2 vec_GF2::get(long i) const
{
   const vec_GF2& v = *this;

   if (i < 0 || i >= v.length()) 
      Error("vec_GF2: subscript out of range");

   long q = i/NTL_BITS_PER_LONG;
   long p = i - q*NTL_BITS_PER_LONG;

   if (v.rep[q] & (1UL << p))
      return to_GF2(1);
   else
      return to_GF2(0);
}

ref_GF2 vec_GF2::operator[](long i)
{
   vec_GF2& v = *this;

   if (i < 0 || i >= v.length()) 
      Error("vec_GF2: subscript out of range");

   long q = i/NTL_BITS_PER_LONG;
   long p =  i - q*NTL_BITS_PER_LONG;
   return ref_GF2(INIT_LOOP_HOLE, &v.rep[q], p);
}



static
void SetBit(vec_GF2& v, long i)
{
   if (i < 0 || i >= v.length())
      Error("vec_GF2: subscript out of range");

   long q = i/NTL_BITS_PER_LONG;
   long p = i - q*NTL_BITS_PER_LONG;

   v.rep[q] |= (1UL << p);
}

static
void ClearBit(vec_GF2& v, long i)
{
   if (i < 0 || i >= v.length())
      Error("vec_GF2: subscript out of range");

   long q = i/NTL_BITS_PER_LONG;
   long p = i - q*NTL_BITS_PER_LONG;

   v.rep[q] &= ~(1UL << p);
}

void vec_GF2::put(long i, GF2 a)
{
   if (a == 1)
      SetBit(*this, i);
   else
      ClearBit(*this, i);
}

void swap(vec_GF2& x, vec_GF2& y)
{
   long xf = x.fixed();
   long yf = y.fixed();

   if (xf != yf || (xf && x.length() != y.length()))
      Error("swap: can't swap these vec_GF2s");

   swap(x.rep, y.rep);
   swap(x._len, y._len);
   swap(x._maxlen, y._maxlen);
}


void append(vec_GF2& v, const GF2& a)
{
   long n = v.length();
   v.SetLength(n+1);
   v.put(n, a);
}

void append(vec_GF2& x, const vec_GF2& a)
{
   long a_len = a.length();
   long x_len = x.length();

   if (a_len == 0) return;
   if (x_len == 0) {
      x = a;
      return;
   }


   x.SetLength(x_len + a_len);
   // new bits are guaranteed zero


   ShiftAdd(x.rep.elts(), a.rep.elts(), a.rep.length(), x_len);
}


long operator==(const vec_GF2& a, const vec_GF2& b)
{
   return a.length() == b.length() && a.rep == b.rep;
}





istream & operator>>(istream& s, vec_GF2& a) 
{   
   static ZZ ival;

   long c;   
   if (!s) Error("bad vec_GF2 input"); 
   
   c = s.peek();  
   while (IsWhiteSpace(c)) {  
      s.get();  
      c = s.peek();  
   }  

   if (c != '[') {  
      Error("bad vec_GF2 input");  
   }  

   vec_GF2 ibuf;  
   
   ibuf.SetLength(0);
      
   s.get();  
   c = s.peek();  
   while (IsWhiteSpace(c)) {  
      s.get();  
      c = s.peek();  
   }  

   while (c != ']' && c != EOF) {   
      if (!(s >> ival)) Error("bad vec_GF2 input");
      append(ibuf, to_GF2(ival));

      c = s.peek();  

      while (IsWhiteSpace(c)) {  
         s.get();  
         c = s.peek();  
      }  
   }   

   if (c == EOF) Error("bad vec_GF2 input");  
   s.get(); 
   
   a = ibuf; 
   return s;   
}    


ostream& operator<<(ostream& s, const vec_GF2& a)   
{   
   long i, n;
   GF2 c;
  
   n = a.length();
   
   s << '[';   
   
   for (i = 0; i < n; i++) {   
      c = a.get(i);
      if (c == 0)
         s << "0";
      else
         s << "1";
      if (i < n-1) s << " ";   
   }   
   
   s << ']';   
      
   return s;   
}   

// math operations:

void mul(vec_GF2& x, const vec_GF2& a, GF2 b)
{
   x = a;
   if (b == 0)
      clear(x);
}

void add(vec_GF2& x, const vec_GF2& a, const vec_GF2& b)
{
   long blen = a.length();

   if (b.length() != blen) Error("vec_GF2 add: length mismatch");

   x.SetLength(blen);

   long wlen = a.rep.length();
   long i;

   _ntl_ulong *xp = x.rep.elts();
   const _ntl_ulong *ap = a.rep.elts();
   const _ntl_ulong *bp = b.rep.elts();

   for (i = 0; i < wlen; i++)
      xp[i] = ap[i] ^ bp[i];
}

void clear(vec_GF2& x)
{
   long wlen = x.rep.length();
   long i;
   _ntl_ulong *xp = x.rep.elts();

   for (i = 0; i < wlen; i++)
      xp[i] = 0;
}


long IsZero(const vec_GF2& x)
{
   long wlen = x.rep.length();
   long i;
   const _ntl_ulong *xp = x.rep.elts();

   for (i = 0; i < wlen; i++)
      if (xp[i] != 0) return 0;

   return 1;
}

vec_GF2 operator+(const vec_GF2& a, const vec_GF2& b)
{
   vec_GF2 res;
   add(res, a, b);
   NTL_OPT_RETURN(vec_GF2, res);
}


vec_GF2 operator-(const vec_GF2& a, const vec_GF2& b)
{
   vec_GF2 res;
   add(res, a, b);
   NTL_OPT_RETURN(vec_GF2, res);
}

static
void ShiftToHigh(vec_GF2& x, const vec_GF2& a, long n)
// assumes 0 <= n < a.length()

{
   long l = a.length();

   x.SetLength(l);

   _ntl_ulong *xp = x.rep.elts();
   const _ntl_ulong *ap = a.rep.elts();

   long wn = n/NTL_BITS_PER_LONG;
   long bn = n - wn*NTL_BITS_PER_LONG;

   long sa = a.rep.length();

   long i;

   if (bn == 0) {
      for (i = sa-1; i >= wn; i--)
         xp[i] = ap[i-wn];
      for (i = wn-1; i >= 0; i--)
         xp[i] = 0; 
   }
   else {
      for (i = sa-1; i >= wn+1; i--)
         xp[i] = (ap[i-wn] << bn) | (ap[i-wn-1] >> (NTL_BITS_PER_LONG-bn));
      xp[wn] = ap[0] << bn;
      for (i = wn-1; i >= 0; i--)
         xp[i] = 0;
   }

   long p = l % NTL_BITS_PER_LONG;

   if (p != 0)
      xp[sa-1] &= (1UL << p) - 1UL;
   
}

static
void ShiftToLow(vec_GF2& x, const vec_GF2& a, long n)
// assumes 0 <= n < a.length()

{
   long l = a.length();

   x.SetLength(l);

   _ntl_ulong *xp = x.rep.elts();
   const _ntl_ulong *ap = a.rep.elts();

   long wn = n/NTL_BITS_PER_LONG;
   long bn = n - wn*NTL_BITS_PER_LONG;

   long sa = a.rep.length();

   long i;

   if (bn == 0) {
      for (i = 0; i < sa-wn; i++)
         xp[i] = ap[i+wn];
   }
   else {
      for (i = 0; i < sa-wn-1; i++)
         xp[i] = (ap[i+wn] >> bn) | (ap[i+wn+1] << (NTL_BITS_PER_LONG - bn));

      xp[sa-wn-1] = ap[sa-1] >> bn;
   }

   for (i = sa-wn; i < sa; i++)
      xp[i] = 0;
}



void shift(vec_GF2& x, const vec_GF2& a, long n)
{
   long l = a.length();

   if (n >= l || n <= -l) {
      x.SetLength(l);
      clear(x);
   }
   else if (n < 0) 
      ShiftToLow(x, a, -n); // |n| < l, so -n won't overflow!
   else
      ShiftToHigh(x, a, n);
}





// This code is simply canibalized from GF2X.c...
// so much for "code re-use" and "modularity"

static _ntl_ulong revtab[256] = {

0UL, 128UL, 64UL, 192UL, 32UL, 160UL, 96UL, 224UL, 16UL, 144UL, 
80UL, 208UL, 48UL, 176UL, 112UL, 240UL, 8UL, 136UL, 72UL, 200UL, 
40UL, 168UL, 104UL, 232UL, 24UL, 152UL, 88UL, 216UL, 56UL, 184UL, 
120UL, 248UL, 4UL, 132UL, 68UL, 196UL, 36UL, 164UL, 100UL, 228UL, 
20UL, 148UL, 84UL, 212UL, 52UL, 180UL, 116UL, 244UL, 12UL, 140UL, 
76UL, 204UL, 44UL, 172UL, 108UL, 236UL, 28UL, 156UL, 92UL, 220UL, 
60UL, 188UL, 124UL, 252UL, 2UL, 130UL, 66UL, 194UL, 34UL, 162UL, 
98UL, 226UL, 18UL, 146UL, 82UL, 210UL, 50UL, 178UL, 114UL, 242UL, 
10UL, 138UL, 74UL, 202UL, 42UL, 170UL, 106UL, 234UL, 26UL, 154UL, 
90UL, 218UL, 58UL, 186UL, 122UL, 250UL, 6UL, 134UL, 70UL, 198UL, 
38UL, 166UL, 102UL, 230UL, 22UL, 150UL, 86UL, 214UL, 54UL, 182UL, 
118UL, 246UL, 14UL, 142UL, 78UL, 206UL, 46UL, 174UL, 110UL, 238UL, 
30UL, 158UL, 94UL, 222UL, 62UL, 190UL, 126UL, 254UL, 1UL, 129UL, 
65UL, 193UL, 33UL, 161UL, 97UL, 225UL, 17UL, 145UL, 81UL, 209UL, 
49UL, 177UL, 113UL, 241UL, 9UL, 137UL, 73UL, 201UL, 41UL, 169UL, 
105UL, 233UL, 25UL, 153UL, 89UL, 217UL, 57UL, 185UL, 121UL, 249UL, 
5UL, 133UL, 69UL, 197UL, 37UL, 165UL, 101UL, 229UL, 21UL, 149UL, 
85UL, 213UL, 53UL, 181UL, 117UL, 245UL, 13UL, 141UL, 77UL, 205UL, 
45UL, 173UL, 109UL, 237UL, 29UL, 157UL, 93UL, 221UL, 61UL, 189UL, 
125UL, 253UL, 3UL, 131UL, 67UL, 195UL, 35UL, 163UL, 99UL, 227UL, 
19UL, 147UL, 83UL, 211UL, 51UL, 179UL, 115UL, 243UL, 11UL, 139UL, 
75UL, 203UL, 43UL, 171UL, 107UL, 235UL, 27UL, 155UL, 91UL, 219UL, 
59UL, 187UL, 123UL, 251UL, 7UL, 135UL, 71UL, 199UL, 39UL, 167UL, 
103UL, 231UL, 23UL, 151UL, 87UL, 215UL, 55UL, 183UL, 119UL, 247UL, 
15UL, 143UL, 79UL, 207UL, 47UL, 175UL, 111UL, 239UL, 31UL, 159UL, 
95UL, 223UL, 63UL, 191UL, 127UL, 255UL  }; 

static inline 
_ntl_ulong rev1(_ntl_ulong a)
{
   return NTL_BB_REV_CODE;
}



void reverse(vec_GF2& c, const vec_GF2& a)
// c = reverse of a

{
   long n = a.length();

   c = a;

   if (n <= 0) {
      return;
   }

   long wn = n/NTL_BITS_PER_LONG;
   long bn = n - wn*NTL_BITS_PER_LONG;

   if (bn != 0) {
      wn++;
      bn = NTL_BITS_PER_LONG - bn;
   }

   _ntl_ulong *cp = c.rep.elts();

   long i;

   if (bn != 0) {
      for (i = wn-1; i >= 1; i--)
         cp[i] = (cp[i] << bn) | (cp[i-1] >> (NTL_BITS_PER_LONG-bn));
      cp[0] = cp[0] << bn;
   }

   for (i = 0; i < wn/2; i++) {
      _ntl_ulong t; t = cp[i]; cp[i] = cp[wn-1-i]; cp[wn-1-i] = t;
   }

   for (i = 0; i < wn; i++)
      cp[i] = rev1(cp[i]);
}

static 
long weight1(_ntl_ulong a)
{
   long res = 0;
   while (a) {
      if (a & 1) res ++;
      a >>= 1;
   }
   return res;
}

long weight(const vec_GF2& a)
{
   long wlen = a.rep.length();
   long res = 0;
   long i;
   for (i = 0; i < wlen; i++)
      res += weight1(a.rep[i]);

   return res;
}

void random(vec_GF2& x, long n)
{
   if (n < 0) Error("random: bad arg");

   x.SetLength(n);

   long wl = x.rep.length();
   long i;

   for (i = 0; i < wl-1; i++) {
      x.rep[i] = RandomWord();
   }

   if (n > 0) {
      long pos = n % NTL_BITS_PER_LONG;
      if (pos == 0) pos = NTL_BITS_PER_LONG;
      x.rep[wl-1] = RandomBits_ulong(pos);
   }
}



void VectorCopy(vec_GF2& x, const vec_GF2& a, long n)
{
   if (n < 0) Error("VectorCopy: negative length");
   if (NTL_OVERFLOW(n, 1, 0)) Error("overflow in VectorCopy");

   long m = min(n, a.length());

   x.SetLength(n);

   long wn = (n + NTL_BITS_PER_LONG - 1)/NTL_BITS_PER_LONG;
   long wm = (m + NTL_BITS_PER_LONG - 1)/NTL_BITS_PER_LONG;

   _ntl_ulong *xp = x.rep.elts();
   const _ntl_ulong *ap = a.rep.elts();

   long i;

   for (i = 0; i < wm; i++)
      xp[i] = ap[i];

   for (i = wm; i < wn; i++)
      xp[i] = 0;

   long p = n % NTL_BITS_PER_LONG;
   if (p != 0) {
      xp[wn-1] &= ((1UL << p) - 1UL);
   }
}


NTL_END_IMPL
