#ifndef NTL_GF2XVec__H
#define NTL_GF2XVec__H

#include <NTL/GF2X.h>

NTL_OPEN_NNS



/*****************************************************************

The class GF2XVec implements vectors of fixed-length GF2X's.
You can allocate a vector of GF2X's of a specified length, where
the maximum size of each GF2X is also specified.
These parameters can be specified once, either with a constructor,
or with SetSize.
It is an error to try to re-size a vector, or store a GF2X that
doesn't fit.
The space can be released with "kill", and then you are free to 
call SetSize again.
If you want more flexible---but less efficient---vectors, 
use vec_GF2X.

*****************************************************************/



class GF2XVec {

private:
   GF2X* v;
   long len;
   long bsize;


public:
   GF2XVec& operator=(const GF2XVec&); 
   GF2XVec(const GF2XVec&); 

   long length() const { return len; }
   long BaseSize() const { return bsize; }
   void SetSize(long n, long d);
   void kill();

   GF2XVec() { v = 0; len = 0; bsize = 0; }
   GF2XVec(long n, long d) { v = 0; len = 0; bsize = 0; SetSize(n, d); }
   ~GF2XVec() { kill(); };

   GF2X* elts() { return v; }
   const GF2X* elts() const { return v; }

   GF2X& operator[](long i) { return v[i]; }
   const GF2X& operator[](long i) const { return v[i]; }

   static void swap_impl(GF2XVec& x, GF2XVec& y);
};

inline void swap(GF2XVec& x, GF2XVec& y)
   { GF2XVec::swap_impl(x, y); }

NTL_CLOSE_NNS

#endif
