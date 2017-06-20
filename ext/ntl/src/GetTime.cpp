#include <NTL/config.h>

#if (defined(NTL_CXX_ONLY) && !defined(__cplusplus))
#error "CXX_ONLY flag set...must use C++ compiler"
#endif

#include <time.h>


#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
extern "C" double _ntl_GetTime();
#endif


double _ntl_GetTime(void)
{
   static clock_t last_clock = 0;
   static double acc = 0;

   clock_t this_clock;
   double delta;

   this_clock = clock();

   delta = (this_clock - last_clock)/((double)CLOCKS_PER_SEC);
   if (delta < 0) delta = 0;

   acc += delta;
   last_clock = this_clock;

   return acc;
}

