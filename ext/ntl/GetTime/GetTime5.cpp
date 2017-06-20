#include <NTL/config.h>

#if (defined(NTL_CXX_ONLY) && !defined(__cplusplus))
#error "CXX_ONLY flag set...must use C++ compiler"
#endif


#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
extern "C" double _ntl_GetTime();
#endif


double _ntl_GetTime(void)
{
   return 0;
}
