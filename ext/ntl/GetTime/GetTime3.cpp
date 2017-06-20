#include <NTL/config.h>

#if (defined(NTL_CXX_ONLY) && !defined(__cplusplus))
#error "CXX_ONLY flag set...must use C++ compiler"
#endif

#include <sys/time.h>
#include <sys/resource.h>
#include <syscall.h>


#if (defined(__cplusplus) && !defined(NTL_CXX_ONLY))
extern "C" double _ntl_GetTime();
#endif


double _ntl_GetTime(void)
{
   struct rusage used;

   syscall(SYS_getrusage, RUSAGE_SELF, &used);
   return (used.ru_utime.tv_sec + used.ru_stime.tv_sec +
      (used.ru_utime.tv_usec + used.ru_stime.tv_usec) / 1e6);
}

