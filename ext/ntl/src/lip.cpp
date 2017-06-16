
#include <NTL/config.h>

#if (defined(NTL_CXX_ONLY) && !defined(__cplusplus))
#error "CXX_ONLY flag set...must use C++ compiler"
#endif

#ifdef NTL_GMP_LIP

#include "g_lip_impl.h"

#else

#include "c_lip_impl.h"

#endif
