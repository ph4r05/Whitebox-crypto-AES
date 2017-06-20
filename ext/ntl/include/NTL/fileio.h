
#ifndef NTL_fileio__H
#define NTL_fileio__H

#include <NTL/tools.h>

#if (defined(NTL_STD_CXX) || defined(NTL_PSTD_NHF))

#include <fstream>                                                              

#else

#include <fstream.h>

#endif

#if 0
namespace foo_bar {

class ofstream;
class ifstream;

}
#endif

NTL_OPEN_NNS


void OpenWrite(NTL_SNS ofstream& s, const char *name);

// opens file for writing...aborts if fails

void OpenRead(NTL_SNS ifstream& s, const char *name);

// opens file for reading

char *FileName(const char* stem, const char *ext);

// builds the name "stem.ext"

char *FileName(const char* stem, const char *ext, long d);

// builds the name stem.ext.d

NTL_CLOSE_NNS

#endif


