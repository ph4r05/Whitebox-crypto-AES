
#include <stdlib.h>

int val_int(int x) { return x; }
unsigned int val_uint(unsigned int x) { return x; }
 
long val_long(long x) { return x; }
unsigned long val_ulong(unsigned long x) { return x; }
 
size_t val_size_t(size_t x) { return x; }

double val_double(double x) { return x; }
 
void touch_int(int* x) {}
void touch_uint(unsigned int* x) {}
 
void touch_long(long* x) {}
void touch_ulong(unsigned long* x) {}

void touch_size_t(size_t* x) {}
 
void touch_double(double* x) {}

