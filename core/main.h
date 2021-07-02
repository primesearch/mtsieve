/* main.h -- (C) Mark Rodenkirch, November 2017.

   Single-threaded CPU sieve application framework.

   For each prime p in 3 <= p0 <= p < p1 < 2^62
     Do something with p
     
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _MAIN_H
#define _MAIN_H 1

#ifdef __linux__
   #define __STDC_FORMAT_MACROS
   #define __STDC_CONSTANT_MACROS
#endif

#include <iostream>
#include <string>
#include <string.h>
#include <math.h>
#include <cinttypes>

using namespace std;

#ifdef WIN32
   #if defined(_MSC_VER) && defined(MEMLEAK)
      #define _CRTDBG_MAP_ALLOC
      #include <stdlib.h>
      #include <crtdbg.h>

      #define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
      #define new DEBUG_NEW
   #endif

   #include "getopt.h"

   #define unlink _unlink
#else
   #include <getopt.h>
   #include <unistd.h>
   #include <errno.h>

   #define Sleep(x) usleep((x)*1000)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// In main.cpp
void     FatalError(const char *fmt, ...) __attribute__ ((format (printf, 1, 2)));
void     AppendLongOpt(struct option *longOpts, const char *name, int has_arg, int *flag, char charSwitch);

void    *xmalloc(size_t size);
void     xfree(void *mem);
uint64_t GetCpuMemoryUsage(void);

#ifdef __cplusplus
}
#endif

#endif /* _MAIN_H */
