/* ErrorChecker.h -- (C) Mark Rodenkirch, February 2012

   This class provides the interface for error checking after calls to OpenCL APIs.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _ERRORCHECKER_H
#define _ERRORCHECKER_H

#include "../core/main.h"

#ifdef __APPLE__
  #include <OpenCL/cl.h>
#else
#ifdef __linux
  #include <CL/cl.h>
#else
  #include <CL/cl.h>
#endif
#endif

class ErrorChecker
{
public:
   // If the status indicates an error, output the error and shutdown
   static void ExitIfError(const char *functionName, cl_int status);
   static void ExitIfError(const char *functionName, cl_int status, const char *fmt, ...) __attribute__ ((format (printf, 3, 4)));

private:
   static const char *GetErrorText(cl_int err);

};

#endif

