/* Kernel.h -- (C) Mark Rodenkirch, February 2012

   This class provides the interface for OpenCL kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _KERNEL_H
#define _KERNEL_H

#include <string.h>

#include "Device.h"
#include "KernelArgument.h"

using namespace std;

#define MAX_KERNEL_ARGUMENTS 50

class Kernel
{
public:
   // The last entry in this array needs to be NULL
   Kernel(Device *device, const char *kernelName, const char *kernelSource[], bool useFMA = false);

   ~Kernel(void);

   // Add an function argument for this kernel
   void       AddArgument(KernelArgument *theArgument);

   // Execute the function (kernel).  If necessary, it will copy host data to the GPU, execute,
   // then copy GPU data back to the host.  The idea is that the function calling this routine
   // doesn't need to know anything about the implementation of the call to the GPU.
   void       Execute(uint32_t workSize);

   int        GetWorkGroupSize(void) { return ii_WorkGroupSize; };

   void       PrintStats(void);

private:
   string            is_KernelName;

   cl_program        im_Program;
   cl_kernel         im_Kernel;
   cl_command_queue  im_CommandQueue;

   int               ii_WorkGroupSize;
   int               ii_ArgumentCount;

   size_t            ii_KernelWorkGroupSize;
   
   Device           *ip_Device;
   KernelArgument   *ip_Arguments[MAX_KERNEL_ARGUMENTS];

   void       SetGPUInput(void);
   void       GetGPUOutput(void);
};

#endif

