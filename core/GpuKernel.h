/* GpuKernel.h -- (C) Mark Rodenkirch, May 2022

   This class provides a generic interface for GPU kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GPU_KERNEL_H
#define _GPU_KERNEL_H

#include <string.h>

#include "GpuDevice.h"

#define MAX_BUFFERS  20

class GpuKernel
{
public:
   // The last entry in this array needs to be NULL
   GpuKernel(GpuDevice *gpuDevice, const char *kernelName, const char *kernelSource, const char *preKernelSources[] = NULL);

   virtual ~GpuKernel(void) = 0;

   virtual void       PrintStatistics(uint64_t bytesPerWorkGroup) = 0;

   // Execute the kernel.  If necessary, it will copy host data to the GPU, execute,
   // then copy GPU data back to the host.
   virtual void       Execute(uint32_t workSize) = 0;
   // This adds an argument to the kernel for memory that the CPU will write to
   // but that the GPU will only read.
   virtual void      *AddCpuArgument(const char *name, uint32_t size, uint32_t count) = 0;
   virtual void      *AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory) = 0;
   
   // This adds an argument to the kernel for memory that the GPU will write to
   // but that the CPU will only read.
   virtual void      *AddGpuArgument(const char *name, uint32_t size, uint32_t count) = 0;
   
   // This adds an argument to the kernel for memory that bot the CPU and GPU
   // can read and write.
   virtual void      *AddSharedArgument(const char *name, uint32_t size, uint32_t count) = 0;

   uint32_t   GetWorkGroupSize(void) { return ii_WorkGroupSize; };

protected:
   GpuDevice     *ip_GpuDevice;

   std::string    is_KernelName;
      
   uint32_t       ii_WorkGroupSize;
};

#endif

