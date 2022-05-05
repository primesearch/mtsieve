/* Kernel.h -- (C) Mark Rodenkirch, Aprile 2022

   This class provides the interface for Metal kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _KERNEL_H
#define _KERNEL_H

#include <string.h>

#include "Device.h"

#define MAX_BUFFERS  20

class Kernel
{
public:
   // The last entry in this array needs to be NULL
   Kernel(Device *device, const char *kernelName, const char *kernelSource, const char *preKernelSources[] = NULL);

   ~Kernel(void);

   Device    *GetDevice(void) { return ip_Device; };

   void       PrintStatistics(uint64_t bytesPerWorkGroup);

   // Execute the function (kernel).  If necessary, it will copy host data to the GPU, execute,
   // then copy GPU data back to the host.  The idea is that the function calling this routine
   // doesn't need to know anything about the implementation of the call to the GPU.
   void       Execute(uint32_t workSize);

   uint32_t   GetWorkGroupSize(void) { return ii_WorkGroupSize; };

   // This adds an argument to the Kernel for memory that the CPU will write to
   // but that the GPU will only read.
   void      *AddCpuArgument(const char *name, uint32_t size, uint32_t count);
   void      *AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory);
   
   // This adds an argument to the Kernel for memory that the GPU will write to
   // but that the CPU will only read.
   void      *AddGpuArgument(const char *name, uint32_t size, uint32_t count);
   
   // This adds an argument to the Kernel for memory that bot the CPU and GPU
   // can read and write.
   void      *AddSharedArgument(const char *name, uint32_t size, uint32_t count);

private:
   Device    *ip_Device;
   
   void      *AddArgument(const char *name, uint32_t size, uint32_t count);
   
   MTL::ComputePipelineState  *ip_ComputePipelineState;
   MTL::CommandQueue          *ip_CommandQueue; 

   MTL::CommandBuffer         *ip_CommandBuffer;
   MTL::ComputeCommandEncoder *ip_ComputeEncoder;
   
   uint32_t             ii_BufferCount;
   MTL::Buffer         *ip_Buffer[MAX_BUFFERS];
   
   uint32_t             ii_WorkGroupSize;
};

#endif

