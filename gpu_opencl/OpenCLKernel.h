/* OpenCLKernel.h -- (C) Mark Rodenkirch, May 2022

   This class provides the interface for OpenCL kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _OPENCL_KERNEL_H
#define _OPENCL_KERNEL_H

#include <string.h>

#include "OpenCLDevice.h"
#include "../core/GpuKernel.h"

#define MAX_KERNEL_ARGUMENTS  20

typedef struct {
   char              name[100];
   uint32_t          size;
   uint32_t          count;
   uint32_t          bytes;
   bool              mustFreeCpuBuffer;
   cl_mem_flags      memFlags;
   void             *cpuBuffer;
   cl_mem            gpuBuffer;
} ka_t;

class OpenCLKernel : public GpuKernel
{
public:
   // The last entry in this array needs to be NULL
   OpenCLKernel(GpuDevice *gpuDevice, const char *kernelName, const char *kernelSource, const char *preKernelSources[] = NULL);

   ~OpenCLKernel(void);

   void       PrintStatistics(uint64_t bytesPerWorkGroup);

   // Execute the function (kernel).  If necessary, it will copy host data to the GPU, execute,
   // then copy GPU data back to the host.  The idea is that the function calling this routine
   // doesn't need to know anything about the implementation of the call to the GPU.
   void       Execute(uint32_t workSize);

   // This adds an argument to the OpenCLKernel for memory that the CPU will write to
   // but that the GPU will only read.
   void      *AddCpuArgument(const char *name, uint32_t size, uint32_t count);
   void      *AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory);
   
   // This adds an argument to the OpenCLKernel for memory that the GPU will write to
   // but that the CPU will only read.
   void      *AddGpuArgument(const char *name, uint32_t size, uint32_t count);
   
   // This adds an argument to the OpenCLKernel for memory that bot the CPU and GPU
   // can read and write.
   void      *AddSharedArgument(const char *name, uint32_t size, uint32_t count);

private:
   void      *AddArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory, cl_mem_flags memFlags);

   void       SetGPUInput(void);
   void       GetGPUOutput(void);
   
   std::string       is_OpenCLKernelName;

   OpenCLDevice     *ip_OpenCLDevice;

   cl_program        im_Program;
   cl_kernel         im_OpenCLKernel;
   cl_command_queue  im_CommandQueue;

   uint32_t          ii_DeviceGlobalMemorySize;
   uint32_t          ii_DeviceLocalMemorySize;
   uint32_t          ii_WorkGroupSizeMultiple;
   uint32_t          ii_LocalMemorySize;
   uint32_t          ii_PrivateMemorySize;

   size_t            ii_KernelWorkGroupSize;
   
   uint32_t          ii_ArgumentCount;
   ka_t             *ip_KernelArguments;
};

#endif

