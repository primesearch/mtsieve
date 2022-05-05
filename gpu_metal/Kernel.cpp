/* Kernel.cpp -- (C) Mark Rodenkirch, April 2022

   This class provides the implementation for Metal kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "Kernel.h"
#include "ErrorChecker.h"
#include "../core/App.h"
#include "../core/Clock.h"

Kernel::Kernel(Device *device, const char *kernelName, const char *kernelSource, const char *preKernelSources[])
{
   NS::Error   *error;
   NS::String  *theSource = NS::String::string("#define USE_METAL", NS::UTF8StringEncoding);

   ip_Device = device;

   if (preKernelSources != NULL)
   {
      uint32_t ii = 0;
      while (preKernelSources[ii])
      {
         theSource = theSource.stringByAppendingString(NS::String::string(preKernelSources[ii], NS::UTF8StringEncoding));
         theSource = theSource.stringByAppendingString("\n", NS::UTF8StringEncoding));
         ii++;
      }
   }

   theSource = theSource.stringByAppendingString(NS::String::string(kernelSource, NS::UTF8StringEncoding));

   auto theLibrary = ip_Device->newLibrary(theSource, nullptr, &error);

   if (theLibrary == nullptr) {
      std::cerr << "Failed to create library with function " << functionName << std::endl;
      std::exit(-1);
   }

   auto strFunctionName = NS::String::string(functionName, NS::ASCIIStringEncoding);
   auto theFunction = theLibrary->newFunction(strFunctionName);

   if (theFunction == nullptr) {
      std::cerr << "Failed to find the function " << functionName << std::endl;
      std::exit(-1);
   }

   ip_ComputePipelineState = ip_Device->newComputePipelineState(theFunction, &error);
   ip_CommandQueue = ip_Device->newCommandQueue();

   ii_WorkGroupSize = ip_ComputePipelineState->maxTotalThreadsPerThreadgroup();

   ip_CommandBuffer = ip_CommandQueue->commandBuffer();
   ip_ComputeEncoder = ip_CommandBuffer->computeCommandEncoder();
}

Kernel::~Kernel(void)
{
}

void  *Kernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count);
}

void  *Kernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory)
{
   void *ptr = AddArgument(name, size, count);
   
   memcpy(ptr, cpuMemory, size * count);
   
   return ptr;
}

void  *Kernel::AddGpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count);
}

void  *Kernel::AddSharedArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count);
}

void  *Kernel::AddArgument(const char *name, uint32_t size, uint32_t count)
{
   ip_Device->IncrementGpuBytes(size * count);
   
   ip_Buffer[ii_BufferCount] = ip_Device->newBuffer(count * size, MTL::ResourceStorageModeShared);

   ii_BufferCount++;

   return ip_Buffer[ii_BufferCount - 1]->contents();
}

void Kernel::PrintStatistics(uint64_t bytesPerWorkGroup)
{
}

void Kernel::Execute(uint32_t workSize)
{
   uint64_t startTime;

   startTime = Clock::GetCurrentMicrosecond();
   
   ip_ComputeEncoder->setComputePipelineState(ip_ComputePipelineState);

   for (uint32_t idx=0; idx<=ii_BufferCount; idx++)
      ip_ComputeEncoder->setBuffer(ip_Buffer[idx], 0, idx);

   MTL::Size gridSize = MTL::Size(count, 1, 1);
   MTL::Size threadGroupSize = MTL::Size(ii_ThreadsPerGroup, 1, 1);

   // This shouldn't happen as the count is expected to be a multiple of ii_ThreadsPerGroup
   if (count < ii_ThreadsPerGroup) {
      threadGroupSize = MTL::Size(count, 1, 1);
   }

   ip_ComputeEncoder->dispatchThreads(gridSize, threadGroupSize);

   ip_ComputeEncoder->endEncoding();

   ip_CommandBuffer->commit();
   ip_CommandBuffer->waitUntilCompleted(); 
   
   ip_Device->AddGpuMicroseconds(Clock::GetCurrentMicrosecond() - startTime);
}
