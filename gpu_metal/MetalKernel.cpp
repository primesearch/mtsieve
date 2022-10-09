/* MetalKernel.cpp -- (C) Mark Rodenkirch, April 2022

   This class provides the implementation for Metal MetalKernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "MetalKernel.h"
#include "../core/App.h"
#include "../core/Clock.h"

MetalKernel::MetalKernel(GpuDevice *gpuDevice, const char *kernelName, const char *kernelSource, const char *preMetalKernelSources[])
   : GpuKernel(gpuDevice, kernelName, kernelSource, preMetalKernelSources)
{
   NS::Error   *error;
   NS::String  *theSource = NS::String::string("#define USE_METAL", NS::UTF8StringEncoding);

   ip_MetalDevice = ((MetalDevice *) gpuDevice)->GetMetalDevice();

   if (preMetalKernelSources != NULL)
   {
      uint32_t ii = 0;
      while (preMetalKernelSources[ii])
      {
         theSource = theSource->stringByAppendingString(NS::String::string(preMetalKernelSources[ii], NS::UTF8StringEncoding));
         theSource = theSource->stringByAppendingString(NS::String::string("\n", NS::UTF8StringEncoding));
         ii++;
      }
   }

   theSource = theSource->stringByAppendingString(NS::String::string(kernelSource, NS::UTF8StringEncoding));

   auto theLibrary = ip_MetalDevice->newLibrary(theSource, nullptr, &error);

   if (theLibrary == nullptr) {
      std::cerr << "Failed to create library with function " << kernelName << std::endl;
      std::exit(-1);
   }

   auto strFunctionName = NS::String::string(kernelName, NS::ASCIIStringEncoding);
   auto theFunction = theLibrary->newFunction(strFunctionName);

   if (theFunction == nullptr) {
      std::cerr << "Failed to find the function " << kernelName << std::endl;
      std::exit(-1);
   }

   ip_ComputePipelineState = ip_MetalDevice->newComputePipelineState(theFunction, &error);
   ip_CommandQueue = ip_MetalDevice->newCommandQueue();

   ii_ThreadsPerGroup = ip_ComputePipelineState->maxTotalThreadsPerThreadgroup();

   ip_CommandBuffer = ip_CommandQueue->commandBuffer();
   ip_ComputeEncoder = ip_CommandBuffer->computeCommandEncoder();
}

MetalKernel::~MetalKernel(void)
{
}

void MetalKernel::PrintStatistics(uint64_t bytesPerWorkGroup)
{
}

void    *MetalKernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, NULL);
}

void    *MetalKernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory)
{
   return AddArgument(name, size, count, cpuMemory);
}

void    *MetalKernel::AddGpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, NULL);
}

void    *MetalKernel::AddSharedArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, NULL);
}

void  *MetalKernel::AddArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory)
{
   if (cpuMemory == NULL)
      cpuMemory = xmalloc(size * count + 1);

   ip_GpuDevice->IncrementGpuBytes(size * count);
   
   ip_Buffer[ii_BufferCount] = ip_MetalDevice->newBuffer(count * size, MTL::ResourceStorageModeShared);

   ii_BufferCount++;

   void *buffer = ip_Buffer[ii_BufferCount - 1]->contents();
   
   if (cpuMemory != NULL)
      memcpy(buffer, cpuMemory, count * size);

   return buffer;
}

void MetalKernel::Execute(uint32_t workSize)
{
   uint64_t startTime;

   startTime = Clock::GetCurrentMicrosecond();
   
   ip_ComputeEncoder->setComputePipelineState(ip_ComputePipelineState);

   for (uint32_t idx=0; idx<=ii_BufferCount; idx++)
      ip_ComputeEncoder->setBuffer(ip_Buffer[idx], 0, idx);

   MTL::Size gridSize = MTL::Size(workSize, 1, 1);
   MTL::Size threadGroupSize = MTL::Size(ii_ThreadsPerGroup, 1, 1);

   // This shouldn't happen as the count is expected to be a multiple of ii_ThreadsPerGroup
   if (workSize < ii_ThreadsPerGroup) {
      threadGroupSize = MTL::Size(workSize, 1, 1);
   }

   ip_ComputeEncoder->dispatchThreads(gridSize, threadGroupSize);

   ip_ComputeEncoder->endEncoding();

   ip_CommandBuffer->commit();
   ip_CommandBuffer->waitUntilCompleted(); 
   
   ip_GpuDevice->AddGpuMicroseconds(Clock::GetCurrentMicrosecond() - startTime);
}
