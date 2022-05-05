/* MetalDevice.cpp -- (C) Mark Rodenkirch, February 2022

   This classs manages Metal devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#include <iostream>
#include "MetalDevice.h"

MetalDevice::MetalDevice(void) 
{
   ip_Pool = NS::AutoreleasePool::alloc()->init();
   ip_Device = MTL::CreateSystemDefaultDevice();

   ii_BufferCount = 0;
}

MetalDevice::~MetalDevice(void)
{
   ip_Pool->release();
}

void MetalDevice::InitializeFunction(const char *source, char *functionName)
{
   NS::Error   *error;
   NS::String  *theSource = NS::String::string(source, NS::UTF8StringEncoding);

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

   ii_ThreadsPerGroup = ip_ComputePipelineState->maxTotalThreadsPerThreadgroup();

   ip_CommandBuffer = ip_CommandQueue->commandBuffer();
   ip_ComputeEncoder = ip_CommandBuffer->computeCommandEncoder();
}

void *MetalDevice::CreateBuffer(uint32_t count, uint32_t size)
{
   ip_Buffer[ii_BufferCount] = ip_Device->newBuffer(count * size, MTL::ResourceStorageModeShared);

   ii_BufferCount++;

   return ip_Buffer[ii_BufferCount - 1]->contents();
}

void MetalDevice::Execute(uint32_t count)
{
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
}
