/* MetalDevice.hpp -- (C) Mark Rodenkirch, April 2022

   This class provides the interface for Metal devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _METALDEVICE_H
#define _METALDEVICE_H

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include "../core/main.h"

#define MAX_BUFFERS  20

class MetalDevice
{
public:
   MetalDevice(void);

   ~MetalDevice(void);

   void       InitializeFunction(const char *source, char *functionName);
   void      *CreateBuffer(uint32_t count, uint32_t size);
   uint32_t   GetThreadsPerGroup(void) { return (uint32_t) ii_ThreadsPerGroup; };
   void       Execute(uint32_t count);

private:
   NS::AutoreleasePool        *ip_Pool;
   MTL::Device                *ip_Device;

   MTL::ComputePipelineState  *ip_ComputePipelineState;
   MTL::CommandQueue          *ip_CommandQueue; 

   MTL::CommandBuffer         *ip_CommandBuffer;
   MTL::ComputeCommandEncoder *ip_ComputeEncoder;

   NS::UInteger         ii_ThreadsPerGroup;

   uint32_t             ii_BufferCount;
   MTL::Buffer         *ip_Buffer[MAX_BUFFERS];
};

#endif
