/* Device.cpp -- (C) Mark Rodenkirch, April 2022

   This classs manages Metal devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdio.h>

#include "../core/main.h"

#include "Device.h"

Device::Device(void)
{
   ip_GpuBytes = new SharedMemoryItem("gpu_bytes");
   ip_GpuMicroseconds = new SharedMemoryItem("gpu_microseconds");
   
   ip_Pool = NS::AutoreleasePool::alloc()->init();
   ip_Device = MTL::CreateSystemDefaultDevice();

   ii_BufferCount = 0;
}

Device::~Device(void)
{
   ip_Pool->release();
}

void  Device::Help(void)
{
}

void  Device::AddCommandLineOptions(string &shortOpts, struct option *longOpts)
{
}

// Returns:
//    0 if the option is OK
//   -1 if the argument is invalid
//   -2 if the argument is out of range
//   99 if the argument is not supported by this module
parse_t Device::ParseOption(int opt, char *arg, const char *source)
{
  return P_UNSUPPORTED;
}
