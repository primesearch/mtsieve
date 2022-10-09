/* MetalDevice.cpp -- (C) Mark Rodenkirch, May 2022

   This classs manages Metal devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <stdio.h>

#include "../core/main.h"

#include "MetalDevice.h"
#include "MetalKernel.h"

MetalDevice::MetalDevice(void)
{
   ip_Pool = NS::AutoreleasePool::alloc()->init();
   ip_MetalDevice = MTL::CreateSystemDefaultDevice();
}

MetalDevice::~MetalDevice(void)
{
   ((NS::AutoreleasePool *) ip_Pool)->release();
}

void *MetalDevice::CreateKernel(const char *kernelName, const char *kernelSource, const char *preKernelSources[])
{
   return new MetalKernel(this, kernelName, kernelSource, preKernelSources);
}

void  MetalDevice::Help(void)
{
   GpuDevice::ParentHelp();
}

void  MetalDevice::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   GpuDevice::ParentAddCommandLineOptions(shortOpts, longOpts);
}

// Returns:
//    0 if the option is OK
//   -1 if the argument is invalid
//   -2 if the argument is out of range
//   99 if the argument is not supported by this module
parse_t MetalDevice::ParseOption(int opt, char *arg, const char *source)
{
  return GpuDevice::ParentParseOption(opt, arg, source);
}

void  MetalDevice::ValidateOptions(void)
{
}
