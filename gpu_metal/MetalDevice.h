 /* Device.h -- (C) Mark Rodenkirch, May 2022

   This class provides the interface for Metal devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _DEVICE_H
#define _DEVICE_H

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include "../core/main.h"
#include "../core/GpuDevice.h"
#include "../core/GpuKernel.h"

class MetalDevice : public GpuDevice
{
public:
   MetalDevice(void);

   ~MetalDevice(void);

   GpuKernel     *CreateKernel(GpuDevice *gpuDevice, const char *kernelName, const char *kernelSource, const char *preKernelSources[]);

   void           Help(void);
   void           AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t        ParseOption(int opt, char *arg, const char *source);
   void           ValidateOptions(void);

private: 
   NS::AutoreleasePool *ip_Pool;
   MTL::Device         *ip_MetalDevice;
};

#endif

