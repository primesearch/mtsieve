/* Device.h -- (C) Mark Rodenkirch, May 2022

   This class provides the interface for Metal devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _DEVICE_H
#define _DEVICE_H

#include "../core/main.h"
#include "../core/GpuDevice.h"
#include "../core/GpuKernel.h"

class MetalDevice : public GpuDevice
{
public:
   MetalDevice(void);

   ~MetalDevice(void);

   void          *CreateKernel(const char *kernelName, const char *kernelSource, const char *preKernelSources[]);

   void           Help(void);
   void           AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   parse_t        ParseOption(int opt, char *arg, const char *source);
   void           ValidateOptions(void);

   void          *GetMetalDevice(void) { return ip_MetalDevice; };

private: 
   void          *ip_Pool;
   void          *ip_MetalDevice;
};

#endif

