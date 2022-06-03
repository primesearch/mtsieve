 /* GpuDevice.h -- (C) Mark Rodenkirch, May 2022

   This class provides a generic interface for GPU devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _GPU_DEVICE_H
#define _GPU_DEVICE_H

#include "main.h"
#include "SharedMemoryItem.h"
#include "Parser.h"

class GpuDevice
{
public:
   GpuDevice(void);

   virtual ~GpuDevice(void) = 0;

   virtual void        *CreateKernel(const char *kernelName, const char *kernelSource, const char *preKernelSources[]) = 0;
   
   virtual void         Help(void) = 0;
   virtual void         AddCommandLineOptions(std::string &shortOpts, struct option *longOpts) = 0;
   virtual parse_t      ParseOption(int opt, char *arg, const char *source) = 0;
   virtual void         ValidateOptions(void) = 0;
   
   void        CleanUp(void);
   void        ParentAddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   void        ParentHelp(void);
   parse_t     ParentParseOption(int opt, char *arg, const char *source);
   void        ParentValidate(void);
               
   void        IncrementGpuBytes(int64_t bytes) { ip_GpuBytes->IncrementValue(bytes); };
               
   int64_t     GetGpuBytes(void) { return ip_GpuBytes->GetValueNoLock(); };
               
   void        AddGpuMicroseconds(int64_t totalMicroseconds) { ip_GpuMicroseconds->IncrementValue(totalMicroseconds); };
   uint64_t    GetGpuMicroseconds(void) { return ip_GpuMicroseconds->GetValueNoLock(); };

private:
   SharedMemoryItem *ip_GpuBytes;
   SharedMemoryItem *ip_GpuMicroseconds;
};

#endif

