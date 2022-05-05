 /* Device.h -- (C) Mark Rodenkirch, April 2022

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

class Device
{
public:
   Device(void);

   ~Device(void);

   void           AddCommandLineOptions(string &shortOpts, struct option *longOpts);
   void           Help(void);
   parse_t        ParseOption(int opt, char *arg, const char *source);
   void           Validate(void);

   void           IncrementGpuBytes(int64_t bytes) { ip_GpuBytes->IncrementValue(bytes); };

   int64_t        GetGpuBytes(void) { return ip_GpuBytes->GetValueNoLock(); };

   void           AddGpuMicroseconds(int64_t totalMicroseconds) { ip_GpuMicroseconds->IncrementValue(totalMicroseconds); };
   uint64_t       GetGpuMicroseconds(void) { return ip_GpuMicroseconds->GetValueNoLock(); };
   
   bool           IsPrintDetails(void) { return ib_PrintDetails; };

private:
   SharedMemoryItem *ip_GpuBytes;
   SharedMemoryItem *ip_GpuMicroseconds;
 
   NS::AutoreleasePool        *ip_Pool;
   MTL::Device                *ip_Device;
};

#endif

