 /* OpenCLDevice.h -- (C) Mark Rodenkirch, May 2022

   This class provides the interface for OpenCL devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _DEVICE_H
#define _DEVICE_H

#ifdef __APPLE__
   #include <OpenCL/cl.h>
   #include <OpenCL/cl_ext.h>
#else
#ifdef __linux
   #include <CL/cl.h>
#else
   #include <CL/cl.h>
#endif
#endif

#include "../core/main.h"
#include "../core/SharedMemoryItem.h"
#include "../core/Parser.h"
#include "../core/GpuDevice.h"
#include "../core/GpuKernel.h"

typedef struct
{
   cl_device_id   deviceId;
   char           vendor[1000];
   char           name[1000];
   cl_uint        maxComputeUnits;
   cl_ulong       localMemSize;
   cl_device_fp_config fpConfig;
} device_t;

typedef struct
{
   cl_platform_id platformId;
   cl_context     context;
   char           vendor[1000];
   char           version[1000];
   char           name[1000];
   cl_uint        deviceCount;
   device_t      *devices;
} platform_t;

class OpenCLDevice : public GpuDevice
{
public:
   OpenCLDevice(void);

   ~OpenCLDevice(void);

   void          *CreateKernel(const char *kernelName, const char *kernelSource, const char *preKernelSources[]);
   
   // List all platforms and devices available.  Note that just because a platform/device is
   // listed, it does not mean that it can be used by the framework
   void           ListAllOpenCLDevices(void);

   void           AddCommandLineOptions(std::string &shortOpts, struct option *longOpts);
   void           Help(void);
   parse_t        ParseOption(int opt, char *arg, const char *source);
   void           ValidateOptions(void);

   cl_context     GetContext(void) { return ip_Platforms[ii_WhichPlatform].context; };
   cl_device_id   GetDeviceId(void) { return ip_Platforms[ii_WhichPlatform].devices[ii_WhichDevice].deviceId; };
   cl_device_id  *GetDeviceIdPtr(void) { return &ip_Platforms[ii_WhichPlatform].devices[ii_WhichDevice].deviceId; };
   cl_uint        GetMaxComputeUnits(void) { return ip_Platforms[ii_WhichPlatform].devices[ii_WhichDevice].maxComputeUnits; };
   
   bool           IsPrintDetails(void) { return ib_PrintDetails; };

private:
   void           GetPlatforms(void);
   void           GetDevicesForPlatform(platform_t *thePlatform);

   platform_t    *ip_Platforms;

   bool           ib_HavePlatform;
   bool           ib_HaveOpenCLDevice;
   bool           ib_PrintDetails;
   
   int32_t        ii_TotalDeviceCount;
   cl_uint        ii_PlatformCount;
   cl_uint        ii_WhichPlatform;
   cl_uint        ii_WhichDevice;
};

#endif

