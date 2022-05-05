/* Device.cpp -- (C) Mark Rodenkirch, February 2012

   This classs manages OpenCL platforms and devices.

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
#include "ErrorChecker.h"

Device::Device(void)
{
   cl_int status;

   ip_GpuBytes = new SharedMemoryItem("gpu_bytes");
   ip_GpuMicroseconds = new SharedMemoryItem("gpu_microseconds");

   ii_PlatformId = -1;
   ii_DeviceId = -1;
   ip_Platforms = 0;
   ii_TotalDeviceCount = 0;
   ib_HavePlatform = false;
   ib_HaveDevice = false;
   ib_PrintDetails = false;

   status = clGetPlatformIDs(0, NULL, &ii_PlatformCount);

   ErrorChecker::ExitIfError("clGetPlatformIDs", status, "%s\n%s", 
      "Please (re)install OpenCL as described at",
      "http://developer.amd.com/gpu/ATIStreamSDK/assets/ATI_Stream_SDK_Installation_Notes.pdf");

   GetPlatforms();

   if (!ii_TotalDeviceCount)
   {
      ListAllDevices();
      FatalError("\nNo devices were found that can run this code\n");
   }
}

Device::~Device(void)
{
   cl_uint  pp;
   cl_uint  ii;

   for (pp=0; pp<ii_PlatformCount; pp++)
      clReleaseContext(ip_Platforms[pp].context);

   delete ip_GpuBytes;
   delete ip_GpuMicroseconds;
   
   for (ii=0; ii<ii_PlatformCount; ii++)
      xfree(ip_Platforms[ii].devices);

   xfree(ip_Platforms);
}

void  Device::Help(void)
{
   printf("-D --platform=D       Use platform D instead of 0\n");
   printf("-d --device=d         Use device d instead of 0\n");
   printf("-H --showgpudetail    Show device and kernel details\n");
   
   ListAllDevices();
}

void  Device::AddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
   shortOpts += "HD:d:";

   AppendLongOpt(longOpts, "showgpudetail", no_argument,       0, 'H');
   AppendLongOpt(longOpts, "platform",      required_argument, 0, 'D');
   AppendLongOpt(longOpts, "device",        required_argument, 0, 'd');
}

// Returns:
//    0 if the option is OK
//   -1 if the argument is invalid
//   -2 if the argument is out of range
//   99 if the argument is not supported by this module
parse_t Device::ParseOption(int opt, char *arg, const char *source)
{
   parse_t      status = P_UNSUPPORTED;

   switch (opt)
   {
      case 'H':
         status = P_SUCCESS;
         ib_PrintDetails = true;
         break;

      case 'D':
         status = Parser::Parse(arg, 0, INT32_MAX, ii_PlatformId);
         ib_HavePlatform = true;
         break;

      case 'd':
         status = Parser::Parse(arg, 0, INT32_MAX, ii_DeviceId);
         ib_HaveDevice = true;
         break;
  }

  return status;
}

void Device::GetPlatforms(void)
{
   cl_int      status;
   cl_platform_id *platforms; 
   cl_uint     ii;

   platforms = (cl_platform_id *) xmalloc(sizeof(cl_platform_id) * ii_PlatformCount);

   status = clGetPlatformIDs(ii_PlatformCount, platforms, NULL);
   ErrorChecker::ExitIfError("clGetPlatformIDs", status, "Unable to get platforms");

   ip_Platforms = (platform_t *) xmalloc(sizeof(platform_t) * ii_PlatformCount);

   for (ii=0; ii<ii_PlatformCount; ii++)
   {
      ip_Platforms[ii].platformId = platforms[ii];
      ip_Platforms[ii].deviceCount = 0;
      ip_Platforms[ii].devices = 0;

      status = clGetPlatformInfo(platforms[ii], CL_PLATFORM_VENDOR, sizeof(ip_Platforms[ii].vendor), ip_Platforms[ii].vendor, NULL);
      ErrorChecker::ExitIfError("clGetPlatformInfo", status, "Unable to get vendor");

      status = clGetPlatformInfo(platforms[ii], CL_PLATFORM_VERSION, sizeof(ip_Platforms[ii].version), ip_Platforms[ii].version, NULL);
      ErrorChecker::ExitIfError("clGetPlatformInfo", status, "Unable to get version");

      status = clGetPlatformInfo(platforms[ii], CL_PLATFORM_NAME, sizeof(ip_Platforms[ii].name), ip_Platforms[ii].name, NULL);
      ErrorChecker::ExitIfError("clGetPlatformInfo", status, "Unable to get version");

      GetDevicesForPlatform(&ip_Platforms[ii]);
   }

   xfree(platforms);
}

void Device::GetDevicesForPlatform(platform_t *thePlatform)
{
   cl_int        status;
   cl_context_properties contextProperties[3];
   cl_uint       numDevices;
   cl_device_id *devices;
   device_t     *theDevice;
   cl_uint       ii;
   size_t        deviceListSize;

   contextProperties[0] = CL_CONTEXT_PLATFORM;
   contextProperties[1] = (cl_context_properties) thePlatform->platformId;
   contextProperties[2] = 0;

   thePlatform->deviceCount = 0;

   // We only want to run this on GPUs
   thePlatform->context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU, NULL, NULL, &status);
   if (status != CL_SUCCESS) return;
   ErrorChecker::ExitIfError("clCreateContextFromType", status, "Unable to create context");

   status = clGetContextInfo(thePlatform->context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize);
   ErrorChecker::ExitIfError("clGetContextInfo", status, "Unable to get device list size");

   status = clGetDeviceIDs(thePlatform->platformId, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
   ErrorChecker::ExitIfError("clGetDeviceIDs", status, "Unable to get number of devices");

   if (deviceListSize == 0 || numDevices == 0)
      return;

   devices = (cl_device_id *) xmalloc(deviceListSize);
 
   status = clGetContextInfo(thePlatform->context, CL_CONTEXT_DEVICES, deviceListSize, devices, NULL);
   ErrorChecker::ExitIfError("clGetContextInfo", status, "Unable to get device list");

   thePlatform->deviceCount = numDevices;
   thePlatform->devices = (device_t *) xmalloc(sizeof(device_t) * numDevices);

   ii_TotalDeviceCount += numDevices;

   for (ii=0; ii<numDevices; ii++)
   {
      thePlatform->devices[ii].deviceId = devices[ii];
      theDevice = &thePlatform->devices[ii];

      status = clGetDeviceInfo(devices[ii], CL_DEVICE_VENDOR, sizeof(theDevice->vendor), theDevice->vendor, NULL);
      ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get device vendor");

      status = clGetDeviceInfo(devices[ii], CL_DEVICE_NAME, sizeof(theDevice->name), theDevice->name, NULL);
      ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get device name");

#ifdef CL_DEVICE_DOUBLE_FP_CONFIG
      theDevice->fpConfig = 0;
      status = clGetDeviceInfo(devices[ii], CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(theDevice->fpConfig), &theDevice->fpConfig, NULL);
      ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_DOUBLE_FP_CONFIG of device");
#endif

      status = clGetDeviceInfo(theDevice->deviceId, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(theDevice->maxComputeUnits), &theDevice->maxComputeUnits, NULL);
      ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_MAX_COMPUTE_UNITS of device");

      status = clGetDeviceInfo(theDevice->deviceId, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(theDevice->localMemSize), &theDevice->localMemSize, NULL);
      ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_LOCAL_MEM_SIZE of device");
   }

   xfree(devices);
}

void Device::ListAllDevices(void)
{
   cl_uint  pp, dd;

   printf("List of available platforms and devices\n");

   for (pp=0; pp<ii_PlatformCount; pp++)
   {
      fprintf(stderr, "Platform %d is a %s %s, version %s\n", pp, ip_Platforms[pp].vendor,
              ip_Platforms[pp].name, ip_Platforms[pp].version);

      if (ip_Platforms[pp].deviceCount == 0)
         fprintf(stderr, "   No devices\n");

      for (dd=0; dd<ip_Platforms[pp].deviceCount; dd++)
         fprintf(stderr, "   Device %d is a %s %s\n", dd, ip_Platforms[pp].devices[dd].vendor,  
                 ip_Platforms[pp].devices[dd].name);
   }
}

void Device::Validate(void)
{
   if (!ib_HavePlatform) ii_PlatformId = 0;

   if (!ib_HaveDevice) ii_DeviceId = 0;

   if (ii_PlatformId >= ii_PlatformCount)
   {
      fprintf(stderr, "Platform %d does not exist.  Here is a list of platforms and devices:", ii_PlatformId);
      ListAllDevices();
      exit(0);
   }

   if (ip_Platforms[ii_PlatformId].deviceCount == 0)
   {
      fprintf(stderr, "Platform %d has no available devices.  Here is a list of platforms and devices:", ii_PlatformId);
      ListAllDevices();
      exit(0);
   }

   if (ii_DeviceId >= ip_Platforms[ii_PlatformId].deviceCount)
   {
      fprintf(stderr, "Platform %d/Device %d does not exist.  Here is a list of platforms and devices:",
              ii_PlatformId, ii_DeviceId);
      ListAllDevices();
      exit(0);
   }
}

bool Device::IsVowel(char ch)
{
  char c = toupper(ch);
 
  return (c == 'A' || c == 'E' || c == 'I' || c == 'O' || c == 'U');
}