/* GpuDevice.cpp -- (C) Mark Rodenkirch, May 2022

   This class provides the generic implementation for GPU devices.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 */

#include "../core/main.h"

#include "GpuDevice.h"

GpuDevice::GpuDevice(void)
{
   ip_GpuBytes = new SharedMemoryItem("gpu_bytes");
   ip_GpuMicroseconds = new SharedMemoryItem("gpu_microseconds");
}

GpuDevice::~GpuDevice(void)
{
}

void  GpuDevice::CleanUp(void)
{
   delete ip_GpuBytes;
   delete ip_GpuMicroseconds;
}

void  GpuDevice::ParentHelp(void)
{
}

void  GpuDevice::ParentAddCommandLineOptions(std::string &shortOpts, struct option *longOpts)
{
}

// Returns:
//    0 if the option is OK
//   -1 if the argument is invalid
//   -2 if the argument is out of range
//   99 if the argument is not supported by this module
parse_t GpuDevice::ParentParseOption(int opt, char *arg, const char *source)
{
   return P_UNSUPPORTED;
}

void GpuDevice::ParentValidate(void)
{
}
