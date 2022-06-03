/* GpuKernel.cpp -- (C) Mark Rodenkirch, May 2022

   This class provides the generic implementation for GPU kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "GpuKernel.h"
#include "../core/App.h"
#include "../core/Clock.h"

GpuKernel::GpuKernel(GpuDevice *gpuDevice, const char *kernelName, const char *kernelSource, const char *preKernelSources[])
{
   ip_GpuDevice = gpuDevice;
}

GpuKernel::~GpuKernel(void)
{
}
