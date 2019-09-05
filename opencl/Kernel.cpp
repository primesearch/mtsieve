/* Kernel.cpp -- (C) Mark Rodenkirch, February 2012

   This class provides the implementation for OpenCL kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "Kernel.h"
#include "ErrorChecker.h"
#include "../core/Clock.h"

Kernel::Kernel(Device *device, const char *kernelName, const char *kernelSource[], bool useFMA)
{
   cl_int  status;
   int     computeUnits, ii, count;
   size_t *sourceSize;
   size_t  len;
   char    buffer[16384];
   string  madEnable = "-cl-mad-enable";
   string  buildOptions = "";

   if (useFMA)
      buildOptions += madEnable;

   ii = 0;
   while (kernelSource[ii])
      ii++;

   count = ii;
   sourceSize = (size_t *) xmalloc(sizeof(size_t) * (ii+1));

   ii = 0;
   while (kernelSource[ii])
   {
      sourceSize[ii] = strlen(kernelSource[ii]);
      ii++;
   }
   sourceSize[ii] = 0;

   ip_Device = device;
   
   is_KernelName = kernelName;
   ii_ArgumentCount = 0;
 
    // Create an OpenCL command queue
#ifdef __APPLE__
   im_CommandQueue = clCreateCommandQueueWithPropertiesAPPLE(ip_Device->GetContext(), ip_Device->GetDeviceId(), 0, &status);
#else
   im_CommandQueue = clCreateCommandQueueWithProperties(ip_Device->GetContext(), ip_Device->GetDeviceId(), 0, &status);
#endif
   ErrorChecker::ExitIfError("clCreateCommandQueue", status);

   im_Program = clCreateProgramWithSource(ip_Device->GetContext(), count, kernelSource, sourceSize, &status);
   ErrorChecker::ExitIfError("clCreateProgramWithSource", status, "kernelSource: %s", kernelSource);

   // create a cl_rogram executable for all the devices specified
   status = clBuildProgram(im_Program, 1, ip_Device->GetDeviceIdPtr(), buildOptions.c_str(), NULL, NULL);

   if (status != CL_SUCCESS)
   {
      clGetProgramBuildInfo(im_Program, ip_Device->GetDeviceId(), CL_PROGRAM_BUILD_LOG, 16384, buffer, &len);
      ErrorChecker::ExitIfError("clBuildProgram", status, buffer);
   }

   // get a kernel object handle for a kernel with the given name
   im_Kernel = clCreateKernel(im_Program, kernelName, &status);
   ErrorChecker::ExitIfError("clCreateKernel", status, "kernelName: %s", kernelName);

   status = clGetKernelWorkGroupInfo(im_Kernel, ip_Device->GetDeviceId(), CL_KERNEL_WORK_GROUP_SIZE,
                                     sizeof(size_t), &ii_KernelWorkGroupSize, NULL);
   ErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_WORK_GROUP_SIZE", kernelName);
    
   computeUnits = ip_Device->GetMaxComputeUnits();

   ii_WorkGroupSize = computeUnits * ii_KernelWorkGroupSize;
   
   xfree(sourceSize);
}

Kernel::~Kernel(void)
{
   clReleaseKernel(im_Kernel);
   clReleaseProgram(im_Program);
   clReleaseCommandQueue(im_CommandQueue);
}

void Kernel::PrintStats(void)
{
   //VerboseOutput("workGroupSize = %d = %d * %d * %d (blocks * workGroupSizeMultiple * deviceComputeUnits)", 
   //              ii_WorkSize, ii_Blocks, (int) ii_Multiplier, ip_Device->GetMaxComputeUnits());
}

void Kernel::AddArgument(KernelArgument *kernelArgument)
{
   ip_Arguments[ii_ArgumentCount++] = kernelArgument;
}

void Kernel::Execute(uint32_t workSize)
{
   cl_int status;
   size_t globalWorkGroupSize[1];
   uint64_t startTime;

   startTime = Clock::GetCurrentMicrosecond();

   SetGPUInput();

   globalWorkGroupSize[0] = workSize;
   status = clEnqueueNDRangeKernel(im_CommandQueue, im_Kernel, 1, NULL, globalWorkGroupSize, NULL, 0, NULL, NULL);
   ErrorChecker::ExitIfError("clEnqueueNDRangeKernel", status, "kernelName: %s  globalworksize %u  localworksize %u", 
                             is_KernelName.c_str(), globalWorkGroupSize[0], ii_KernelWorkGroupSize);


   GetGPUOutput();

   ip_Device->AddGPUClockTime(Clock::GetCurrentMicrosecond() - startTime);
}

void Kernel::SetGPUInput(void)
{
   cl_int  status;
   int     aa;

   for (aa=0; aa<ii_ArgumentCount; aa++)
   {
      status = clSetKernelArg(im_Kernel, aa, sizeof(cl_mem), ip_Arguments[aa]->GetAddress());
      ErrorChecker::ExitIfError("clSetKernelArg", status, "kernelName: %s  index %d  argument: %s  size: %d",
                                is_KernelName.c_str(), aa, ip_Arguments[aa]->GetName().c_str(), ip_Arguments[aa]->GetSize());

      ip_Arguments[aa]->WriteToGPU(im_CommandQueue);
   }
}

void Kernel::GetGPUOutput(void)
{
   int     aa;

   for (aa=0; aa<ii_ArgumentCount; aa++)
      ip_Arguments[aa]->ReadFromGPU(im_CommandQueue);
}

