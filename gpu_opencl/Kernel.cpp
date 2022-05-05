/* Kernel.cpp -- (C) Mark Rodenkirch, February 2012

   This class provides the implementation for OpenCL kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "Kernel.h"
#include "ErrorChecker.h"
#include "../core/App.h"
#include "../core/Clock.h"

Kernel::Kernel(Device *device, const char *kernelName, const char *kernelSource, const char *preKernelSources[])
{
   cl_int    status;
   int       computeUnits, ii;
   const char     *sources[3];
   size_t    sourcesSize[3];
   size_t    len;
   size_t    workGroupSizeMultiple;
   cl_ulong  deviceGlobalMemorySize;
   cl_ulong  deviceLocalMemorySize;
   cl_ulong  localMemorySize;
   cl_ulong  privateMemorySize;
   char      buffer[16384];
   char     *tempSource;
   string    madEnable = "-cl-mad-enable";
   string    buildOptions = "";

   tempSource = (char *) xmalloc(50000);

   sources[0] = tempSource;
   sources[1] = kernelSource;
   sources[2] = NULL;
   
   sprintf(tempSource, "#define USE_OPENCL\n");
   
   if (preKernelSources != NULL)
   {
      ii = 0;
      while (preKernelSources[ii])
      {
         strcat(tempSource, preKernelSources[ii]);
         strcat(tempSource, "\n");
         ii++;
      }
   }

   sourcesSize[0] = strlen(sources[0]);
   sourcesSize[1] = strlen(sources[1]);
   sourcesSize[2] = 0;

   ip_Device = device;
   
   ip_KernelArguments = (ka_t *) xmalloc(sizeof(ka_t) * MAX_KERNEL_ARGUMENTS);
   is_KernelName = kernelName;
   ii_ArgumentCount = 0;
 
    // Create an OpenCL command queue
#ifdef __APPLE__
   const cl_queue_properties_APPLE queue_properties = 0;

   im_CommandQueue = clCreateCommandQueueWithPropertiesAPPLE(ip_Device->GetContext(), ip_Device->GetDeviceId(), &queue_properties, &status);
#else
   im_CommandQueue = clCreateCommandQueueWithProperties(ip_Device->GetContext(), ip_Device->GetDeviceId(), 0, &status);
#endif
   ErrorChecker::ExitIfError("clCreateCommandQueue", status);

   im_Program = clCreateProgramWithSource(ip_Device->GetContext(), 2, sources, sourcesSize, &status);
   ErrorChecker::ExitIfError("clCreateProgramWithSource", status, "kernelName: %s", kernelName);

   xfree(tempSource);

   // create a cl_rogram executable for all the devices specified
   status = clBuildProgram(im_Program, 1, ip_Device->GetDeviceIdPtr(), buildOptions.c_str(), NULL, NULL);

   if (status != CL_SUCCESS)
   {
      clGetProgramBuildInfo(im_Program, ip_Device->GetDeviceId(), CL_PROGRAM_BUILD_LOG, (uint32_t) sizeof(buffer), buffer, &len);
      ErrorChecker::ExitIfError("clBuildProgram", status, "%s", buffer);
   }


   status = clGetDeviceInfo(ip_Device->GetDeviceId(), CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(deviceGlobalMemorySize), &deviceGlobalMemorySize, NULL);
   ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_GLOBAL_MEM_SIZE of device");
   
   status = clGetDeviceInfo(ip_Device->GetDeviceId(), CL_DEVICE_LOCAL_MEM_SIZE, sizeof(deviceLocalMemorySize), &deviceLocalMemorySize, NULL);
   ErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_LOCAL_MEM_SIZE of device");
      
   // get a kernel object handle for a kernel with the given name
   im_Kernel = clCreateKernel(im_Program, kernelName, &status);
   ErrorChecker::ExitIfError("clCreateKernel", status, "kernelName: %s", kernelName);

   status = clGetKernelWorkGroupInfo(im_Kernel, ip_Device->GetDeviceId(), CL_KERNEL_WORK_GROUP_SIZE,
                                     sizeof(ii_KernelWorkGroupSize), &ii_KernelWorkGroupSize, NULL);
   ErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_WORK_GROUP_SIZE", kernelName);
    
   status = clGetKernelWorkGroupInfo(im_Kernel, ip_Device->GetDeviceId(), CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                     sizeof(workGroupSizeMultiple), &workGroupSizeMultiple, NULL);
   ErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE", kernelName);
   
   status = clGetKernelWorkGroupInfo(im_Kernel, ip_Device->GetDeviceId(), CL_KERNEL_LOCAL_MEM_SIZE,
                                     sizeof(privateMemorySize), &localMemorySize, NULL);
   ErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_LOCAL_MEM_SIZE", kernelName);
   
   status = clGetKernelWorkGroupInfo(im_Kernel, ip_Device->GetDeviceId(), CL_KERNEL_PRIVATE_MEM_SIZE,
                                     sizeof(privateMemorySize), &privateMemorySize, NULL);
   ErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_PRIVATE_MEM_SIZE", kernelName);
   
   ii_DeviceGlobalMemorySize = deviceGlobalMemorySize;
   ii_DeviceLocalMemorySize = deviceLocalMemorySize;
   ii_LocalMemorySize = localMemorySize;
   ii_PrivateMemorySize = privateMemorySize;
   ii_WorkGroupSizeMultiple = workGroupSizeMultiple;
      
   computeUnits = ip_Device->GetMaxComputeUnits();

   ii_WorkGroupSize = computeUnits * ii_KernelWorkGroupSize;
}

Kernel::~Kernel(void)
{
   ka_t      *ka;
   uint32_t   aa;
   
   clReleaseKernel(im_Kernel);
   clReleaseProgram(im_Program);
   clReleaseCommandQueue(im_CommandQueue);

   for (aa=0; aa<ii_ArgumentCount; aa++)
   {
      ka = &ip_KernelArguments[aa];
      clReleaseMemObject(ka->gpuBuffer);
      xfree(ka->cpuBuffer);
   }
   
   xfree(ip_KernelArguments);
}

void  *Kernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, CL_MEM_READ_ONLY);
}

void  *Kernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory)
{
   void *ptr = AddArgument(name, size, count, CL_MEM_READ_ONLY);
   
   memcpy(ptr, cpuMemory, size * count);
   
   return ptr;
}

void  *Kernel::AddGpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, CL_MEM_WRITE_ONLY);
}

void  *Kernel::AddSharedArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, CL_MEM_READ_WRITE);
}

void  *Kernel::AddArgument(const char *name, uint32_t size, uint32_t count, cl_mem_flags memFlags)
{
   ka_t      *ka = &ip_KernelArguments[ii_ArgumentCount];
   cl_int     status;
   
   ii_ArgumentCount++;

   strcpy(ka->name, name);
   ka->size = size;
   ka->count = count;
   ka->memFlags = memFlags;
   ka->cpuBuffer = xmalloc(size * (count + 1));
   ka->bytes = size * count;

   // Don't know why I get "out of resources" if I don't make this slightly larger
   ka->gpuBuffer = clCreateBuffer(ip_Device->GetContext(), memFlags, ka->bytes, NULL, &status);
   
   if (status != CL_SUCCESS)
      ErrorChecker::ExitIfError("clCreateBuffer", status, "bytes: %d", ka->bytes);
   
   return ka->cpuBuffer;
}

void Kernel::PrintStatistics(uint64_t bytesPerWorkGroup)
{
   if (ip_Device->IsPrintDetails())
   {
      App   *theApp = get_app();

      uint64_t privateBytes = bytesPerWorkGroup;
      
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_MAX_COMPUTE_UNITS = %u", ip_Device->GetMaxComputeUnits());
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_GLOBAL_MEM_SIZE = %u", ii_DeviceGlobalMemorySize);
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_LOCAL_MEM_SIZE = %u", ii_DeviceLocalMemorySize);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_WORK_GROUP_SIZE = %u", (uint32_t) ii_KernelWorkGroupSize);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE = %u", ii_WorkGroupSizeMultiple);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_LOCAL_MEM_SIZE = %u", ii_LocalMemorySize);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_PRIVATE_MEM_SIZE = %u", ii_PrivateMemorySize);
      
      theApp->WriteToConsole(COT_OTHER, "GPU global bytes allocated = %" PRIu64"", ip_Device->GetGpuBytes());
      
      if (privateBytes > 0)
         theApp->WriteToConsole(COT_OTHER, "GPU private bytes allocated = %" PRIu64"", bytesPerWorkGroup * ii_KernelWorkGroupSize);
   }
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
                             is_KernelName.c_str(), (uint32_t) globalWorkGroupSize[0], (uint32_t) ii_KernelWorkGroupSize);

   ErrorChecker::ExitIfError("clFinish", status, "kernelName: %s",  is_KernelName.c_str());
                             
   GetGPUOutput();
   
   status =  clFinish(im_CommandQueue);

   ip_Device->AddGpuMicroseconds(Clock::GetCurrentMicrosecond() - startTime);
}

void Kernel::SetGPUInput(void)
{
   ka_t      *ka;
   cl_int     status;
   uint32_t   aa;

   for (aa=0; aa<ii_ArgumentCount; aa++)
   {
      ka = &ip_KernelArguments[aa];
      status = clSetKernelArg(im_Kernel, aa, sizeof(cl_mem), &ka->gpuBuffer);
      
      ErrorChecker::ExitIfError("clSetKernelArg", status, "kernelName: %s  index %d  argument: %s  size: %d",
                                is_KernelName.c_str(), aa, ka->name, ka->bytes);

      if (ka->memFlags == CL_MEM_WRITE_ONLY)
         continue;
      
      status = clEnqueueWriteBuffer(im_CommandQueue, ka->gpuBuffer, CL_TRUE, 0, ka->bytes, ka->cpuBuffer, 0, NULL, NULL);

      ErrorChecker::ExitIfError("clEnqueueWriteBuffer", status, "argument: %s", ka->name);
   }
}

void Kernel::GetGPUOutput(void)
{
   ka_t      *ka;
   cl_int     status;
   uint32_t   aa;

   for (aa=0; aa<ii_ArgumentCount; aa++)
   {
      ka = &ip_KernelArguments[aa];
            
      if (ka->memFlags == CL_MEM_READ_ONLY)
         continue;
           
      status = clEnqueueReadBuffer(im_CommandQueue, ka->gpuBuffer, CL_TRUE, 0, ka->bytes, ka->cpuBuffer, 0, NULL, NULL);

      ErrorChecker::ExitIfError("clEnqueueReadBuffer", status, "argument: %s", ka->name);
   }
}

