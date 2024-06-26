/* OpenCLKernel.cpp -- (C) Mark Rodenkirch, May 2022
   This class provides the implementation for OpenCL kernels.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "OpenCLKernel.h"
#include "OpenCLErrorChecker.h"
#include "../core/App.h"
#include "../core/Clock.h"

OpenCLKernel::OpenCLKernel(GpuDevice *device, const char *kernelName, const char *kernelSource, const char *preOpenCLKernelSources[]) 
   : GpuKernel(device, kernelName, kernelSource, preOpenCLKernelSources)
{
   cl_int    status;
   int       computeUnits, ii;
   const char     *sources[3];
   size_t    sourcesSize[3];
   size_t    len;
   size_t    workGroupSizeMultiple;
   cl_ulong  deviceGlobalMemorySize;
   cl_ulong  deviceMaxMemAllocSize;
   cl_ulong  deviceLocalMemorySize;
   cl_ulong  deviceMaxConstantBufferSize;
   cl_ulong  localMemorySize;
   cl_ulong  privateMemorySize;
   char     *tempSource;
   std::string    madEnable = "-cl-mad-enable";
   std::string    buildOptions = "";

   tempSource = (char *) xmalloc(500000, 1, "tempSource");

   sources[0] = tempSource;
   sources[1] = kernelSource;
   sources[2] = NULL;
   
   snprintf(tempSource, 500000, "#define USE_OPENCL\n");
   
   if (preOpenCLKernelSources != NULL)
   {
      ii = 0;
      while (preOpenCLKernelSources[ii])
      {
         strcat(tempSource, preOpenCLKernelSources[ii]);
         strcat(tempSource, "\n");
         ii++;
      }
   }

   sourcesSize[0] = strlen(sources[0]);
   sourcesSize[1] = strlen(sources[1]);
   sourcesSize[2] = 0;

   ip_OpenCLDevice = (OpenCLDevice *) device;
   
   ip_KernelArguments = (ka_t *) xmalloc(MAX_KERNEL_ARGUMENTS, sizeof(ka_t), "kernelArguments");
   is_KernelName = kernelName;
   ii_ArgumentCount = 0;
 
    // Create an OpenCL command queue
#ifdef __APPLE__
   const cl_queue_properties_APPLE queue_properties = 0;

   im_CommandQueue = clCreateCommandQueueWithPropertiesAPPLE(ip_OpenCLDevice->GetContext(), ip_OpenCLDevice->GetDeviceId(), &queue_properties, &status);
#else
   im_CommandQueue = clCreateCommandQueueWithProperties(ip_OpenCLDevice->GetContext(), ip_OpenCLDevice->GetDeviceId(), 0, &status);
#endif
   OpenCLErrorChecker::ExitIfError("clCreateCommandQueue", status);

   im_Program = clCreateProgramWithSource(ip_OpenCLDevice->GetContext(), 2, sources, sourcesSize, &status);
   OpenCLErrorChecker::ExitIfError("clCreateProgramWithSource", status, "kernelName: %s", kernelName);

   xfree(tempSource);

   // create a cl_rogram executable for all the devices specified
   status = clBuildProgram(im_Program, 1, ip_OpenCLDevice->GetDeviceIdPtr(), buildOptions.c_str(), NULL, NULL);

   if (status != CL_SUCCESS)
   {
      char      buffer[100000];
      
      clGetProgramBuildInfo(im_Program, ip_OpenCLDevice->GetDeviceId(), CL_PROGRAM_BUILD_LOG, (uint32_t) sizeof(buffer), buffer, &len);
      OpenCLErrorChecker::ExitIfError("clBuildProgram", status, "%s", buffer);
   }

   status = clGetDeviceInfo(ip_OpenCLDevice->GetDeviceId(), CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(deviceGlobalMemorySize), &deviceGlobalMemorySize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_GLOBAL_MEM_SIZE of device");
   
   status = clGetDeviceInfo(ip_OpenCLDevice->GetDeviceId(), CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(deviceMaxMemAllocSize), &deviceMaxMemAllocSize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_MAX_MEM_ALLOC_SIZE of device");
   
   status = clGetDeviceInfo(ip_OpenCLDevice->GetDeviceId(), CL_DEVICE_LOCAL_MEM_SIZE, sizeof(deviceLocalMemorySize), &deviceLocalMemorySize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_LOCAL_MEM_SIZE of device");
      
   status = clGetDeviceInfo(ip_OpenCLDevice->GetDeviceId(), CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(deviceMaxConstantBufferSize), &deviceMaxConstantBufferSize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetDeviceInfo", status, "Unable to get CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE of device");
   
   // get a kernel object handle for a kernel with the given name
   im_OpenCLKernel = clCreateKernel(im_Program, kernelName, &status);
   OpenCLErrorChecker::ExitIfError("clCreateOpenCLKernel", status, "kernelName: %s", kernelName);

   status = clGetKernelWorkGroupInfo(im_OpenCLKernel, ip_OpenCLDevice->GetDeviceId(), CL_KERNEL_WORK_GROUP_SIZE,
                                     sizeof(ii_KernelWorkGroupSize), &ii_KernelWorkGroupSize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_WORK_GROUP_SIZE", kernelName);
    
   status = clGetKernelWorkGroupInfo(im_OpenCLKernel, ip_OpenCLDevice->GetDeviceId(), CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,
                                     sizeof(workGroupSizeMultiple), &workGroupSizeMultiple, NULL);
   OpenCLErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE", kernelName);
   
   status = clGetKernelWorkGroupInfo(im_OpenCLKernel, ip_OpenCLDevice->GetDeviceId(), CL_KERNEL_LOCAL_MEM_SIZE,
                                     sizeof(privateMemorySize), &localMemorySize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_LOCAL_MEM_SIZE", kernelName);
   
   status = clGetKernelWorkGroupInfo(im_OpenCLKernel, ip_OpenCLDevice->GetDeviceId(), CL_KERNEL_PRIVATE_MEM_SIZE,
                                     sizeof(privateMemorySize), &privateMemorySize, NULL);
   OpenCLErrorChecker::ExitIfError("clGetKernelWorkGroupInfo", status, "kernelName: %s  argument CL_KERNEL_PRIVATE_MEM_SIZE", kernelName);
   
   il_DeviceGlobalMemorySize = deviceGlobalMemorySize;
   il_DeviceMaxMemAllocSize = deviceMaxMemAllocSize;
   ii_DeviceLocalMemorySize = deviceLocalMemorySize;
   ii_DeviceMaxConstantBufferSize = deviceMaxConstantBufferSize;
   ii_LocalMemorySize = localMemorySize;
   ii_PrivateMemorySize = privateMemorySize;
   ii_WorkGroupSizeMultiple = workGroupSizeMultiple;
      
   computeUnits = ip_OpenCLDevice->GetMaxComputeUnits();

   ii_WorkGroupSize = computeUnits * ii_KernelWorkGroupSize;
}

OpenCLKernel::~OpenCLKernel(void)
{
   ka_t      *ka;
   uint32_t   aa;
   
   clReleaseKernel(im_OpenCLKernel);
   clReleaseProgram(im_Program);
   clReleaseCommandQueue(im_CommandQueue);

   for (aa=0; aa<ii_ArgumentCount; aa++)
   {
      ka = &ip_KernelArguments[aa];
      
      if (ka->mustFreeGpuBuffer)
         clReleaseMemObject(ka->gpuBuffer);
      
      if (ka->mustFreeCpuBuffer)
         xfree(ka->cpuBuffer);
   }
   
   xfree(ip_KernelArguments);
}

void  *OpenCLKernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, NULL, CL_MEM_READ_ONLY);
}

void  *OpenCLKernel::AddCpuArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory)
{
   return AddArgument(name, size, count, cpuMemory, CL_MEM_READ_ONLY);
}

void  *OpenCLKernel::AddGpuArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, NULL, CL_MEM_WRITE_ONLY);
}

void  *OpenCLKernel::AddSharedArgument(const char *name, uint32_t size, uint32_t count)
{
   return AddArgument(name, size, count, NULL, CL_MEM_READ_WRITE);
}

void  *OpenCLKernel::AddArgument(const char *name, GpuKernel *other)
{
   OpenCLKernel *kernel = (OpenCLKernel *) other;
   
   ka_t      *ka = &ip_KernelArguments[ii_ArgumentCount];
   ka_t      *otherKa = kernel->GetKernelArgument(name);
   
   ii_ArgumentCount++;
   
   strcpy(ka->name, otherKa->name);
   
   ka->size = otherKa->size;
   ka->count = otherKa->count;
   ka->bytes = otherKa->bytes;
   ka->memFlags = otherKa->memFlags;
   
   memcpy(&ka->gpuBuffer, &otherKa->gpuBuffer, sizeof(otherKa->gpuBuffer));
   ka->cpuBuffer = otherKa->cpuBuffer;
   
   ka->mustFreeGpuBuffer = false;
   ka->mustFreeCpuBuffer = false;
  
   return ka->cpuBuffer;
}

ka_t  *OpenCLKernel::GetKernelArgument(const char *name)
{
   ka_t      *ka;
   
   for (uint32_t idx=0; idx<ii_ArgumentCount; idx++)
   {
      ka = &ip_KernelArguments[idx];
      
      if (!strcmp(ka->name, name))
         return ka;
   }
   
   FatalError("Could not find KernelArgument with the name %s", name);
   return NULL;
}

void  *OpenCLKernel::AddArgument(const char *name, uint32_t size, uint32_t count, void *cpuMemory, cl_mem_flags memFlags)
{
   ka_t      *ka = &ip_KernelArguments[ii_ArgumentCount];
   cl_int     status;
   
   ii_ArgumentCount++;

   strcpy(ka->name, name);
   ka->size = size;
   ka->count = count;
   ka->memFlags = memFlags;
   ka->mustFreeGpuBuffer = true;
   
   if (cpuMemory == NULL)
   {
      ka->cpuBuffer = xmalloc(count + 1, size, name);
      ka->mustFreeCpuBuffer = true;
   }
   else
   {
      ka->cpuBuffer = cpuMemory;
      ka->mustFreeCpuBuffer = false;
   }
   
   ka->bytes = size * count;

   ka->gpuBuffer = clCreateBuffer(ip_OpenCLDevice->GetContext(), memFlags, ka->bytes, NULL, &status);
   
   ip_OpenCLDevice->IncrementGpuBytes(ka->bytes);
   
   if (status != CL_SUCCESS)
      OpenCLErrorChecker::ExitIfError("clCreateBuffer", status, "bytes: %d", ka->bytes);
   
   return ka->cpuBuffer;
}

void OpenCLKernel::PrintStatistics(uint64_t bytesPerWorkGroup)
{
   if (ip_OpenCLDevice->IsPrintDetails())
   {
      App   *theApp = get_app();

      uint64_t privateBytes = bytesPerWorkGroup;
      
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_MAX_COMPUTE_UNITS = %u", ip_OpenCLDevice->GetMaxComputeUnits());
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_GLOBAL_MEM_SIZE = %" PRIu64"", il_DeviceGlobalMemorySize);
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_MAX_MEM_ALLOC_SIZE = %" PRIu64"", il_DeviceMaxMemAllocSize);
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_LOCAL_MEM_SIZE = %u", ii_DeviceLocalMemorySize);
      theApp->WriteToConsole(COT_OTHER, "CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE = %u", ii_DeviceMaxConstantBufferSize);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_WORK_GROUP_SIZE = %u", (uint32_t) ii_KernelWorkGroupSize);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE = %u", ii_WorkGroupSizeMultiple);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_LOCAL_MEM_SIZE = %u", ii_LocalMemorySize);
      theApp->WriteToConsole(COT_OTHER, "CL_KERNEL_PRIVATE_MEM_SIZE = %u", ii_PrivateMemorySize);
      
      theApp->WriteToConsole(COT_OTHER, "GPU global bytes allocated = %" PRIu64"", ip_OpenCLDevice->GetGpuBytes());
      
      if (privateBytes > 0)
         theApp->WriteToConsole(COT_OTHER, "GPU private bytes allocated = %" PRIu64"", bytesPerWorkGroup * ii_KernelWorkGroupSize);
   }
}

void OpenCLKernel::Execute(uint32_t workSize)
{
   cl_int status;
   size_t globalWorkGroupSize[1];
   uint64_t startTime;

   startTime = Clock::GetCurrentMicrosecond();

   SetGPUInput();

   globalWorkGroupSize[0] = workSize;

   status = clEnqueueNDRangeKernel(im_CommandQueue, im_OpenCLKernel, 1, NULL, globalWorkGroupSize, NULL, 0, NULL, NULL);

   OpenCLErrorChecker::ExitIfError("clEnqueueNDRangeOpenCLKernel", status, "kernelName: %s  globalworksize %u  privatememsize %u", 
                             is_KernelName.c_str(), (uint32_t) globalWorkGroupSize[0], (uint32_t) ii_PrivateMemorySize);

   GetGPUOutput();
   
   status =  clFinish(im_CommandQueue);

   OpenCLErrorChecker::ExitIfError("clFinish", status, "kernelName: %s",  is_KernelName.c_str());
   
   ip_OpenCLDevice->AddGpuMicroseconds(Clock::GetCurrentMicrosecond() - startTime);
}

void OpenCLKernel::SetGPUInput(void)
{
   ka_t      *ka;
   cl_int     status;
   uint32_t   aa;

   for (aa=0; aa<ii_ArgumentCount; aa++)
   {
      ka = &ip_KernelArguments[aa];
      status = clSetKernelArg(im_OpenCLKernel, aa, sizeof(cl_mem), &ka->gpuBuffer);
      
      OpenCLErrorChecker::ExitIfError("clSetOpenCLKernelArg", status, "kernelName: %s  index %d  argument: %s  size: %d",
                                is_KernelName.c_str(), aa, ka->name, ka->bytes);

      if (ka->memFlags == CL_MEM_WRITE_ONLY)
         continue;
      
      status = clEnqueueWriteBuffer(im_CommandQueue, ka->gpuBuffer, CL_TRUE, 0, ka->bytes, ka->cpuBuffer, 0, NULL, NULL);

      OpenCLErrorChecker::ExitIfError("clEnqueueWriteBuffer", status, "argument: %s", ka->name);
   }
}

void OpenCLKernel::GetGPUOutput(void)
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

      OpenCLErrorChecker::ExitIfError("clEnqueueReadBuffer", status, "argument: %s", ka->name);
   }
}
