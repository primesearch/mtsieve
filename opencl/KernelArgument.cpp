/* KernelArgument.cpp -- (C) Mark Rodenkirch, February 2012

   This class provides the implementation for OpenCL kernel arguments.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#include "KernelArgument.h"
#include "ErrorChecker.h"

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int8_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint8_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int16_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint16_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int32_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint32_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int64_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint64_t *argument, int32_t count)
{
   is_ArgumentName = argumentName;
   ip_Device = theDevice;
   Initialize(direction, argument, sizeof(*argument), count);
}

KernelArgument::~KernelArgument(void)
{
  if (im_GPUBuffer) clReleaseMemObject(im_GPUBuffer);
}

void KernelArgument::Initialize(gpu_dir_t direction, void *argument, int32_t size, int32_t count)
{
   cl_int         status;
   cl_mem_flags   memFlags;

   if (direction & KA_HOST_TO_GPU) memFlags = CL_MEM_READ_ONLY;
   if (direction & KA_GPU_TO_HOST) memFlags = CL_MEM_WRITE_ONLY;
   if (direction & KA_BIDIRECTIONAL) memFlags = CL_MEM_READ_WRITE;

   ie_Direction = direction;
   ip_HostMemory = argument;
   ii_Count = count;
   ii_Bytes = ii_Count * size;

   im_GPUBuffer = clCreateBuffer(ip_Device->GetContext(), memFlags, ii_Bytes*2, NULL, &status);
   
   if (status != CL_SUCCESS)
      ErrorChecker::ExitIfError("clCreateBuffer", status, "bytes: %d", ii_Bytes);

   // printf("%u for %s (%u %u) %x %x\n", ii_Bytes, is_ArgumentName.c_str(), size, count, (void *) im_GPUBuffer, ip_HostMemory);

   ip_Device->IncrementGpuBytes(ii_Bytes);
}

void KernelArgument::WriteToGPU(cl_command_queue commandQueue)
{
   cl_int status;

   if (im_GPUBuffer && (ie_Direction & KA_HOST_TO_GPU))
   {
      status = clEnqueueWriteBuffer(commandQueue, im_GPUBuffer, CL_TRUE,
                                    0, ii_Bytes, ip_HostMemory, 0, NULL, NULL);

      ErrorChecker::ExitIfError("clEnqueueWriteBuffer", status, "argument: %s", is_ArgumentName.c_str());
   }
}

void KernelArgument::ReadFromGPU(cl_command_queue commandQueue)
{
   cl_int status;

   if (im_GPUBuffer && (ie_Direction & KA_GPU_TO_HOST))
   {
      // printf("reading %s %x %x %u, bytes %u\n", is_ArgumentName.c_str(), (void *) im_GPUBuffer, ip_HostMemory, ii_Count, ii_Bytes);
      
      status = clEnqueueReadBuffer(commandQueue, im_GPUBuffer, CL_TRUE,
                                   0, ii_Bytes, ip_HostMemory, 0, NULL, NULL);

      ErrorChecker::ExitIfError("clEnqueueReadBuffer", status, "argument: %s", is_ArgumentName.c_str());
   }
}

