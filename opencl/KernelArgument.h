/* KernelArgument.h -- (C) Mark Rodenkirch, February 2012

   This class provides the interface for OpenCL kernel arguments.
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
*/

#ifndef _KERNELARGUMENT_H
#define _KERNELARGUMENT_H

#include "../core/main.h"
#include "Device.h"

using namespace std;

typedef enum { KA_HOST_TO_GPU = 0x01, KA_GPU_TO_HOST = 0x02, KA_BIDIRECTIONAL = 0x03 } gpu_dir_t;

class KernelArgument
{
public:
   // Various constructors, each one handles a different datatype
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int8_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint8_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int16_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint16_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int32_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint32_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, int64_t *argument, int32_t count);
   KernelArgument(Device *theDevice, const char *argumentName, gpu_dir_t direction, uint64_t *argument, int32_t count);

   ~KernelArgument(void);

   // Copy from host memory to GPU memory
   void              WriteToGPU(cl_command_queue commandQueue);

   // Copy from GPU memory to host memory
   void              ReadFromGPU(cl_command_queue commandQueue);

   string            GetName(void)    { return is_ArgumentName; };
   cl_mem           *GetAddress(void) { return &im_GPUBuffer; };
   int               GetSize(void)    { return ii_Bytes; };

   void              SetHostMemory(void *ptr) { ip_HostMemory = ptr; };
   
private:
   Device           *ip_Device;
   string            is_ArgumentName;
   gpu_dir_t         ie_Direction;
   void             *ip_HostMemory;
   int               ii_Bytes;
   int               ii_Count;
   
   cl_mem            im_GPUBuffer;

   void              Initialize(gpu_dir_t direction, void *argument, int32_t size, int32_t count);
   void              AllocateBuffer(void);
};

#endif

